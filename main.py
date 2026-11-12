"""
RestrictionMap - find where restriction enzymes cut a DNA sequence.

Functional skeleton: every feature works, but pages use plain default Drafter
components with no styling yet.

File layout:
    1. Dataclasses (RestrictionEnzyme, SavedAnalysis, State)
    2. Data: global constants, the enzyme catalog and sample sequences
    3. Analysis helper functions (no Drafter) + their unit tests
    4. Routes (Drafter pages)
    5. Unit tests for every route
    6. start_server

Conventions used everywhere in this file:
    - Python strings are 0-based, but every number SHOWN to the user is 1-based
      (the first base of the sequence is "base 1").
    - A "cut position" is the number of bases to the LEFT of the cut on the top
      strand. That same number is also the 1-based position of the last base
      before the cut, so a cut position of 4 means "cuts after base 4".
    - DNA is treated as linear (not a circular plasmid).
"""

from dataclasses import dataclass
from drafter import *


# ---------------------------------------------------------------------------
# 1. Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class RestrictionEnzyme:
    """
    One restriction enzyme from the catalog.

    name: e.g. "EcoRI"
    recognition_sequence: the DNA site it recognizes, e.g. "GAATTC"
    cut_offset: how many bases into the recognition sequence the top strand is
        cut. EcoRI cuts G^AATTC, so its offset is 1.
    description: a short note shown on the enzyme selection page
    """
    name: str
    recognition_sequence: str
    cut_offset: int
    description: str


@dataclass
class SavedAnalysis:
    """A finished analysis the user chose to keep on the Saved Analyses page."""
    sequence_name: str
    dna_sequence: str
    enzyme_names: list[str]
    cut_positions: list[int]
    fragment_lengths: list[int]


@dataclass
class State:
    """
    Everything the site remembers between pages.

    sequence_name: set by validate_sequence_page, load_sample, reopen_analysis
    dna_sequence: set by validate_sequence_page, load_sample, reopen_analysis
    selected_enzymes: the RestrictionEnzyme instances the user ticked;
        set by analysis_page and reopen_analysis, cleared by clear_results
    cut_positions: set by analysis_page and reopen_analysis
    fragment_lengths: set by analysis_page and reopen_analysis
    saved_analyses: added to by save_analysis, removed from by delete_analysis
    """
    sequence_name: str
    dna_sequence: str
    selected_enzymes: list[RestrictionEnzyme]
    cut_positions: list[int]
    fragment_lengths: list[int]
    saved_analyses: list[SavedAnalysis]


def empty_state() -> State:
    """Return a brand-new State with nothing entered yet."""
    return State("", "", [], [], [], [])


def clear_results(state: State) -> State:
    """
    Forget the last analysis (used whenever the DNA sequence changes).

    Changes the given State and also returns it, so it is easy to test.
    """
    state.selected_enzymes = []
    state.cut_positions = []
    state.fragment_lengths = []
    return state


# ---------------------------------------------------------------------------
# 2. Data
# ---------------------------------------------------------------------------

# Global constants (never changed while the site runs).
DNA_BASES = "ATGC"
WHITESPACE = " \t\n\r"
MAP_WIDTH = 60
ENZYME_TABLE_HEADER = ["Enzyme", "Site", "Sites found", "Site starts at base",
                       "Cuts after base"]
FRAGMENT_TABLE_HEADER = ["Fragment", "From base", "To base", "Length (bp)"]
COMPARE_TABLE_HEADER = ["Enzyme", "Site", "Cuts", "Fragments",
                        "Smallest (bp)", "Largest (bp)"]


# Every enzyme here has a PALINDROMIC site: the site reads the same on the
# bottom strand (read 5'->3'). That is why searching only the top strand finds
# every site. A non-palindromic enzyme (e.g. BsaI) would also need a search of
# the reverse complement, which this project does not do.
#
# TODO: double-check every site and cut offset against NEB's published data.
def get_enzyme_catalog() -> list[RestrictionEnzyme]:
    """Return every enzyme the user can pick from, in display order."""
    return [
        RestrictionEnzyme("EcoRI", "GAATTC", 1,
                          "G^AATTC. Leaves a 4-base 5' overhang (sticky)."),
        RestrictionEnzyme("BamHI", "GGATCC", 1,
                          "G^GATCC. Leaves a 4-base 5' overhang (sticky)."),
        RestrictionEnzyme("HindIII", "AAGCTT", 1,
                          "A^AGCTT. Leaves a 4-base 5' overhang (sticky)."),
        RestrictionEnzyme("XbaI", "TCTAGA", 1,
                          "T^CTAGA. Leaves a 4-base 5' overhang (sticky)."),
        RestrictionEnzyme("SalI", "GTCGAC", 1,
                          "G^TCGAC. Leaves a 4-base 5' overhang (sticky)."),
        RestrictionEnzyme("XhoI", "CTCGAG", 1,
                          "C^TCGAG. Leaves a 4-base 5' overhang (sticky)."),
        RestrictionEnzyme("PstI", "CTGCAG", 5,
                          "CTGCA^G. Leaves a 4-base 3' overhang."),
        RestrictionEnzyme("KpnI", "GGTACC", 5,
                          "GGTAC^C. Leaves a 4-base 3' overhang."),
        RestrictionEnzyme("SacI", "GAGCTC", 5,
                          "GAGCT^C. Leaves a 4-base 3' overhang."),
        RestrictionEnzyme("SphI", "GCATGC", 5,
                          "GCATG^C. Leaves a 4-base 3' overhang."),
        RestrictionEnzyme("SmaI", "CCCGGG", 3,
                          "CCC^GGG. Cuts straight across (blunt end)."),
        RestrictionEnzyme("EcoRV", "GATATC", 3,
                          "GAT^ATC. Cuts straight across (blunt end)."),
        RestrictionEnzyme("NotI", "GCGGCCGC", 2,
                          "GC^GGCCGC. 8-base site, so it cuts rarely."),
    ]


def get_sample_sequences() -> dict[str, str]:
    """Return the Sequence Library samples as {sample name: DNA sequence}."""
    return {
        # Hand-checked example from the design doc: EcoRI sites start at
        # bases 4 and 12, so EcoRI cuts after bases 4 and 12.
        "Design doc example": "ATGGAATTCCGGAATTC",
        # pUC19's multiple cloning site: 10 enzyme sites in a row.
        # TODO: confirm against the pUC19 map on Addgene or NEB.
        "pUC19 polylinker": "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGG"
                            + "CATGCAAGCTT",
        # Made-up sequence with no sites for any catalog enzyme.
        "No-sites example": "ATATATATTTTTAAAAACCCCATATGGGTTTAAA",
    }


# ---------------------------------------------------------------------------
# 3. Analysis helper functions (no Drafter, so they can be tested on their own)
# ---------------------------------------------------------------------------

def get_enzyme_names_of(enzymes: list[RestrictionEnzyme]) -> list[str]:
    """Return the name of each enzyme in the list, in the same order."""
    names = []
    for enzyme in enzymes:
        names.append(enzyme.name)
    return names


def get_enzyme_names() -> list[str]:
    """Return the name of every enzyme in the catalog."""
    return get_enzyme_names_of(get_enzyme_catalog())


def get_enzymes_by_names(names: list[str]) -> list[RestrictionEnzyme]:
    """Return the catalog enzymes whose names are in `names` (catalog order)."""
    chosen = []
    for enzyme in get_enzyme_catalog():
        if enzyme.name in names:
            chosen.append(enzyme)
    return chosen


def join_with_commas(words: list[str]) -> str:
    """Join words into one string, separated by a comma and a space."""
    text = ""
    for i in range(len(words)):
        if i > 0:
            text = text + ", "
        text = text + words[i]
    return text


def numbers_to_text(numbers: list[int]) -> str:
    """Turn a list of numbers into one comma-separated string."""
    words = []
    for number in numbers:
        words.append(str(number))
    return join_with_commas(words)


def to_one_based_text(zero_based_positions: list[int]) -> str:
    """Turn 0-based string indexes into comma-separated 1-based positions."""
    one_based = []
    for position in zero_based_positions:
        one_based.append(position + 1)
    return numbers_to_text(one_based)


def describe_sequence(sequence_name: str, dna_sequence: str) -> str:
    """Return a short label like 'Demo (17 bp)'."""
    return sequence_name + " (" + str(len(dna_sequence)) + " bp)"


def describe_enzyme(enzyme: RestrictionEnzyme) -> str:
    """Return the label shown next to an enzyme's checkbox."""
    return (enzyme.name + " (" + enzyme.recognition_sequence + "): "
            + enzyme.description)


def clean_sequence(raw_text: str) -> str:
    """Remove spaces, tabs and line breaks, and convert to uppercase."""
    cleaned = ""
    for character in raw_text:
        if character not in WHITESPACE:
            cleaned = cleaned + character.upper()
    return cleaned


def find_invalid_characters(sequence: str) -> list[str]:
    """Return each character that is not A, T, G or C (each listed once)."""
    invalid = []
    for character in sequence:
        if character not in DNA_BASES and character not in invalid:
            invalid.append(character)
    return invalid


def is_valid_dna(sequence: str) -> bool:
    """A valid sequence is non-empty and contains only A, T, G and C."""
    return len(sequence) > 0 and len(find_invalid_characters(sequence)) == 0


def find_site_starts(dna: str, recognition_sequence: str) -> list[int]:
    """
    Return the 0-based start index of every match, overlapping ones included.

    The window slides one base at a time, so overlapping sites in repetitive
    DNA are all found.
    """
    starts = []
    site_length = len(recognition_sequence)
    for i in range(len(dna) - site_length + 1):
        if dna[i:i + site_length] == recognition_sequence:
            starts.append(i)
    return starts


def find_cut_positions(dna: str, enzyme: RestrictionEnzyme) -> list[int]:
    """Return where one enzyme cuts: each site start plus its cut offset."""
    cuts = []
    for start in find_site_starts(dna, enzyme.recognition_sequence):
        cuts.append(start + enzyme.cut_offset)
    return cuts


def combine_cut_positions(dna: str,
                          enzymes: list[RestrictionEnzyme]) -> list[int]:
    """
    Return the cuts from ALL the given enzymes together, sorted, no repeats.

    This models a multi-enzyme digest: every enzyme cuts in the same tube.
    """
    all_cuts = []
    for enzyme in enzymes:
        for cut in find_cut_positions(dna, enzyme):
            if cut not in all_cuts:
                all_cuts.append(cut)
    return sorted(all_cuts)


def calculate_fragment_lengths(sequence_length: int,
                               cut_positions: list[int]) -> list[int]:
    """
    Return fragment lengths for LINEAR DNA, from left to right.

    cut_positions must be sorted. Fragments run from the start to the first
    cut, between each pair of cuts, and from the last cut to the end, so
    n cuts always give n + 1 fragments. No cuts gives one uncut fragment.
    """
    fragments = []
    previous_cut = 0
    for cut in cut_positions:
        fragments.append(cut - previous_cut)
        previous_cut = cut
    fragments.append(sequence_length - previous_cut)
    return fragments


def draw_restriction_map(sequence_length: int, cut_positions: list[int],
                         width: int) -> str:
    """
    Return a text map `width` characters wide: '-' is DNA, '|' is a cut.

    Each cut is scaled into a column so every sequence fits the same width.
    Cuts closer together than one column share a single '|'.
    """
    cut_columns = []
    for cut in cut_positions:
        cut_columns.append(cut * width // sequence_length)
    map_text = ""
    for column in range(width):
        if column in cut_columns:
            map_text = map_text + "|"
        else:
            map_text = map_text + "-"
    return map_text


def make_map_text(sequence_length: int, cut_positions: list[int]) -> str:
    """Return the restriction map with base 1 and the last base labeled."""
    return ("1 " + draw_restriction_map(sequence_length, cut_positions,
                                        MAP_WIDTH)
            + " " + str(sequence_length))


def make_enzyme_row(dna: str, enzyme: RestrictionEnzyme) -> list[str]:
    """Return one results-table row: where this enzyme's sites and cuts are."""
    starts = find_site_starts(dna, enzyme.recognition_sequence)
    if len(starts) == 0:
        return [enzyme.name, enzyme.recognition_sequence, "0", "no sites", "-"]
    return [enzyme.name, enzyme.recognition_sequence, str(len(starts)),
            to_one_based_text(starts),
            numbers_to_text(find_cut_positions(dna, enzyme))]


def make_enzyme_rows(dna: str,
                     enzymes: list[RestrictionEnzyme]) -> list[list[str]]:
    """Return one results-table row per enzyme (called on State's list)."""
    rows = []
    for enzyme in enzymes:
        rows.append(make_enzyme_row(dna, enzyme))
    return rows


def make_fragment_rows(fragment_lengths: list[int]) -> list[list[str]]:
    """Return one row per fragment: its number, first base, last base, size."""
    rows = []
    start_base = 1
    for i in range(len(fragment_lengths)):
        end_base = start_base + fragment_lengths[i] - 1
        rows.append([str(i + 1), str(start_base), str(end_base),
                     str(fragment_lengths[i])])
        start_base = end_base + 1
    return rows


def fragment_check_message(fragment_lengths: list[int],
                           sequence_length: int) -> str:
    """Check that the fragments add back up to the whole sequence."""
    total = sum(fragment_lengths)
    if total == sequence_length:
        return "Check passed: fragments add up to " + str(total) + " bp."
    return ("Warning: fragments add up to " + str(total)
            + " bp, but the sequence is " + str(sequence_length) + " bp.")


def make_comparison_row(dna: str, enzyme: RestrictionEnzyme) -> list[str]:
    """Return one comparison-table row: cuts and smallest/largest fragment."""
    cuts = find_cut_positions(dna, enzyme)
    fragments = calculate_fragment_lengths(len(dna), cuts)
    return [enzyme.name, enzyme.recognition_sequence, str(len(cuts)),
            str(len(fragments)), str(min(fragments)), str(max(fragments))]


def describe_saved_analysis(saved: SavedAnalysis, number: int) -> str:
    """Return the one-line summary shown on the Saved Analyses page."""
    return ("#" + str(number) + " "
            + describe_sequence(saved.sequence_name, saved.dna_sequence)
            + ", enzymes: " + join_with_commas(saved.enzyme_names)
            + ", fragments (bp): " + numbers_to_text(saved.fragment_lengths))


# --- Unit tests for the helpers (expected values worked out by hand) --------

DEMO_DNA = "ATGGAATTCCGGAATTC"
# This 17 bp sequence has GAATTC (EcoRI) starting at indexes 3 and 11.

assert_equal(empty_state(), State("", "", [], [], [], []))

assert_equal(clear_results(State("Demo", DEMO_DNA,
                                 get_enzymes_by_names(["EcoRI"]), [4, 12],
                                 [4, 8, 5], [])),
             State("Demo", DEMO_DNA, [], [], [], []))

assert_equal(get_enzyme_names()[0], "EcoRI")
assert_equal(len(get_enzyme_names()), 13)
assert_equal(get_enzyme_names_of([]), [])
assert_equal(get_enzyme_names_of(get_enzymes_by_names(["SmaI", "BamHI"])),
             ["BamHI", "SmaI"])
assert_equal(get_enzymes_by_names([]), [])
assert_equal(get_enzymes_by_names(["NotARealEnzyme"]), [])

assert_equal(join_with_commas([]), "")
assert_equal(join_with_commas(["EcoRI"]), "EcoRI")
assert_equal(join_with_commas(["EcoRI", "BamHI"]), "EcoRI, BamHI")
assert_equal(numbers_to_text([4, 8, 5]), "4, 8, 5")
assert_equal(numbers_to_text([]), "")
assert_equal(to_one_based_text([3, 11]), "4, 12")
assert_equal(describe_sequence("Demo", DEMO_DNA), "Demo (17 bp)")
assert_equal(describe_enzyme(get_enzymes_by_names(["EcoRI"])[0]),
             "EcoRI (GAATTC): G^AATTC. Leaves a 4-base 5' overhang (sticky).")

assert_equal(clean_sequence("at gc\nAA\t"), "ATGCAA")
assert_equal(clean_sequence(""), "")
assert_equal(find_invalid_characters("ATGC"), [])
assert_equal(find_invalid_characters("ATGXNX"), ["X", "N"])
assert_equal(is_valid_dna("ATGC"), True)
assert_equal(is_valid_dna(""), False)
assert_equal(is_valid_dna("ATGU"), False)

assert_equal(find_site_starts(DEMO_DNA, "GAATTC"), [3, 11])
assert_equal(find_site_starts("ATGC", "GAATTC"), [])
assert_equal(find_site_starts("AAAA", "AA"), [0, 1, 2])
assert_equal(find_site_starts("GAATTC", "GAATTC"), [0])

# EcoRI cuts 1 base into its site (G^AATTC); PstI cuts 5 in (CTGCA^G).
assert_equal(find_cut_positions(DEMO_DNA, get_enzymes_by_names(["EcoRI"])[0]),
             [4, 12])
assert_equal(find_cut_positions("AACTGCAGTT",
                                get_enzymes_by_names(["PstI"])[0]), [7])

# BamHI cuts at 1, EcoRI at 9: the combined list comes back sorted.
assert_equal(combine_cut_positions("GGATCCAAGAATTC",
                                   get_enzymes_by_names(["EcoRI", "BamHI"])),
             [1, 9])
assert_equal(combine_cut_positions("ATGC", get_enzyme_catalog()), [])
# The no-sites sample really has no sites for any catalog enzyme.
assert_equal(combine_cut_positions(
    get_sample_sequences()["No-sites example"], get_enzyme_catalog()), [])
# The pUC19 polylinker has one site each for 10 different catalog enzymes.
assert_equal(len(combine_cut_positions(
    get_sample_sequences()["pUC19 polylinker"], get_enzyme_catalog())), 10)

# Design doc example: cuts at 200 and 700 in 1,000 bp give 200, 500, 300.
assert_equal(calculate_fragment_lengths(1000, [200, 700]), [200, 500, 300])
assert_equal(calculate_fragment_lengths(17, [4, 12]), [4, 8, 5])
assert_equal(calculate_fragment_lengths(10, []), [10])

assert_equal(draw_restriction_map(10, [5], 10), "-----|----")
assert_equal(draw_restriction_map(20, [], 5), "-----")
# 2 * 5 // 10 = 1 and 3 * 5 // 10 = 1, so both cuts share column 1.
assert_equal(draw_restriction_map(10, [2, 3], 5), "-|---")
# 4 * 60 // 17 = 14 and 12 * 60 // 17 = 42, so '|' at columns 14 and 42.
assert_equal(make_map_text(17, [4, 12]),
             "1 --------------|---------------------------|"
             + "----------------- 17")

assert_equal(make_enzyme_row(DEMO_DNA, get_enzymes_by_names(["EcoRI"])[0]),
             ["EcoRI", "GAATTC", "2", "4, 12", "4, 12"])
assert_equal(make_enzyme_row(DEMO_DNA, get_enzymes_by_names(["BamHI"])[0]),
             ["BamHI", "GGATCC", "0", "no sites", "-"])
assert_equal(make_enzyme_rows(DEMO_DNA, []), [])
assert_equal(make_enzyme_rows(DEMO_DNA,
                              get_enzymes_by_names(["EcoRI", "BamHI"])),
             [["EcoRI", "GAATTC", "2", "4, 12", "4, 12"],
              ["BamHI", "GGATCC", "0", "no sites", "-"]])

assert_equal(make_fragment_rows([4, 8, 5]),
             [["1", "1", "4", "4"], ["2", "5", "12", "8"],
              ["3", "13", "17", "5"]])
assert_equal(make_fragment_rows([17]), [["1", "1", "17", "17"]])

assert_equal(fragment_check_message([4, 8, 5], 17),
             "Check passed: fragments add up to 17 bp.")
assert_equal(fragment_check_message([4, 8], 17),
             "Warning: fragments add up to 12 bp, but the sequence is 17 bp.")

assert_equal(make_comparison_row(DEMO_DNA,
                                 get_enzymes_by_names(["EcoRI"])[0]),
             ["EcoRI", "GAATTC", "2", "3", "4", "8"])
assert_equal(make_comparison_row(DEMO_DNA,
                                 get_enzymes_by_names(["BamHI"])[0]),
             ["BamHI", "GGATCC", "0", "1", "17", "17"])

assert_equal(describe_saved_analysis(
    SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5]), 1),
    "#1 Demo (17 bp), enzymes: EcoRI, fragments (bp): 4, 8, 5")


# ---------------------------------------------------------------------------
# 4. Routes
# ---------------------------------------------------------------------------

@route
def index(state: State) -> Page:
    """Home page. (Drafter requires the first page to be named `index`.)"""
    return Page(state, [
        Header("RestrictionMap"),
        "Find where restriction enzymes cut a DNA sequence and what "
        + "fragments they make.",
        Button("Start a new analysis", sequence_input_page),
        Button("Sequence library", sequence_library_page),
        Button("Compare two enzymes", comparison_page),
        Button("Saved analyses", saved_analyses_page),
        Button("About", about_page),
    ])


@route
def sequence_input_page(state: State) -> Page:
    """
    Form for naming and pasting a DNA sequence.

    In Drafter, EVERY button (or link) on a page with form fields sends those
    fields to its route, and Drafter warns when that route has no parameter
    for them. So form pages only get the button whose route uses the fields;
    the pages after it have the Home and navigation buttons.
    """
    return Page(state, [
        Header("Enter a DNA sequence"),
        "Sequence name:",
        TextBox("sequence_name", state.sequence_name),
        "DNA sequence (A, T, G and C only; spaces and line breaks are "
        + "ignored):",
        TextArea("dna_sequence", state.dna_sequence),
        Button("Check sequence", validate_sequence_page),
    ])


@route
def validate_sequence_page(state: State, sequence_name: str,
                           dna_sequence: str) -> Page:
    """Clean and check the submitted sequence, then save it to State."""
    state.sequence_name = sequence_name.strip()
    if state.sequence_name == "":
        state.sequence_name = "Untitled sequence"
    state.dna_sequence = clean_sequence(dna_sequence)
    clear_results(state)

    if state.dna_sequence == "":
        return Page(state, [
            Header("Sequence is empty"),
            "Paste or type a DNA sequence, or load a sample one.",
            Button("Back to sequence input", sequence_input_page),
            Button("Load a sample sequence", sequence_library_page),
            Button("Home", index),
        ])

    invalid = find_invalid_characters(state.dna_sequence)
    if len(invalid) > 0:
        return Page(state, [
            Header("Sequence has invalid characters"),
            "Only A, T, G and C are allowed. Found: "
            + join_with_commas(invalid),
            Button("Fix sequence", sequence_input_page),
            Button("Home", index),
        ])

    return Page(state, [
        Header("Sequence is valid"),
        describe_sequence(state.sequence_name, state.dna_sequence),
        Button("Choose enzymes", enzyme_selection_page),
        Button("Edit sequence", sequence_input_page),
        Button("Home", index),
    ])


@route
def sequence_library_page(state: State) -> Page:
    """List the built-in sample sequences, each with a Load button."""
    samples = get_sample_sequences()
    content = [
        Header("Sequence library"),
        "Load a sample sequence into the input page.",
    ]
    for sample_name in samples:
        content.append(describe_sequence(sample_name, samples[sample_name]))
        content.append(Button("Load " + sample_name, load_sample,
                              [Argument("sample_name", sample_name)]))
    content.append(Button("Back to sequence input", sequence_input_page))
    return Page(state, content)


@route
def load_sample(state: State, sample_name: str) -> Page:
    """Copy a sample into State, then show the input page to check it."""
    samples = get_sample_sequences()
    if sample_name in samples:
        state.sequence_name = sample_name
        state.dna_sequence = samples[sample_name]
        clear_results(state)
    return sequence_input_page(state)


@route
def enzyme_selection_page(state: State) -> Page:
    """Show every enzyme with a checkbox, plus an Analyze button."""
    if not is_valid_dna(state.dna_sequence):
        return Page(state, [
            "Enter a valid DNA sequence before choosing enzymes.",
            Button("Enter a sequence", sequence_input_page),
            Button("Home", index),
        ])

    content = [
        Header("Choose restriction enzymes"),
        "Sequence: " + describe_sequence(state.sequence_name,
                                         state.dna_sequence),
    ]
    for enzyme in get_enzyme_catalog():
        # RelatedCheckBox (a Drafter component): every box shares the name
        # "enzyme_names", and the values of the ticked boxes arrive at
        # analysis_page together as one list[str], e.g. ["EcoRI", "BamHI"].
        already_ticked = enzyme in state.selected_enzymes
        content.append(RelatedCheckBox("enzyme_names", enzyme.name,
                                       already_ticked))
        content.append(describe_enzyme(enzyme))
    # Only one button: it is the route that uses the checkboxes.
    content.append(Button("Analyze", analysis_page))
    return Page(state, content)


@route
def analysis_page(state: State, enzyme_names: list[str]) -> Page:
    """Run the digest for the ticked enzymes and store the results."""
    chosen = get_enzymes_by_names(enzyme_names)
    if len(chosen) == 0:
        return Page(state, [
            "Please select at least one enzyme.",
            Button("Back to enzymes", enzyme_selection_page),
            Button("Edit sequence", sequence_input_page),
            Button("Home", index),
        ])

    state.selected_enzymes = chosen
    state.cut_positions = combine_cut_positions(state.dna_sequence, chosen)
    state.fragment_lengths = calculate_fragment_lengths(
        len(state.dna_sequence), state.cut_positions)
    return results_page(state)


@route
def results_page(state: State) -> Page:
    """Show each enzyme's sites, the fragment table, and the text map."""
    if len(state.selected_enzymes) == 0:
        return Page(state, [
            "No analysis has been run yet.",
            Button("Choose enzymes", enzyme_selection_page),
        ])

    sequence_length = len(state.dna_sequence)
    content = [
        Header("Results: " + state.sequence_name),
        "Sequence length: " + str(sequence_length) + " bp",
        # Loops over State's list of RestrictionEnzyme instances.
        Table(make_enzyme_rows(state.dna_sequence, state.selected_enzymes),
              header=ENZYME_TABLE_HEADER),
    ]

    if len(state.cut_positions) == 0:
        content.append("None of the selected enzymes cut this sequence, so "
                       + "it stays as one uncut piece of "
                       + str(sequence_length) + " bp.")
    else:
        content.append(Header("Fragments", 2))
        content.append(Table(make_fragment_rows(state.fragment_lengths),
                             header=FRAGMENT_TABLE_HEADER))
        content.append(fragment_check_message(state.fragment_lengths,
                                              sequence_length))

    content.append(Header("Restriction map", 2))
    # PreformattedText (a Drafter component) shows text in a fixed-width
    # font, so each '-' and '|' in the map takes up the same space.
    content.append(PreformattedText(make_map_text(sequence_length,
                                                  state.cut_positions)))
    content.append(Button("Save this analysis", save_analysis))
    content.append(Button("Change enzymes", enzyme_selection_page))
    content.append(Button("Compare two enzymes", comparison_page))
    content.append(Button("Home", index))
    return Page(state, content)


@route
def comparison_page(state: State) -> Page:
    """Two dropdowns for picking the enzymes to compare."""
    if not is_valid_dna(state.dna_sequence):
        return Page(state, [
            "Enter a valid DNA sequence before comparing enzymes.",
            Button("Enter a sequence", sequence_input_page),
            Button("Home", index),
        ])

    names = get_enzyme_names()
    return Page(state, [
        Header("Compare two enzymes"),
        "Sequence: " + describe_sequence(state.sequence_name,
                                         state.dna_sequence),
        "First enzyme:",
        SelectBox("first_enzyme", names, names[0]),
        "Second enzyme:",
        SelectBox("second_enzyme", names, names[1]),
        # Only one button: it is the route that uses the two dropdowns.
        Button("Compare", compare_results_page),
    ])


@route
def compare_results_page(state: State, first_enzyme: str,
                         second_enzyme: str) -> Page:
    """Side-by-side table: number of cuts and smallest/largest fragment."""
    if first_enzyme == second_enzyme:
        return Page(state, [
            "Pick two different enzymes to compare.",
            Button("Back to comparison", comparison_page),
            Button("Home", index),
        ])

    rows = []
    for enzyme in get_enzymes_by_names([first_enzyme, second_enzyme]):
        rows.append(make_comparison_row(state.dna_sequence, enzyme))

    return Page(state, [
        Header("Comparison: " + first_enzyme + " vs " + second_enzyme),
        "Sequence: " + describe_sequence(state.sequence_name,
                                         state.dna_sequence),
        Table(rows, header=COMPARE_TABLE_HEADER),
        Button("Compare others", comparison_page),
        Button("Home", index),
    ])


@route
def save_analysis(state: State) -> Page:
    """Save the current results to the history, unless already saved."""
    if len(state.selected_enzymes) > 0:
        new_save = SavedAnalysis(state.sequence_name, state.dna_sequence,
                                 get_enzyme_names_of(state.selected_enzymes),
                                 state.cut_positions, state.fragment_lengths)
        if new_save not in state.saved_analyses:
            state.saved_analyses.append(new_save)
    return saved_analyses_page(state)


@route
def saved_analyses_page(state: State) -> Page:
    """List every saved analysis with Reopen and Delete buttons."""
    content = [Header("Saved analyses")]
    if len(state.saved_analyses) == 0:
        content.append("Nothing saved yet. Run an analysis and click Save "
                       + "on the results page.")
    for position in range(len(state.saved_analyses)):
        number = position + 1
        content.append(describe_saved_analysis(state.saved_analyses[position],
                                               number))
        content.append(Button("Reopen #" + str(number), reopen_analysis,
                              [Argument("position", position)]))
        content.append(Button("Delete #" + str(number), delete_analysis,
                              [Argument("position", position)]))
    content.append(Button("Home", index))
    return Page(state, content)


@route
def reopen_analysis(state: State, position: int) -> Page:
    """Load a saved analysis back into State and show its results."""
    if position >= 0 and position < len(state.saved_analyses):
        saved = state.saved_analyses[position]
        state.sequence_name = saved.sequence_name
        state.dna_sequence = saved.dna_sequence
        state.selected_enzymes = get_enzymes_by_names(saved.enzyme_names)
        state.cut_positions = saved.cut_positions
        state.fragment_lengths = saved.fragment_lengths
        return results_page(state)
    return saved_analyses_page(state)


@route
def delete_analysis(state: State, position: int) -> Page:
    """Remove one saved analysis from the history."""
    if position >= 0 and position < len(state.saved_analyses):
        state.saved_analyses.pop(position)
    return saved_analyses_page(state)


@route
def about_page(state: State) -> Page:
    """Explain what the tool does and what it can't tell you."""
    return Page(state, [
        Header("About RestrictionMap"),
        "RestrictionMap shows the reasoning behind a restriction digest: "
        + "where each enzyme recognizes the DNA, where it cuts, and how the "
        + "fragment lengths add back up to the full sequence.",
        Header("What it does not model", 2),
        BulletedList([
            "Circular DNA such as whole plasmids (sequences are linear here).",
            "Non-palindromic enzymes, or sites with ambiguous bases like N.",
            "Methylation that blocks cutting, buffers, or partial digests.",
            "Whether small fragments would actually be visible on a gel.",
        ]),
        Button("Home", index),
    ])


# ---------------------------------------------------------------------------
# 5. Unit tests for every route
# ---------------------------------------------------------------------------
# A few tests compare one route's page to another route's page (for example,
# load_sample should show exactly the sequence input page). Those other
# routes have their own full tests above them.

# index
assert_equal(index(empty_state()), Page(empty_state(), [
    Header("RestrictionMap"),
    "Find where restriction enzymes cut a DNA sequence and what "
    + "fragments they make.",
    Button("Start a new analysis", sequence_input_page),
    Button("Sequence library", sequence_library_page),
    Button("Compare two enzymes", comparison_page),
    Button("Saved analyses", saved_analyses_page),
    Button("About", about_page),
]))

# sequence_input_page
assert_equal(
    sequence_input_page(State("Demo", DEMO_DNA, [], [], [], [])),
    Page(State("Demo", DEMO_DNA, [], [], [], []), [
        Header("Enter a DNA sequence"),
        "Sequence name:",
        TextBox("sequence_name", "Demo"),
        "DNA sequence (A, T, G and C only; spaces and line breaks are "
        + "ignored):",
        TextArea("dna_sequence", DEMO_DNA),
        Button("Check sequence", validate_sequence_page),
    ]))

# validate_sequence_page: valid, empty, and invalid input
assert_equal(
    validate_sequence_page(empty_state(), " Demo ", "atgg aattc\ncggaattc"),
    Page(State("Demo", DEMO_DNA, [], [], [], []), [
        Header("Sequence is valid"),
        "Demo (17 bp)",
        Button("Choose enzymes", enzyme_selection_page),
        Button("Edit sequence", sequence_input_page),
        Button("Home", index),
    ]))
assert_equal(
    validate_sequence_page(empty_state(), "", "  "),
    Page(State("Untitled sequence", "", [], [], [], []), [
        Header("Sequence is empty"),
        "Paste or type a DNA sequence, or load a sample one.",
        Button("Back to sequence input", sequence_input_page),
        Button("Load a sample sequence", sequence_library_page),
        Button("Home", index),
    ]))
assert_equal(
    validate_sequence_page(empty_state(), "Bad", "ATGXNX"),
    Page(State("Bad", "ATGXNX", [], [], [], []), [
        Header("Sequence has invalid characters"),
        "Only A, T, G and C are allowed. Found: X, N",
        Button("Fix sequence", sequence_input_page),
        Button("Home", index),
    ]))

# sequence_library_page
assert_equal(sequence_library_page(empty_state()), Page(empty_state(), [
    Header("Sequence library"),
    "Load a sample sequence into the input page.",
    "Design doc example (17 bp)",
    Button("Load Design doc example", load_sample,
           [Argument("sample_name", "Design doc example")]),
    "pUC19 polylinker (57 bp)",
    Button("Load pUC19 polylinker", load_sample,
           [Argument("sample_name", "pUC19 polylinker")]),
    "No-sites example (34 bp)",
    Button("Load No-sites example", load_sample,
           [Argument("sample_name", "No-sites example")]),
    Button("Back to sequence input", sequence_input_page),
]))

# load_sample: a real sample, and an unknown name (State unchanged)
assert_equal(
    load_sample(empty_state(), "Design doc example"),
    sequence_input_page(
        State("Design doc example", DEMO_DNA, [], [], [], [])))
assert_equal(load_sample(empty_state(), "Nope"),
             sequence_input_page(empty_state()))

# enzyme_selection_page: blocked until a valid sequence is entered
assert_equal(enzyme_selection_page(empty_state()), Page(empty_state(), [
    "Enter a valid DNA sequence before choosing enzymes.",
    Button("Enter a sequence", sequence_input_page),
    Button("Home", index),
]))

# results_page: no analysis yet, enzymes that cut, and enzymes that don't
assert_equal(results_page(State("Demo", DEMO_DNA, [], [], [], [])),
             Page(State("Demo", DEMO_DNA, [], [], [], []), [
                 "No analysis has been run yet.",
                 Button("Choose enzymes", enzyme_selection_page),
             ]))
assert_equal(
    results_page(State("Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]),
                       [4, 12], [4, 8, 5], [])),
    Page(State("Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]),
               [4, 12], [4, 8, 5], []), [
        Header("Results: Demo"),
        "Sequence length: 17 bp",
        Table([["EcoRI", "GAATTC", "2", "4, 12", "4, 12"]],
              header=ENZYME_TABLE_HEADER),
        Header("Fragments", 2),
        Table([["1", "1", "4", "4"], ["2", "5", "12", "8"],
               ["3", "13", "17", "5"]], header=FRAGMENT_TABLE_HEADER),
        "Check passed: fragments add up to 17 bp.",
        Header("Restriction map", 2),
        PreformattedText("1 --------------|---------------------------|"
                         + "----------------- 17"),
        Button("Save this analysis", save_analysis),
        Button("Change enzymes", enzyme_selection_page),
        Button("Compare two enzymes", comparison_page),
        Button("Home", index),
    ]))
assert_equal(
    results_page(State("Demo", DEMO_DNA, get_enzymes_by_names(["BamHI"]),
                       [], [17], [])),
    Page(State("Demo", DEMO_DNA, get_enzymes_by_names(["BamHI"]),
               [], [17], []), [
        Header("Results: Demo"),
        "Sequence length: 17 bp",
        Table([["BamHI", "GGATCC", "0", "no sites", "-"]],
              header=ENZYME_TABLE_HEADER),
        "None of the selected enzymes cut this sequence, so it stays as "
        + "one uncut piece of 17 bp.",
        Header("Restriction map", 2),
        PreformattedText("1 ------------------------------"
                         + "------------------------------ 17"),
        Button("Save this analysis", save_analysis),
        Button("Change enzymes", enzyme_selection_page),
        Button("Compare two enzymes", comparison_page),
        Button("Home", index),
    ]))

# analysis_page: nothing ticked, then EcoRI ticked
assert_equal(
    analysis_page(State("Demo", DEMO_DNA, [], [], [], []), []),
    Page(State("Demo", DEMO_DNA, [], [], [], []), [
        "Please select at least one enzyme.",
        Button("Back to enzymes", enzyme_selection_page),
        Button("Edit sequence", sequence_input_page),
        Button("Home", index),
    ]))
assert_equal(
    analysis_page(State("Demo", DEMO_DNA, [], [], [], []), ["EcoRI"]),
    results_page(State("Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]),
                       [4, 12], [4, 8, 5], [])))

# comparison_page: blocked without a sequence, then the two dropdowns
assert_equal(comparison_page(empty_state()), Page(empty_state(), [
    "Enter a valid DNA sequence before comparing enzymes.",
    Button("Enter a sequence", sequence_input_page),
    Button("Home", index),
]))
assert_equal(
    comparison_page(State("Demo", DEMO_DNA, [], [], [], [])),
    Page(State("Demo", DEMO_DNA, [], [], [], []), [
        Header("Compare two enzymes"),
        "Sequence: Demo (17 bp)",
        "First enzyme:",
        SelectBox("first_enzyme", get_enzyme_names(), "EcoRI"),
        "Second enzyme:",
        SelectBox("second_enzyme", get_enzyme_names(), "BamHI"),
        Button("Compare", compare_results_page),
    ]))

# compare_results_page: same enzyme twice, then EcoRI vs BamHI
assert_equal(
    compare_results_page(State("Demo", DEMO_DNA, [], [], [], []),
                         "EcoRI", "EcoRI"),
    Page(State("Demo", DEMO_DNA, [], [], [], []), [
        "Pick two different enzymes to compare.",
        Button("Back to comparison", comparison_page),
        Button("Home", index),
    ]))
assert_equal(
    compare_results_page(State("Demo", DEMO_DNA, [], [], [], []),
                         "EcoRI", "BamHI"),
    Page(State("Demo", DEMO_DNA, [], [], [], []), [
        Header("Comparison: EcoRI vs BamHI"),
        "Sequence: Demo (17 bp)",
        Table([["EcoRI", "GAATTC", "2", "3", "4", "8"],
               ["BamHI", "GGATCC", "0", "1", "17", "17"]],
              header=COMPARE_TABLE_HEADER),
        Button("Compare others", comparison_page),
        Button("Home", index),
    ]))

# saved_analyses_page: empty, then one saved analysis
assert_equal(saved_analyses_page(empty_state()), Page(empty_state(), [
    Header("Saved analyses"),
    "Nothing saved yet. Run an analysis and click Save on the results page.",
    Button("Home", index),
]))
assert_equal(
    saved_analyses_page(State("", "", [], [], [], [
        SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])])),
    Page(State("", "", [], [], [], [
        SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])]), [
        Header("Saved analyses"),
        "#1 Demo (17 bp), enzymes: EcoRI, fragments (bp): 4, 8, 5",
        Button("Reopen #1", reopen_analysis, [Argument("position", 0)]),
        Button("Delete #1", delete_analysis, [Argument("position", 0)]),
        Button("Home", index),
    ]))

# save_analysis: saves once, and saving the same thing again adds nothing
assert_equal(
    save_analysis(State("Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]),
                        [4, 12], [4, 8, 5], [])),
    saved_analyses_page(State(
        "Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]), [4, 12], [4, 8, 5],
        [SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])])))
assert_equal(
    save_analysis(State(
        "Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]), [4, 12], [4, 8, 5],
        [SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])])),
    saved_analyses_page(State(
        "Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]), [4, 12], [4, 8, 5],
        [SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])])))

# reopen_analysis: a saved analysis comes back as its results page
assert_equal(
    reopen_analysis(State("", "", [], [], [], [
        SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])]), 0),
    results_page(State(
        "Demo", DEMO_DNA, get_enzymes_by_names(["EcoRI"]), [4, 12], [4, 8, 5],
        [SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])])))
assert_equal(reopen_analysis(empty_state(), 3),
             saved_analyses_page(empty_state()))

# delete_analysis
assert_equal(
    delete_analysis(State("", "", [], [], [], [
        SavedAnalysis("Demo", DEMO_DNA, ["EcoRI"], [4, 12], [4, 8, 5])]), 0),
    saved_analyses_page(empty_state()))

# about_page
assert_equal(about_page(empty_state()), Page(empty_state(), [
    Header("About RestrictionMap"),
    "RestrictionMap shows the reasoning behind a restriction digest: "
    + "where each enzyme recognizes the DNA, where it cuts, and how the "
    + "fragment lengths add back up to the full sequence.",
    Header("What it does not model", 2),
    BulletedList([
        "Circular DNA such as whole plasmids (sequences are linear here).",
        "Non-palindromic enzymes, or sites with ambiguous bases like N.",
        "Methylation that blocks cutting, buffers, or partial digests.",
        "Whether small fragments would actually be visible on a gel.",
    ]),
    Button("Home", index),
]))


# ---------------------------------------------------------------------------
# 6. Start the site
# ---------------------------------------------------------------------------

start_server(empty_state())
