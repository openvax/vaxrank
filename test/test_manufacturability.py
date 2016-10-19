from vaxrank.manufacturability import compute_manufacturability_scores

def test_c_terminal_proline():
    scores = compute_manufacturability_scores("A" * 6 + "P")
    assert scores.c_terminal_proline

    scores = compute_manufacturability_scores("A" * 7)
    assert not scores.c_terminal_proline

def test_n_terminal_cysteine():
    scores = compute_manufacturability_scores("C" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = compute_manufacturability_scores(7 * "A")
    assert not scores.difficult_n_terminal_residue

def test_n_terminal_glutamic_acid():
    scores = compute_manufacturability_scores("E" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = compute_manufacturability_scores(7 * "A")
    assert not scores.difficult_n_terminal_residue

def test_n_terminal_glutamine():
    scores = compute_manufacturability_scores("Q" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = compute_manufacturability_scores(7 * "A")
    assert not scores.difficult_n_terminal_residue

def test_asp_pro_bond_count():
    scores = compute_manufacturability_scores("A" * 7)
    assert scores.asparagine_proline_bond_count == 0

    scores = compute_manufacturability_scores("NP" + "A" * 7 + "NP")
    assert scores.asparagine_proline_bond_count == 2

def test_cysteine_count():
    scores = compute_manufacturability_scores("C" * 7)
    assert scores.cysteine_count == 7

def cterm_7mer_gravy_score():
    scores = compute_manufacturability_scores("QLFY" + "A" * 7)
    # hydropathy of alanine is 1.8 from Kyte & Doolittle 1982
    assert scores.cterm_7mer_gravy_score == 1.8

def max_7mer_gravy_score():
    scores = compute_manufacturability_scores("H" * 3 + "A" * 7)
    # hydropathy of alanine is 1.8, histidine is -3.2
    # from Kyte & Doolittle 1982
    assert scores.max_7mer_gravy_score == 1.8
