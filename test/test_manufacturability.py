from vaxrank.manufacturability import ManufacturabilityScores

def test_c_terminal_proline():
    scores = ManufacturabilityScores.from_amino_acids("A" * 6 + "P")
    assert scores.c_terminal_proline

    scores = ManufacturabilityScores.from_amino_acids("A" * 7)
    assert not scores.c_terminal_proline

def test_n_terminal_cysteine():
    scores = ManufacturabilityScores.from_amino_acids("C" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = ManufacturabilityScores.from_amino_acids(7 * "A")
    assert not scores.difficult_n_terminal_residue

def test_n_terminal_glutamic_acid():
    scores = ManufacturabilityScores.from_amino_acids("E" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = ManufacturabilityScores.from_amino_acids(7 * "A")
    assert not scores.difficult_n_terminal_residue

def test_n_terminal_glutamine():
    scores = ManufacturabilityScores.from_amino_acids("Q" + 6 * "A")
    assert scores.difficult_n_terminal_residue

    scores = ManufacturabilityScores.from_amino_acids(7 * "A")
    assert not scores.difficult_n_terminal_residue

def test_asp_pro_bond_count():
    scores = ManufacturabilityScores.from_amino_acids("A" * 7)
    assert scores.asparagine_proline_bond_count == 0

    scores = ManufacturabilityScores.from_amino_acids("NP" + "A" * 7 + "NP")
    assert scores.asparagine_proline_bond_count == 2

def test_cysteine_count():
    scores = ManufacturabilityScores.from_amino_acids("C" * 7)
    assert scores.cysteine_count == 7

def cterm_7mer_gravy_score():
    scores = ManufacturabilityScores.from_amino_acids("QLFY" + "A" * 7)
    # hydropathy of alanine is 1.8 from Kyte & Doolittle 1982
    assert scores.cterm_7mer_gravy_score == 1.8

def max_7mer_gravy_score():
    scores = ManufacturabilityScores.from_amino_acids("H" * 3 + "A" * 7)
    # hydropathy of alanine is 1.8, histidine is -3.2
    # from Kyte & Doolittle 1982
    assert scores.max_7mer_gravy_score == 1.8
