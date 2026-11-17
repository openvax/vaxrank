"""CTA and tissue-risk resolution use Topiary's shared transcript rule."""

from dataclasses import replace

import pytest
from topiary import representative_transcript

from vaxrank.risk_ligand import (
    PROTEIN_RESOLUTION_ALL_ISOFORMS, PROTEIN_RESOLUTION_REPRESENTATIVE,
    resolve_tissue_risk_protein_sequences,
)
from .test_cta_identity import Annotation, ALTERNATE
from .test_cta_expression import run
from .test_risk_ligand import _tissue_risk


@pytest.mark.parametrize('case', ['complete', 'length', 'tie', 'alphabet', 'nmd'])
def test_cta_and_risk_sources_choose_the_same_representative(tmp_path, case):
    genome = Annotation(112)
    ordered = sorted(genome.gene_by_id(ALTERNATE).transcripts,
                     key=lambda t: (-len(t.protein_sequence), t.id))
    transcripts = [ordered[0], ordered[-1]]
    genome.gene_by_id(ALTERNATE).transcripts = transcripts
    long, short = transcripts
    long.complete = short.complete = True
    if case == 'complete':
        long.complete = False
    elif case == 'tie':
        long.protein_sequence = short.protein_sequence
    elif case == 'alphabet':
        long.complete = short.complete = False
        long.protein_sequence = 'X' + long.protein_sequence[1:]
    elif case == 'nmd':
        long.biotype = 'nonsense_mediated_decay'
        long.is_protein_coding = False
    expected = representative_transcript(transcripts)
    if case in ('complete', 'alphabet', 'nmd'):
        assert expected is short
    elif case == 'length':
        assert expected is long
    else:
        assert expected.id == min(t.id for t in transcripts)

    cta = run(tmp_path, genome, [(ALTERNATE, 7)]).admitted_antigens[0]
    assert cta.transcript_ids == (expected.transcript_id,)
    assert cta.amino_acids == expected.protein_sequence
    assert dict(cta.source_metadata)['sequence_selection'] == 'annotation_representative_transcript'

    tissue_risk = _tissue_risk()
    tissue_risk = replace(tissue_risk, proteins=(replace(
        tissue_risk.proteins[0], gene_id=ALTERNATE, gene_name='PRAME'),))
    selected = resolve_tissue_risk_protein_sequences(
        tissue_risk, genome, resolution_policy=PROTEIN_RESOLUTION_REPRESENTATIVE)
    assert len(selected.proteins) == 1
    assert selected.proteins[0].amino_acids == cta.amino_acids
    assert expected.id in selected.proteins[0].transcript_ids
    all_isoforms = resolve_tissue_risk_protein_sequences(
        tissue_risk, genome, resolution_policy=PROTEIN_RESOLUTION_ALL_ISOFORMS)
    assert {p.amino_acids for p in all_isoforms.proteins} == {
        t.protein_sequence for t in transcripts if t.is_protein_coding and 'X' not in t.protein_sequence}


def test_historical_risk_resolution_policy_is_not_rewritten():
    from vaxrank.risk_ligand import PROTEIN_RESOLUTION_LONGEST_ISOFORM
    from .test_risk_ligand import _protein_sequences
    saved = replace(_protein_sequences(), resolution_policy=PROTEIN_RESOLUTION_LONGEST_ISOFORM)
    assert type(saved).from_json(saved.to_json()).resolution_policy == PROTEIN_RESOLUTION_LONGEST_ISOFORM
