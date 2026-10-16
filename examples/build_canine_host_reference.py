"""Build the experimental canine CTA host reference from pinned public assets."""

import argparse

from oncoref import canine
from orfdata import Proteome, ReferenceContext


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', help='ORFdata proteome snapshot path')
    parser.add_argument('--cache-root', default=None, help='OncoRef asset cache directory')
    args = parser.parse_args()
    snapshot = canine.cta_snapshot()
    identity = snapshot.as_dict()['references'][0]['identity']
    registered = next(r for r in canine.references.registered_references()
                      if r['identity']['reference_id'] == identity['reference_id'])
    paths = {role: canine.references.ensure(registered['assets'][role], cache_root=args.cache_root)
             for role in ('protein', 'feature_table')}
    reference = ReferenceContext(
        identity['taxon_id'], 'dog', identity['assembly_accession'],
        identity['annotation_source'], identity['annotation_release'],
        assembly_accession=identity['assembly_accession'])
    proteome = Proteome.from_refseq(paths['protein'], paths['feature_table'], reference=reference,
        expected_sha256={role: identity['content_provenance'][role] for role in paths})
    print(proteome.save(args.output))


if __name__ == '__main__':
    main()
