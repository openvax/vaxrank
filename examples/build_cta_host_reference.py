"""Build an annotated CTA host reference from OncoRef's pinned source inventory."""

import argparse

from vaxrank.cta_reference_adapter import build_cta_host_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('species', help='Species name or NCBI taxonomy ID')
    parser.add_argument('output', help='ORFdata proteome snapshot path')
    parser.add_argument('--evidence-snapshot', help='Frozen OncoRef CTA snapshot')
    parser.add_argument('--cache-root', help='OncoRef asset cache directory')
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--max-download-bytes', type=int, default=0,
                        help='Total missing-source transfer budget; default: cached assets only')
    args = parser.parse_args()
    host = build_cta_host_reference(args.species, evidence_snapshot=args.evidence_snapshot,
        cache_root=args.cache_root, offline=args.offline, max_bytes=args.max_download_bytes)
    print(host.save(args.output))


if __name__ == '__main__':
    main()
