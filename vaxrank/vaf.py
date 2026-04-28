# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""DNA VAF extraction from VCF FORMAT/sample data.

varcode's ``load_vcf`` keeps per-variant INFO and per-sample FORMAT fields
on ``VariantCollection.source_to_metadata_dict``. This module reads that
metadata and returns a ``{Variant: float}`` mapping of DNA VAF values.
"""

import logging

logger = logging.getLogger(__name__)


def _parse_af(value, alt_allele_index=0):
    """Parse a FORMAT 'AF' value, which may be a float, list, or string.

    For multiallelic records, AF is conventionally one value per ALT
    allele; ``alt_allele_index`` (0-based into the alt-only list) selects
    the right entry. Scalar values are returned as-is regardless of
    ``alt_allele_index``.
    """
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if alt_allele_index < 0 or alt_allele_index >= len(value):
            return None
        value = value[alt_allele_index]
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _vaf_from_ad(value, alt_allele_index):
    """
    Compute VAF from a FORMAT 'AD' (allelic depth) array.

    AD is conventionally [ref_depth, alt1_depth, alt2_depth, ...].
    ``alt_allele_index`` is the 0-based index into the alt-only list (so
    AD index = alt_allele_index + 1). VAF is reported as
    ``alt_K / sum(AD)`` (matches GATK's convention for FORMAT/AF —
    fraction of total reads supporting this specific alt).

    Returns None if AD is absent, contains a None entry (filtering
    those out would shift indices and silently change the meaning of
    ``alt_allele_index``), or sums to zero.
    """
    if not isinstance(value, (list, tuple)):
        return None
    if any(x is None for x in value):
        return None
    try:
        depths = [int(x) for x in value]
    except (TypeError, ValueError):
        return None
    if not depths or sum(depths) == 0:
        return None
    alt_idx = alt_allele_index + 1
    if alt_idx >= len(depths):
        return None
    return depths[alt_idx] / sum(depths)


# Sentinel returned when a tumor_sample_name was supplied but doesn't appear in
# the VCF's FORMAT samples — distinct from "no usable sample" so the caller can
# warn about typos rather than silently dropping every variant.
_SAMPLE_NAME_NOT_FOUND = object()


def _pick_tumor_sample(sample_info, tumor_sample_name):
    """Pick the sample dict to read AF/AD from.

    Returns ``None`` when no usable sample is found (no FORMAT data, or
    a multi-sample VCF without ``--tumor-sample-name``). Returns
    ``_SAMPLE_NAME_NOT_FOUND`` when a name was supplied but doesn't
    appear in the FORMAT dict.
    """
    if not sample_info:
        return None
    if tumor_sample_name:
        if tumor_sample_name not in sample_info:
            return _SAMPLE_NAME_NOT_FOUND
        return sample_info[tumor_sample_name]
    if len(sample_info) == 1:
        return next(iter(sample_info.values()))
    return None


def extract_dna_vaf_by_variant(variant_collection, tumor_sample_name=None):
    """
    Build a ``{Variant: dna_vaf}`` mapping by reading per-sample FORMAT
    data attached to ``variant_collection.source_to_metadata_dict``.

    Prefers the FORMAT/AF field; falls back to AD when AF is absent.
    Variants without usable depth/frequency data are omitted.

    Parameters
    ----------
    variant_collection : varcode.VariantCollection
    tumor_sample_name : str, optional
        Name of the FORMAT-column sample to read. Required for multi-sample
        VCFs; auto-picked when only one sample is present.

    Returns
    -------
    dict[varcode.Variant, float]
    """
    metadata = getattr(variant_collection, 'source_to_metadata_dict', None) or {}
    if not metadata:
        return {}

    dna_vaf = {}
    skipped_missing = 0
    skipped_ambiguous = 0
    skipped_unmatched_name = 0
    available_sample_names = set()
    for _source, per_variant in metadata.items():
        for variant, info in per_variant.items():
            sample_info = info.get('sample_info') or {}
            available_sample_names.update(sample_info)
            sample = _pick_tumor_sample(sample_info, tumor_sample_name)
            if sample is _SAMPLE_NAME_NOT_FOUND:
                skipped_unmatched_name += 1
                continue
            if sample is None:
                if sample_info and len(sample_info) > 1:
                    skipped_ambiguous += 1
                continue
            alt_idx = info.get('alt_allele_index', 0)
            vaf = _parse_af(sample.get('AF'), alt_idx)
            if vaf is None:
                vaf = _vaf_from_ad(sample.get('AD'), alt_idx)
            if vaf is None:
                skipped_missing += 1
                continue
            dna_vaf[variant] = vaf

    if skipped_unmatched_name:
        logger.warning(
            "DNA VAF: --tumor-sample-name '%s' did not match any FORMAT "
            "sample in the VCF; %d variant(s) skipped. Available samples: %s.",
            tumor_sample_name,
            skipped_unmatched_name,
            ", ".join(sorted(available_sample_names)) or "(none)",
        )
    if skipped_ambiguous:
        logger.warning(
            "DNA VAF: %d variant(s) skipped because the VCF has multiple "
            "samples but --tumor-sample-name was not provided.",
            skipped_ambiguous,
        )
    if skipped_missing:
        logger.info(
            "DNA VAF: %d variant(s) lacked a usable AF/AD field.",
            skipped_missing,
        )
    return dna_vaf
