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

"""Public API, imported lazily on first access to keep CLI startup lightweight."""

from importlib import import_module
from typing import TYPE_CHECKING

from .version import __version__

if TYPE_CHECKING:
    from .amino_acids import (
        STANDARD_AMINO_ACIDS as STANDARD_AMINO_ACIDS,
        has_only_standard_amino_acids as has_only_standard_amino_acids,
        validate_amino_acid_sequence as validate_amino_acid_sequence,
    )
    from .antigen_input import (
        AntigenDataset as AntigenDataset,
        predict_antigens as predict_antigens,
    )
    from .antigen_origin import (
        AntigenOrigin as AntigenOrigin,
        CodingReference as CodingReference,
        GenomicEndpoint as GenomicEndpoint,
        JunctionOrigin as JunctionOrigin,
        VariantOrigin as VariantOrigin,
    )
    from .candidate_epitope import (
        CandidateEpitope as CandidateEpitope,
        Peptide as Peptide,
    )
    from .cleavage_inference import (
        audit_pepsickle_inputs as audit_pepsickle_inputs,
        audit_peptidase_inputs as audit_peptidase_inputs,
    )
    from .cleavage_profile import (
        CleavageIntervalEvidence as CleavageIntervalEvidence,
        CleavageProfile as CleavageProfile,
        profile_from_cleavage_result as profile_from_cleavage_result,
        profile_from_residue_scores as profile_from_residue_scores,
    )
    from .construct_audit import (
        ConstructAudit as ConstructAudit,
        audit_construct_sequence as audit_construct_sequence,
        load_construct_sequences as load_construct_sequences,
        save_construct_sequences as save_construct_sequences,
        write_construct_audits as write_construct_audits,
    )
    from .construct_sequence import (
        ConstructChemicalModification as ConstructChemicalModification,
        ConstructEvidence as ConstructEvidence,
        ConstructPlacement as ConstructPlacement,
        ConstructSequence as ConstructSequence,
        ConstructSequenceEdit as ConstructSequenceEdit,
    )
    from .core_logic import (
        run_vaxrank as run_vaxrank,
        vaccine_peptides_for_antigen as vaccine_peptides_for_antigen,
    )
    from .cta_admission import (
        CTAAdmissionAssessment as CTAAdmissionAssessment,
        CTAAdmissionError as CTAAdmissionError,
        CTAAdmissionPolicy as CTAAdmissionPolicy,
        CTAOverrideEvidence as CTAOverrideEvidence,
        CTAReferenceEvidence as CTAReferenceEvidence,
        CTAReferenceResolution as CTAReferenceResolution,
        PatientTumorExpressionEvidence as PatientTumorExpressionEvidence,
        assess_cta_antigen as assess_cta_antigen,
        build_cta_vaccine_antigen as build_cta_vaccine_antigen,
        resolve_cta_reference_evidence as resolve_cta_reference_evidence,
    )
    from .epitope_config import EpitopeConfig as EpitopeConfig
    from .epitope_logic import predict_epitopes as predict_epitopes
    from .host_reference import HostReference as HostReference
    from .identifiers import (
        normalize_ensembl_gene_id as normalize_ensembl_gene_id,
        normalize_mhc_allele as normalize_mhc_allele,
    )
    from .isovar_antigens import antigens_from_isovar as antigens_from_isovar
    from .isovar_fusion import (
        IsovarFusionAntigens as IsovarFusionAntigens,
        fusion_antigens_from_isovar as fusion_antigens_from_isovar,
    )
    from .isovar_input import (
        load_isovar_input as load_isovar_input,
        reconstruct_isovar_antigens as reconstruct_isovar_antigens,
    )
    from .isovar_options import (
        IsovarAdmission as IsovarAdmission,
        IsovarSourceOptions as IsovarSourceOptions,
    )
    from .near_self import (
        BLOSUM62 as BLOSUM62,
        DEFAULT_NEAR_SELF_COMPARATOR as DEFAULT_NEAR_SELF_COMPARATOR,
        AminoAcidSubstitutionMatrix as AminoAcidSubstitutionMatrix,
        HammingRiskComparator as HammingRiskComparator,
        NearSelfAssessment as NearSelfAssessment,
        NearSelfError as NearSelfError,
        NearSelfQuery as NearSelfQuery,
        NearestRiskLigandHit as NearestRiskLigandHit,
        SimilarityComparatorProvenance as SimilarityComparatorProvenance,
        SubstitutionEvidence as SubstitutionEvidence,
        assess_near_self as assess_near_self,
        assess_near_self_queries as assess_near_self_queries,
    )
    from .patient_background import PatientBackgroundPairs as PatientBackgroundPairs
    from .prediction_input import (
        finite_prediction_value as finite_prediction_value,
        prediction_integer as prediction_integer,
    )
    from .proteome import to_proteome as to_proteome
    from .reference_proteome import (
        clear_reference_proteome_caches as clear_reference_proteome_caches,
        ensembl_dataset_cache_identity as ensembl_dataset_cache_identity,
        self_reference_matches as self_reference_matches,
    )
    from .risk_ligand import (
        PROTEIN_RESOLUTION_ALL_ISOFORMS as PROTEIN_RESOLUTION_ALL_ISOFORMS,
        PROTEIN_RESOLUTION_LONGEST_ISOFORM as PROTEIN_RESOLUTION_LONGEST_ISOFORM,
        PROTEIN_RESOLUTION_REPRESENTATIVE as PROTEIN_RESOLUTION_REPRESENTATIVE,
        PatientHLARiskLigandIndex as PatientHLARiskLigandIndex,
        ProteinResolutionCoverage as ProteinResolutionCoverage,
        ResolvedRiskProtein as ResolvedRiskProtein,
        RiskCoverageGap as RiskCoverageGap,
        RiskLigand as RiskLigand,
        RiskLigandCoverage as RiskLigandCoverage,
        RiskLigandError as RiskLigandError,
        RiskLigandIndexProvenance as RiskLigandIndexProvenance,
        RiskLigandPrediction as RiskLigandPrediction,
        RiskLigandSelectionPolicy as RiskLigandSelectionPolicy,
        RiskLigandSource as RiskLigandSource,
        RiskPercentileThreshold as RiskPercentileThreshold,
        RiskProteinSequenceSet as RiskProteinSequenceSet,
        build_patient_hla_risk_ligand_index as build_patient_hla_risk_ligand_index,
        resolve_tissue_risk_protein_sequences as resolve_tissue_risk_protein_sequences,
        risk_ligand_index_from_prediction_frame as risk_ligand_index_from_prediction_frame,
    )
    from .rna_support import RNAFeatureEvidence as RNAFeatureEvidence
    from .safety_assessment import (
        AntigenSafetyAssessment as AntigenSafetyAssessment,
        ConstructBoundary as ConstructBoundary,
        EmittedSafetyLigand as EmittedSafetyLigand,
        MHCRequest as MHCRequest,
        SafetyAssessmentError as SafetyAssessmentError,
        SafetyPrediction as SafetyPrediction,
        SafetyPredictionCoverage as SafetyPredictionCoverage,
        WindowSafetyAssessment as WindowSafetyAssessment,
        assess_antigen_safety as assess_antigen_safety,
        assess_vaccine_antigen_window as assess_vaccine_antigen_window,
        near_self_queries_for_window as near_self_queries_for_window,
        safety_assessment_from_prediction_frame as safety_assessment_from_prediction_frame,
    )
    from .tissue_risk import (
        HPAEvidenceCoverage as HPAEvidenceCoverage,
        HPAReferenceProvenance as HPAReferenceProvenance,
        SafetyTissueGroup as SafetyTissueGroup,
        TissueRiskAssessment as TissueRiskAssessment,
        TissueRiskDerivationError as TissueRiskDerivationError,
        TissueRiskEvidence as TissueRiskEvidence,
        TissueRiskPolicy as TissueRiskPolicy,
        TissueRiskProtein as TissueRiskProtein,
        build_hpa_reference_provenance as build_hpa_reference_provenance,
        derive_tissue_risk_proteins as derive_tissue_risk_proteins,
        resolve_tissue_risk_policy as resolve_tissue_risk_policy,
    )
    from .vaccine_antigen import (
        ANTIGEN_KIND_CTA as ANTIGEN_KIND_CTA,
        ANTIGEN_KIND_ERV as ANTIGEN_KIND_ERV,
        ANTIGEN_KIND_FUSION as ANTIGEN_KIND_FUSION,
        ANTIGEN_KIND_SPLICE as ANTIGEN_KIND_SPLICE,
        AminoAcidInterval as AminoAcidInterval,
        SelfReferenceMatch as SelfReferenceMatch,
        SelfReferenceSource as SelfReferenceSource,
        TargetableMask as TargetableMask,
        TumorSpecificityAttestation as TumorSpecificityAttestation,
        TumorSpecificityEvidence as TumorSpecificityEvidence,
        VaccineAntigen as VaccineAntigen,
    )
    from .vaccine_config import VaccineConfig as VaccineConfig
    from .vaccine_library import (
        antigen_construct_name as antigen_construct_name,
        select_antigen_window as select_antigen_window,
    )
    from .vaccine_peptide import VaccinePeptide as VaccinePeptide


_LAZY_EXPORTS = {
    "run_vaxrank": (".core_logic", "run_vaxrank"),
    "vaccine_peptides_for_antigen": (".core_logic", "vaccine_peptides_for_antigen"),
    "STANDARD_AMINO_ACIDS": (".amino_acids", "STANDARD_AMINO_ACIDS"),
    "has_only_standard_amino_acids": (".amino_acids", "has_only_standard_amino_acids"),
    "validate_amino_acid_sequence": (".amino_acids", "validate_amino_acid_sequence"),
    "normalize_ensembl_gene_id": (".identifiers", "normalize_ensembl_gene_id"),
    "normalize_mhc_allele": (".identifiers", "normalize_mhc_allele"),
    "finite_prediction_value": (".prediction_input", "finite_prediction_value"),
    "prediction_integer": (".prediction_input", "prediction_integer"),
    "EpitopeConfig": (".epitope_config", "EpitopeConfig"),
    "predict_epitopes": (".epitope_logic", "predict_epitopes"),
    "CandidateEpitope": (".candidate_epitope", "CandidateEpitope"),
    "Peptide": (".candidate_epitope", "Peptide"),
    "VaccineConfig": (".vaccine_config", "VaccineConfig"),
    "AntigenOrigin": (".antigen_origin", "AntigenOrigin"),
    "VariantOrigin": (".antigen_origin", "VariantOrigin"),
    "CodingReference": (".antigen_origin", "CodingReference"),
    "GenomicEndpoint": (".antigen_origin", "GenomicEndpoint"),
    "JunctionOrigin": (".antigen_origin", "JunctionOrigin"),
    "ANTIGEN_KIND_CTA": (".vaccine_antigen", "ANTIGEN_KIND_CTA"),
    "ANTIGEN_KIND_ERV": (".vaccine_antigen", "ANTIGEN_KIND_ERV"),
    "ANTIGEN_KIND_FUSION": (".vaccine_antigen", "ANTIGEN_KIND_FUSION"),
    "ANTIGEN_KIND_SPLICE": (".vaccine_antigen", "ANTIGEN_KIND_SPLICE"),
    "AminoAcidInterval": (".vaccine_antigen", "AminoAcidInterval"),
    "SelfReferenceMatch": (".vaccine_antigen", "SelfReferenceMatch"),
    "SelfReferenceSource": (".vaccine_antigen", "SelfReferenceSource"),
    "TargetableMask": (".vaccine_antigen", "TargetableMask"),
    "TumorSpecificityAttestation": (".vaccine_antigen", "TumorSpecificityAttestation"),
    "TumorSpecificityEvidence": (".vaccine_antigen", "TumorSpecificityEvidence"),
    "VaccineAntigen": (".vaccine_antigen", "VaccineAntigen"),
    "VaccinePeptide": (".vaccine_peptide", "VaccinePeptide"),
    "AntigenDataset": (".antigen_input", "AntigenDataset"),
    "predict_antigens": (".antigen_input", "predict_antigens"),
    "HostReference": (".host_reference", "HostReference"),
    "IsovarFusionAntigens": (".isovar_fusion", "IsovarFusionAntigens"),
    "fusion_antigens_from_isovar": (".isovar_fusion", "fusion_antigens_from_isovar"),
    "antigens_from_isovar": (".isovar_antigens", "antigens_from_isovar"),
    "load_isovar_input": (".isovar_input", "load_isovar_input"),
    "reconstruct_isovar_antigens": (".isovar_input", "reconstruct_isovar_antigens"),
    "IsovarAdmission": (".isovar_options", "IsovarAdmission"),
    "IsovarSourceOptions": (".isovar_options", "IsovarSourceOptions"),
    "PatientBackgroundPairs": (".patient_background", "PatientBackgroundPairs"),
    "ConstructChemicalModification": (".construct_sequence", "ConstructChemicalModification"),
    "ConstructEvidence": (".construct_sequence", "ConstructEvidence"),
    "ConstructPlacement": (".construct_sequence", "ConstructPlacement"),
    "ConstructSequence": (".construct_sequence", "ConstructSequence"),
    "ConstructSequenceEdit": (".construct_sequence", "ConstructSequenceEdit"),
    "ConstructAudit": (".construct_audit", "ConstructAudit"),
    "audit_construct_sequence": (".construct_audit", "audit_construct_sequence"),
    "load_construct_sequences": (".construct_audit", "load_construct_sequences"),
    "save_construct_sequences": (".construct_audit", "save_construct_sequences"),
    "write_construct_audits": (".construct_audit", "write_construct_audits"),
    "CleavageProfile": (".cleavage_profile", "CleavageProfile"),
    "CleavageIntervalEvidence": (".cleavage_profile", "CleavageIntervalEvidence"),
    "profile_from_cleavage_result": (".cleavage_profile", "profile_from_cleavage_result"),
    "profile_from_residue_scores": (".cleavage_profile", "profile_from_residue_scores"),
    "audit_pepsickle_inputs": (".cleavage_inference", "audit_pepsickle_inputs"),
    "audit_peptidase_inputs": (".cleavage_inference", "audit_peptidase_inputs"),
    "antigen_construct_name": (".vaccine_library", "antigen_construct_name"),
    "select_antigen_window": (".vaccine_library", "select_antigen_window"),
    "AntigenSafetyAssessment": (".safety_assessment", "AntigenSafetyAssessment"),
    "ConstructBoundary": (".safety_assessment", "ConstructBoundary"),
    "EmittedSafetyLigand": (".safety_assessment", "EmittedSafetyLigand"),
    "MHCRequest": (".safety_assessment", "MHCRequest"),
    "SafetyAssessmentError": (".safety_assessment", "SafetyAssessmentError"),
    "SafetyPrediction": (".safety_assessment", "SafetyPrediction"),
    "SafetyPredictionCoverage": (".safety_assessment", "SafetyPredictionCoverage"),
    "WindowSafetyAssessment": (".safety_assessment", "WindowSafetyAssessment"),
    "assess_antigen_safety": (".safety_assessment", "assess_antigen_safety"),
    "assess_vaccine_antigen_window": (".safety_assessment", "assess_vaccine_antigen_window"),
    "near_self_queries_for_window": (".safety_assessment", "near_self_queries_for_window"),
    "safety_assessment_from_prediction_frame": (
        ".safety_assessment",
        "safety_assessment_from_prediction_frame",
    ),
    "HPAEvidenceCoverage": (".tissue_risk", "HPAEvidenceCoverage"),
    "HPAReferenceProvenance": (".tissue_risk", "HPAReferenceProvenance"),
    "SafetyTissueGroup": (".tissue_risk", "SafetyTissueGroup"),
    "TissueRiskAssessment": (".tissue_risk", "TissueRiskAssessment"),
    "TissueRiskDerivationError": (".tissue_risk", "TissueRiskDerivationError"),
    "TissueRiskEvidence": (".tissue_risk", "TissueRiskEvidence"),
    "TissueRiskPolicy": (".tissue_risk", "TissueRiskPolicy"),
    "TissueRiskProtein": (".tissue_risk", "TissueRiskProtein"),
    "build_hpa_reference_provenance": (".tissue_risk", "build_hpa_reference_provenance"),
    "derive_tissue_risk_proteins": (".tissue_risk", "derive_tissue_risk_proteins"),
    "resolve_tissue_risk_policy": (".tissue_risk", "resolve_tissue_risk_policy"),
    "PROTEIN_RESOLUTION_ALL_ISOFORMS": (".risk_ligand", "PROTEIN_RESOLUTION_ALL_ISOFORMS"),
    "PROTEIN_RESOLUTION_LONGEST_ISOFORM": (".risk_ligand", "PROTEIN_RESOLUTION_LONGEST_ISOFORM"),
    "PROTEIN_RESOLUTION_REPRESENTATIVE": (".risk_ligand", "PROTEIN_RESOLUTION_REPRESENTATIVE"),
    "PatientHLARiskLigandIndex": (".risk_ligand", "PatientHLARiskLigandIndex"),
    "ProteinResolutionCoverage": (".risk_ligand", "ProteinResolutionCoverage"),
    "ResolvedRiskProtein": (".risk_ligand", "ResolvedRiskProtein"),
    "RiskCoverageGap": (".risk_ligand", "RiskCoverageGap"),
    "RiskLigand": (".risk_ligand", "RiskLigand"),
    "RiskLigandCoverage": (".risk_ligand", "RiskLigandCoverage"),
    "RiskLigandError": (".risk_ligand", "RiskLigandError"),
    "RiskLigandIndexProvenance": (".risk_ligand", "RiskLigandIndexProvenance"),
    "RiskLigandPrediction": (".risk_ligand", "RiskLigandPrediction"),
    "RiskLigandSelectionPolicy": (".risk_ligand", "RiskLigandSelectionPolicy"),
    "RiskLigandSource": (".risk_ligand", "RiskLigandSource"),
    "RiskPercentileThreshold": (".risk_ligand", "RiskPercentileThreshold"),
    "RiskProteinSequenceSet": (".risk_ligand", "RiskProteinSequenceSet"),
    "build_patient_hla_risk_ligand_index": (
        ".risk_ligand",
        "build_patient_hla_risk_ligand_index",
    ),
    "resolve_tissue_risk_protein_sequences": (
        ".risk_ligand",
        "resolve_tissue_risk_protein_sequences",
    ),
    "risk_ligand_index_from_prediction_frame": (
        ".risk_ligand",
        "risk_ligand_index_from_prediction_frame",
    ),
    "BLOSUM62": (".near_self", "BLOSUM62"),
    "DEFAULT_NEAR_SELF_COMPARATOR": (".near_self", "DEFAULT_NEAR_SELF_COMPARATOR"),
    "AminoAcidSubstitutionMatrix": (".near_self", "AminoAcidSubstitutionMatrix"),
    "HammingRiskComparator": (".near_self", "HammingRiskComparator"),
    "NearSelfAssessment": (".near_self", "NearSelfAssessment"),
    "NearSelfError": (".near_self", "NearSelfError"),
    "NearSelfQuery": (".near_self", "NearSelfQuery"),
    "NearestRiskLigandHit": (".near_self", "NearestRiskLigandHit"),
    "SimilarityComparatorProvenance": (".near_self", "SimilarityComparatorProvenance"),
    "SubstitutionEvidence": (".near_self", "SubstitutionEvidence"),
    "assess_near_self": (".near_self", "assess_near_self"),
    "assess_near_self_queries": (".near_self", "assess_near_self_queries"),
    "CTAAdmissionAssessment": (".cta_admission", "CTAAdmissionAssessment"),
    "CTAAdmissionError": (".cta_admission", "CTAAdmissionError"),
    "CTAAdmissionPolicy": (".cta_admission", "CTAAdmissionPolicy"),
    "CTAOverrideEvidence": (".cta_admission", "CTAOverrideEvidence"),
    "CTAReferenceEvidence": (".cta_admission", "CTAReferenceEvidence"),
    "CTAReferenceResolution": (".cta_admission", "CTAReferenceResolution"),
    "PatientTumorExpressionEvidence": (".cta_admission", "PatientTumorExpressionEvidence"),
    "assess_cta_antigen": (".cta_admission", "assess_cta_antigen"),
    "build_cta_vaccine_antigen": (".cta_admission", "build_cta_vaccine_antigen"),
    "resolve_cta_reference_evidence": (".cta_admission", "resolve_cta_reference_evidence"),
    "clear_reference_proteome_caches": (
        ".reference_proteome",
        "clear_reference_proteome_caches",
    ),
    "ensembl_dataset_cache_identity": (
        ".reference_proteome",
        "ensembl_dataset_cache_identity",
    ),
    "self_reference_matches": (".reference_proteome", "self_reference_matches"),
    "to_proteome": (".proteome", "to_proteome"),
    "RNAFeatureEvidence": (".rna_support", "RNAFeatureEvidence"),
}

__all__ = ["__version__", *list(_LAZY_EXPORTS)]


def __getattr__(name):
    """PEP 562: load public Vaxrank objects on first attribute access."""
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name, package=__name__), attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    """Include lazy public exports in interactive discovery before first use."""
    return sorted(set(globals()) | set(__all__))
