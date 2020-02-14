from serializable import Serializable

class PatientInfo(Serializable):
    def __init__(
            self,
            patient_id,
            vcf_paths,
            bam_path,
            mhc_alleles,
            num_somatic_variants,
            num_coding_effect_variants,
            num_variants_with_rna_support,
            num_variants_with_vaccine_peptides):
        self.patient_id = patient_id
        self.vcf_paths = vcf_paths
        self.bam_path = bam_path
        self.mhc_alleles = mhc_alleles
        self.num_somatic_variants = num_somatic_variants
        self.num_coding_effect_variants = num_coding_effect_variants
        self.num_variants_with_rna_support = num_variants_with_rna_support
        self.num_variants_with_vaccine_peptides = num_variants_with_vaccine_peptides
