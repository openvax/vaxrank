


class VaxrankCoreLogic(object):
    def __init__(
            self,
            variants,
            reads_generator,
            mhc_predictor,
            vaccine_peptide_length,
            padding_around_mutation,
            max_vaccine_peptides_per_variant,
            min_alt_rna_reads,
            min_variant_sequence_coverage,
            variant_sequence_assembly,
            num_mutant_epitopes_to_keep=10000,
            min_epitope_score=0.0,
            gene_pathway_check=None):
        """
        Parameters
        ----------
        variants : VariantCollection
            Variant objects to evaluate for vaccine inclusion

        reads_generator : generator
            Yields sequence of varcode.Variant objects paired with sequences of
            AlleleRead objects that support that variant.

        mhc_predictor : mhctools.BasePredictor
            Object with predict_peptides method, used for making pMHC binding
            predictions

        vaccine_peptide_length : int
            Length of vaccine SLP to construct

        padding_around_mutation : int
            Number of off-center windows around the mutation to consider as vaccine
            peptides.

        max_vaccine_peptides_per_variant : int
            Number of vaccine peptides to generate for each mutation.

        min_alt_rna_reads : int
            Drop variant sequences at loci with fewer than this number of reads
            supporting the alt allele.

        min_variant_sequence_coverage : int
            Trim variant sequences to positions supported by at least this number
            of RNA reads.

        variant_sequence_assembly : int
            If True, then assemble variant cDNA sequences based on overlap of RNA
            reads. If False, then variant cDNA sequences must be fully spanned and
            contained within RNA reads.

        num_mutant_epitopes_to_keep : int, optional
            Number of top-ranking epitopes for each vaccine peptide to include in
            computation.

        min_epitope_score : float, optional
            Ignore peptides with binding predictions whose normalized score is less
            than this.

        gene_pathway_check : GenePathwayCheck, optional
            If provided, will check against known pathways/gene sets/variant sets and
            include the info in the all-variants output file.
        """
        self.variants = variants
        self.reads_generator = reads_generator
        self.mhc_predictor = mhc_predictor
        self.vaccine_peptide_length = vaccine_peptide_length
        self.padding_around_mutation = padding_around_mutation
        self.max_vaccine_peptides_per_variant = max_vaccine_peptides_per_variant
        self.min_alt_rna_reads = min_alt_rna_reads
        self.min_variant_sequence_coverage = min_variant_sequence_coverage
        self.variant_sequence_assembly = variant_sequence_assembly
        self.num_mutant_epitopes_to_keep = num_mutant_epitopes_to_keep
        self.min_epitope_score = min_epitope_score
        self.gene_pathway_check = gene_pathway_check

        # will be a dictionary: varcode.Variant -> isovar protein sequence object
        self._isovar_protein_sequence_dict = None

        # will be a dictionary: varcode.Variant -> list(VaccinePeptide)
        self._vaccine_peptides_dict = None

