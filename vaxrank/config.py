import msgspec

class Config(msgspec.Struct):

    """Parameters for score, filtering, and ranking both epitopes and vaccine peptides"""
    logistic_epitope_score_midpoint : float = 350.0
    logistic_epitope_score_width : float = 150.0
    
    min_epitope_score : float = 0.0
    binding_affinity_cutoff : float = 5000.0
    