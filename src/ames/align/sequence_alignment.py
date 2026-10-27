from Bio import Align

aligner = Align.PairwiseAligner()
blosum62 = Align.substitution_matrices.load("BLOSUM62")
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -11 
aligner.extend_gap_score = -1

def sequence_identity(query, target):
    alignment = aligner.align(query, target)
    a1 = alignment[0]
    aligned_residues = a1.counts().identities
    sequence_identity = aligned_residues / a1.length
    return round(sequence_identity, 3)


