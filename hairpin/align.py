from Bio.Align.substitution_matrices import Array
from Bio.Align import PairwiseAligner
import numpy as np


class Alignment:
    alphabet = "ACGTN"
    rule_matrix_1 = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    rule_matrix_2 = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    def __init__(self, sub_matrix, mode='global', open_gap_score=-3, extend_gap_score=-2):
        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score
        if sub_matrix == 'rule_matrix_1':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_1)
        elif sub_matrix == 'rule_matrix_2':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_2)
        else:
            raise ValueError(f"Invalid sub_matrix argument: {sub_matrix}")
        self.aligner = aligner
        self.type = sub_matrix


    def align(self, seq1, seq2):
        alignments = self.aligner.align(seq1[::-1],seq2[::-1])
        align_read1, align_read2 = alignments[0][0][::-1], alignments[0][1][::-1]
        if self.type == 'rule_matrix_1':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): -20, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
            }

            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
            align_read1 = align_read1[:max_index]
            align_read2 = align_read2[:max_index]

            result = [c2 if ((c1 == c2) or (c1 == 'T' and c2 == 'C')) else 'N' for c1, c2 in zip(align_read1, align_read2)]

        elif self.type == 'rule_matrix_2':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): 1, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10,  ('-', '-'): -10
            }
            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
            align_read1 = align_read1[:max_index]
            align_read2 = align_read2[:max_index]

            result = [c2 if ((c1 == c2) or (c1 == 'T' and c2 == 'C') ) else c1 if (c1 == 'G' and c2 == 'A') else 'N' for
                      c1, c2 in zip(align_read1, align_read2)]
        else:
            raise ValueError(f"Invalid sub_matrix argument: {self.type}")

        return ''.join(result), align_read1, align_read2









