import time
from typing import List

IUPAC_DNA_A = "A"
IUPAC_DNA_C = "C"
IUPAC_DNA_G = "G"
IUPAC_DNA_T = "TU"
IUPAC_DNA_U = "UT"
IUPAC_DNA_R = "AG"
IUPAC_DNA_Y = "CTU"
IUPAC_DNA_S = "GC"
IUPAC_DNA_W = "ATU"
IUPAC_DNA_K = "GTU"
IUPAC_DNA_M = "AC"
IUPAC_DNA_B = "CGTU"
IUPAC_DNA_D = "AGTU"
IUPAC_DNA_H = "ACTU"
IUPAC_DNA_V = "ACG"
IUPAC_DNA_N = "AGCTU"
IUPAC_DNA_DOT = "."

IUPAC_ALL = "ACGTURYSWKMBDHVN."
IUPAC_WILDCARD = "RYSWKMBDHVN"
IUPAC_NO_WILDCARD = "AGCTU"
IUPAC_DNA_BASE = 'AGCT'
IUPAC_RNA_BASE = 'AGCU'

IUPAC_XNA_MAP = {'A': IUPAC_DNA_A, 'C': IUPAC_DNA_C, 'G': IUPAC_DNA_G, 'T': IUPAC_DNA_T, 'U': IUPAC_DNA_U,
                 'R': IUPAC_DNA_R,
                 'Y': IUPAC_DNA_Y, 'S': IUPAC_DNA_S, 'W': IUPAC_DNA_W, 'K': IUPAC_DNA_K, 'M': IUPAC_DNA_M,
                 'B': IUPAC_DNA_B,
                 'D': IUPAC_DNA_D, 'H': IUPAC_DNA_H, 'V': IUPAC_DNA_V, 'N': IUPAC_DNA_N, '.': IUPAC_DNA_DOT}


# Algorithm definitions
# match - match from sequence one and sequence two
# mismatch - delete from sequence one and sequence two
# insert - insert from sequence one (gap on sequence two)
# delete - insert from sequence two (gap on sequence one)
def is_valid_sequence(sequence: str, inc_wildcard: bool = True) -> bool:
    legal_chars = IUPAC_ALL if inc_wildcard else IUPAC_NO_WILDCARD
    for c in sequence.upper():
        if c not in legal_chars:
            return False
    return True


def common_dna_code(dna_one: str, dna_two: str) -> List[str]:
    code_one = IUPAC_XNA_MAP.get(dna_one.upper())
    code_two = IUPAC_XNA_MAP.get(dna_two.upper())
    if code_one is None:
        raise ValueError("{} is not a know IUPAC DNA code".format(dna_one))
    if code_two is None:
        raise ValueError("{} is not a know IUPAC DNA code".format(dna_two))
    return [c for c in code_one if c in code_two]


def agree(char_one, char_two, match_score_N, match_score, mismatch_penalty):
    if common_dna_code(char_one, char_two):
        if char_two.upper() == 'N':
            return match_score_N
        else:
            return match_score
    return mismatch_penalty


# Object containing functions that control alignment scores for each action (match, deletion and insertion)
# Each function receives four objects: sequence one with it's current index and sequence two with it's current index
class SequenceAlignmentScore:
    def __init__(self, minmax_func=max, match_func=lambda s1, i1, s2, i2: agree(s1[i1], s2[i2], 1, 1, -1),
                 delete_func=lambda s1, i1, s2, i2: 0, insert_func=lambda s1, i1, s2, i2: 0):
        self.minmax_func = minmax_func
        self.match_func = match_func
        self.delete_func = delete_func
        self.insert_func = insert_func


def align_iupac_dna_sequence(seq_one, seq_two, sequence_alignment_score=SequenceAlignmentScore()):
    score_matrix = [([0] * (len(seq_two) + 1)) for i in range(0, len(seq_one) + 1)]
    for one_index in range(0, len(seq_one) + 1):
        for two_index in range(0, len(seq_two) + 1):
            if one_index == 0 and two_index == 0:
                continue
            choose_list = []
            # Match or mismatch (score decided by agree)
            if one_index > 0 and two_index > 0:
                score = sequence_alignment_score.match_func(seq_one, one_index - 1, seq_two, two_index - 1)
                choose_list.append(score_matrix[one_index - 1][two_index - 1] + score)
            # insert from seq 1
            if one_index > 0:
                insert_penalty = sequence_alignment_score.insert_func(seq_one, one_index - 1, seq_two,
                                                                      max(two_index - 1, 0))
                choose_list.append(score_matrix[one_index - 1][two_index] + insert_penalty)
            # "delete" from seq 2
            if two_index > 0:
                # N is wildcard, deletion of it should be different
                delete_penalty = sequence_alignment_score.delete_func(seq_one, max(one_index - 1, 0), seq_two,
                                                                      two_index - 1)
                choose_list.append(score_matrix[one_index][two_index - 1] + delete_penalty)
            score_matrix[one_index][two_index] = sequence_alignment_score.minmax_func(choose_list)
    return score_matrix


def get_best_score(score_matrix):
    return score_matrix[len(score_matrix) - 1][len(score_matrix[0]) - 1]


class IUPACAlignmentError(Exception):
    def __init__(self, message):
        self.message = message


# returns all optimal alignments
def generate_optimal_alignments(seq_one, seq_two, score_matrix=None, sequence_alignment_score=SequenceAlignmentScore(),
                                return_single=False):
    if score_matrix is None:
        score_matrix = align_iupac_dna_sequence(seq_one, seq_two, sequence_alignment_score)
    results = []
    alignment_stack = [("", "", (len(seq_one), len(seq_two)))]
    while len(alignment_stack) > 0:
        align_one, align_two, matrix_index = alignment_stack.pop()
        # print("Running from stack: {}\n{}\n{}".format(matrix_index, align_one, align_two))
        while matrix_index[0] > 0 or matrix_index[1] > 0:
            # print("start round: {}\n{}\n{}".format(matrix_index, align_one, align_two))
            matched_this_round = False
            round_align_one = None
            round_align_two = None
            round_matrix_index = None
            curr_score = score_matrix[matrix_index[0]][matrix_index[1]]
            # mismatch / match
            if matrix_index[0] > 0 and matrix_index[1] > 0:
                change = sequence_alignment_score.match_func(seq_one, matrix_index[0] - 1, seq_two, matrix_index[1] - 1)
                if curr_score == score_matrix[matrix_index[0] - 1][matrix_index[1] - 1] + change:
                    round_align_one = seq_one[matrix_index[0] - 1] + align_one
                    round_align_two = seq_two[matrix_index[1] - 1] + align_two
                    matched_this_round = True
                    round_matrix_index = (matrix_index[0] - 1, matrix_index[1] - 1)
            # insert
            if matrix_index[0] > 0:
                change = sequence_alignment_score.insert_func(seq_one, matrix_index[0] - 1, seq_two,
                                                              max(matrix_index[1] - 1, 0))
                if curr_score == score_matrix[matrix_index[0] - 1][matrix_index[1]] + change:
                    if matched_this_round and not return_single:
                        alignment_stack.insert(0, (seq_one[matrix_index[0] - 1] + align_one, '-' + align_two,
                                                   (matrix_index[0] - 1, matrix_index[1])))
                    elif not matched_this_round:
                        round_align_one = seq_one[matrix_index[0] - 1] + align_one
                        round_align_two = '-' + align_two
                        round_matrix_index = (matrix_index[0] - 1, matrix_index[1])
                        matched_this_round = True
            # delete
            if matrix_index[1] > 0:
                change = sequence_alignment_score.delete_func(seq_one, max(matrix_index[0] - 1, 0), seq_two,
                                                               matrix_index[1] - 1)
                if curr_score == score_matrix[matrix_index[0]][matrix_index[1] - 1] + change:
                    if matched_this_round and not return_single:
                        alignment_stack.insert(0, ('-' + align_one, seq_two[matrix_index[1] - 1] + align_two,
                                                   (matrix_index[0], matrix_index[1] - 1)))
                    elif not matched_this_round:
                        round_align_one = '-' + align_one
                        round_align_two = seq_two[matrix_index[1] - 1] + align_two
                        round_matrix_index = (matrix_index[0], matrix_index[1] - 1)
                        matched_this_round = True
            if not matched_this_round:
                # Should NEVER arrive here, to check that we selected one of the options
                raise IUPACAlignmentError("Score matrix, no match! {}\n{}".format((matrix_index[0], matrix_index[1])
                                                                                  , score_matrix))
            align_one = round_align_one
            align_two = round_align_two
            matrix_index = round_matrix_index
        # print("End stack:\n{}\n{}".format(align_one, align_two))
        results.append((align_one, align_two))
    return results


if __name__ == "__main__":
    def test(i: int, seq_one: str, seq_two: str, sequence_alignment_score: SequenceAlignmentScore):
        print("Test {}\n===========================".format(i))
        start_time = time.process_time()
        my_score_matrix = align_iupac_dna_sequence(seq_one, seq_two, sequence_alignment_score=sequence_alignment_score)
        print("Score matrix ({}) took {} seconds:".format(get_best_score(my_score_matrix),
                                                          (time.process_time() - start_time) / 1000))
        if print_score_matrix:
            print(my_score_matrix)
        start_time = time.process_time()
        alignment_list = generate_optimal_alignments(seq_one, seq_two, score_matrix=my_score_matrix,
                                                     return_single=single_res,
                                                     sequence_alignment_score=sequence_alignment_score)
        print("Generating {} optimal alignment{}: took {} seconds".format('single' if single_res else 'all',
                                                                          '' if single_res else 's',
                                                                          (time.process_time() - start_time) / 1000))
        for index, alignment_str in zip(range(0, len(alignment_list)), alignment_list):
            print("alignment {}:\n{}\n{}".format(index, alignment_str[0], alignment_str[1]))
        print('===========================\n')
    '''
    seq_a = "AGGUAGGCAUCGCGCGTTCGUACTAGTCGATCAUTGCATCGACGGGATUCGAAGCA"    # AGGUAGGCA---UCGCGCG-TT-C-GUAC-TAG--TCGATC--AUTGC-ATCGACGGGATUCGAAGCA-
    seq_b = "NNNNNNNRRRRNNNNNGGUUGCCUACGUAGGGCUCCCGCCNNNNNNNKKHHHHNNNCAU" # NNNNNNN-RRRRNNNNN-GGUUGCC-UACGUAGGG-C--UCCC---GCCNNNNNNNKKHHHHNNN-CAU
    seq_a = "AGCGUACGCGCGAAGAAAACCCUCGCCGA" # AGCGUACGCGCGAAGAAAACCCUCGC--CGA---
    seq_b = "NNNNNGGGNNNNAAACCNNNNNUUHHHHH" # NNNNN--G-G-GNNNNAAACCNNNNNUUH-HHHH
    seq_a = "ACCGGUC" # ACCGGUC
    seq_b = "NGGN"    # N--GGN-
    seq_a = "AGGUAGGCAUCGCGCGTTCGUACTAGTCGATCAUTGCATCGACGGGATUCGAAGCA"  # AGGUAGGCA---UCGCGCG-TT-C-GUAC-TAG--TCGATC--AUTGC-ATCGACGGGATUCGAAGCA-
    seq_b = "NNNNNNNRRRRNNNNNGGUUGCCUACGUAGGGCUCCCGCCNNNNNNNKKHHHHNNNCAU"  # NNNNNNN-RRRRNNNNN-GGUUGCC-UACGUAGGG-C--UCCC---GCCNNNNNNNKKHHHHNNN-CAU
    seq_a = "UGGGCCCGUUUGGGACGCCAAUUCAGCGUUGCGUUACUACCACCGCGCGUAGGGCGGUGAUCGGGCCCA"
    seq_b = "NNNNNNNNUNNNNNNNNNNNNNNNNNNNNNNNNUNNNUNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNN"
    seq_a = ''
    seq_b = 'NNGCNN'
    '''
    count = 0
    seq_a = "AGGUAGGCAUCGCGCGTTCGUACTAGTCGATCAUTGCATCGACGGGATUCGAAGCA"  # AGGUAGGCA---UCGCGCG-TT-C-GUAC-TAG--TCGATC--AUTGC-ATCGACGGGATUCGAAGCA-
    seq_b = "NNNNNNNRRRRNNNNNGGUUGCCUACGUAGGGCUCCCGCCNNNNNNNKKHHHHNNNCAU"  # NNNNNNN-RRRRNNNNN-GGUUGCC-UACGUAGGG-C--UCCC---GCCNNNNNNNKKHHHHNNN-CAU
    single_res = True
    print_score_matrix = False

    print("starting alignment:\n{}\n{}".format(seq_a, seq_b))
    # Test 0
    sequence_alignment_object = SequenceAlignmentScore()
    test(0, seq_a, seq_b, sequence_alignment_object)

    # Test 1
    sequence_alignment_object = SequenceAlignmentScore(
        delete_func=lambda s1, i1, s2, i2: -50,
        match_func=lambda s1, i1, s2, i2: agree(s1[i1], s2[i2], 1, 0, -100))
    test(1, seq_a, seq_b, sequence_alignment_object)

    # Test 2
    sequence_alignment_object = SequenceAlignmentScore(
        delete_func=lambda s1, i1, s2, i2: -1 if s2[i2].upper() == 'N' else -100,
        insert_func=lambda s1, i1, s2, i2: -1,
        match_func=lambda s1, i1, s2, i2: agree(s1[i1], s2[i2], 0, 0, -100))
    test(2, seq_a, seq_b, sequence_alignment_object)

    # Test 3
    sequence_alignment_object = SequenceAlignmentScore(
        match_func=lambda s1, i1, s2, i2: agree(s1[i1], s2[i2], 1, 100, 0))
    test(3, seq_a, seq_b, sequence_alignment_object)

    # Test 4
    sequence_alignment_object = SequenceAlignmentScore(
        delete_func=lambda s1, i1, s2, i2: 1 if s2[i2].upper() == 'N' else 1000,
        insert_func=lambda s1, i1, s2, i2: 1,
        match_func=lambda s1, i1, s2, i2: agree(s1[i1], s2[i2], 0, 0, 10000),
        minmax_func=min)
    test(4, seq_a, seq_b, sequence_alignment_object)

    # Test 5 should be same as 4
    def del_lower(s1, i1, s2, i2):
        if s2[i2] == 'N':
            return 1
        elif s2[i2] == 'n':
            if 0 < i2 < len(s2) - 1 and s2[i2 - 1] > 'Z' and s2[i2 + 1] > 'Z':
                return 20
        else:
            return 1000

    sequence_alignment_object = SequenceAlignmentScore(
        delete_func=del_lower,
        insert_func=lambda s1, i1, s2, i2: 20 if 0 < i2 < len(s2) - 1 and s2[i2 - 1] > 'Z' and s2[i2] > 'Z' else 1,
        match_func=lambda s1, i1, s2, i2: agree(s1[i1], s2[i2], 0, 0, 1000),
        minmax_func=min)
    test(5, seq_a, seq_b, sequence_alignment_object)

    # Test 6 - motifs
    seq_b = "NNNNNNNrrrrNNNNNGGUUGCCUACGUAGGGCUCCCGCCNNNNNNNKkhhhhnnncaU"
    test(6, seq_a, seq_b, sequence_alignment_object)

    # Test 7 - bug
    seq_a = 'CAGUGU'
    seq_b = 'auunng'
    test(7, seq_a, seq_b, sequence_alignment_object)
