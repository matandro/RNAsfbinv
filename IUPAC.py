import time

IUPAC_DNA_A = "A"
IUPAC_DNA_C = "C"
IUPAC_DNA_G = "G"
IUPAC_DNA_T = "TU"
IUPAC_DNA_U = "TU"
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
IUPAC_DNA_BASE = 'AGCT'
IUPAC_RNA_BASE = 'AGCU'

IUPAC_DNA_MAP = {'A': IUPAC_DNA_A, 'C': IUPAC_DNA_C, 'G': IUPAC_DNA_G, 'T': IUPAC_DNA_T, 'U': IUPAC_DNA_U,
                 'R': IUPAC_DNA_R,
                 'Y': IUPAC_DNA_Y, 'S': IUPAC_DNA_S, 'W': IUPAC_DNA_W, 'K': IUPAC_DNA_K, 'M': IUPAC_DNA_M,
                 'B': IUPAC_DNA_B,
                 'D': IUPAC_DNA_D, 'H': IUPAC_DNA_H, 'V': IUPAC_DNA_V, 'N': IUPAC_DNA_N, '.': IUPAC_DNA_DOT}


# Algorithm definitions
# match - match from sequence one and sequence two
# mismatch - delete from sequence one and sequence two
# insert - insert from sequence one (gap on sequence two)
# delete - insert from sequence two (gap on sequence one)


def check_valid_sequence(sequence):
    for c in sequence.upper():
        if c not in IUPAC_ALL:
            return False
    return True


def common_dna_code(dna_one, dna_two):
    code_one = IUPAC_DNA_MAP.get(dna_one)
    code_two = IUPAC_DNA_MAP.get(dna_two)
    if code_one is None:
        raise ValueError("{} is not a know IUPAC DNA code".format(dna_one))
    if code_two is None:
        raise ValueError("{} is not a know IUPAC DNA code".format(dna_two))
    return [c for c in code_one if c in code_two]


def agree(char_one, char_two, match_score_N, match_score, mismatch_penalty):
    if len(common_dna_code(char_one, char_two)) > 0:
        if char_two == 'N':
            return match_score_N
        else:
            return match_score
    return mismatch_penalty


def align_iupac_dna_sequence(seq_one, seq_two, match_score=1, delete_penalty=0, insert_penalty=0, mismatch_penalty=-1,
                             panelty_del_N=0, match_score_N=1, minmax_func=max):
    score_matrix = [([0] * (len(seq_two) + 1)) for i in range(0, len(seq_one) + 1)]
    for one_index in range(0, len(seq_one) + 1):
        for two_index in range(0, len(seq_two) + 1):
            if one_index == 0 and two_index == 0:
                continue
            choose_list = []
            # Match or mismatch (score decided by agree)
            if one_index > 0 and two_index > 0:
                score = agree(seq_one[one_index - 1], seq_two[two_index - 1],
                              match_score_N, match_score, mismatch_penalty)
                choose_list.append(score_matrix[one_index - 1][two_index - 1] + score)
            # insert from seq 1
            if one_index > 0:
                choose_list.append(score_matrix[one_index - 1][two_index] + insert_penalty)
            # "delete" from seq 2
            if two_index > 0:
                # N is wildcard, deletion of it should be different
                paneltry = panelty_del_N if seq_two[two_index - 1] == 'N' else delete_penalty
                choose_list.append(score_matrix[one_index][two_index - 1] + paneltry)
            score_matrix[one_index][two_index] = minmax_func(choose_list)
    return score_matrix


def get_best_score(score_matrix):
    return score_matrix[len(score_matrix) - 1][len(score_matrix[0]) - 1]


class IUPACAlignmentError(Exception):
    def __init__(self, message):
        self.message = message


# returns all optimal alignments
def generate_optimal_alignments(seq_one, seq_two, score_matrix=None, match_score=1, panelty_del_N=0, match_score_N=1,
                                delete_penalty=0, insert_penalty=0, mismatch_penalty=-1, return_single=False,
                                minmax_func=max):
    if score_matrix is None:
        score_matrix = align_iupac_dna_sequence(seq_one, seq_two, match_score=match_score, match_score_N=match_score_N,
                                                delete_penalty=delete_penalty, insert_penalty=insert_penalty,
                                                mismatch_penalty=mismatch_penalty, panelty_del_N=panelty_del_N,
                                                minmax_func=minmax_func)
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
                temp_match_score = match_score_N if seq_two[matrix_index[1] - 1] == 'N' else match_score
                diagonal_score = score_matrix[matrix_index[0] - 1][matrix_index[1] - 1]
                # Match
                if curr_score == diagonal_score + agree(seq_one[matrix_index[0] - 1], seq_two[matrix_index[1] - 1],
                                                        match_score_N, match_score, mismatch_penalty):
                    round_align_one = seq_one[matrix_index[0] - 1] + align_one
                    round_align_two = seq_two[matrix_index[1] - 1] + align_two
                    matched_this_round = True
                    round_matrix_index = (matrix_index[0] - 1, matrix_index[1] - 1)
            # insert
            if matrix_index[0] > 0 and curr_score == score_matrix[matrix_index[0] - 1][matrix_index[1]] + \
                    insert_penalty:
                if matched_this_round and not return_single:
                    alignment_stack.insert(0, (seq_one[matrix_index[0] - 1] + align_one, '-' + align_two,
                                               (matrix_index[0] - 1, matrix_index[1])))
                elif not matched_this_round:
                    round_align_one = seq_one[matrix_index[0] - 1] + align_one
                    round_align_two = '-' + align_two
                    round_matrix_index = (matrix_index[0] - 1, matrix_index[1])
                    matched_this_round = True
            # delete - special case for deleting N
            if matrix_index[1] > 0:
                temp_del_penalty = panelty_del_N if seq_two[matrix_index[1] - 1] == 'N' else delete_penalty
                if curr_score == score_matrix[matrix_index[0]][matrix_index[1] - 1] + temp_del_penalty:
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
                raise IUPACAlignmentError("Score matrix, no match! {}\n{}".format((matrix_index[0], matrix_index[1] - 1)
                                                                                  , score_matrix))
            align_one = round_align_one
            align_two = round_align_two
            matrix_index = round_matrix_index
        # print("End stack:\n{}\n{}".format(align_one, align_two))
        results.append((align_one, align_two))
    return results


if __name__ == "__main__":
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
    seq_a = "AGGUAGGCAUCGCGCGTTCGUACTAGTCGATCAUTGCATCGACGGGATUCGAAGCA"  # AGGUAGGCA---UCGCGCG-TT-C-GUAC-TAG--TCGATC--AUTGC-ATCGACGGGATUCGAAGCA-
    seq_b = "NNNNNNNRRRRNNNNNGGUUGCCUACGUAGGGCUCCCGCCNNNNNNNKKHHHHNNNCAU"  # NNNNNNN-RRRRNNNNN-GGUUGCC-UACGUAGGG-C--UCCC---GCCNNNNNNNKKHHHHNNN-CAU
    single_res = True

    print("starting alignment:\n{}\n{}".format(seq_a, seq_b))
    start_time = time.clock()
    my_score_matrix = align_iupac_dna_sequence(seq_a, seq_b, delete_penalty=-50, mismatch_penalty=-100, match_score=0)
    print("Score matrix took {} seconds:\n{}".format((time.clock() - start_time) / 1000, my_score_matrix))
    start_time = time.clock()
    alignment_list = generate_optimal_alignments(seq_a, seq_b, score_matrix=my_score_matrix, return_single=single_res
                                                 , delete_penalty=-50, mismatch_penalty=-100, match_score=0)
    print("Generating {} optimal alignment{}: took {} seconds".format('single' if single_res else 'all',
                                                                      '' if single_res else 's',
                                                                      (time.clock() - start_time) / 1000))
    for index, alignment_str in zip(range(0, len(alignment_list)), alignment_list):
        print("alignment {}:\n{}\n{}".format(index, alignment_str[0], alignment_str[1]))
    # Test 2
    start_time = time.clock()
    my_score_matrix = align_iupac_dna_sequence(seq_a, seq_b, delete_penalty=-100, mismatch_penalty=-100,
                                               panelty_del_N=-1, match_score=0, insert_penalty=-1)
    print("Score matrix took {} seconds:\n{}".format((time.clock() - start_time) / 1000, my_score_matrix))
    start_time = time.clock()
    alignment_list = generate_optimal_alignments(seq_a, seq_b, score_matrix=my_score_matrix, return_single=single_res
                                                 , delete_penalty=-100, mismatch_penalty=-100,
                                                 panelty_del_N=-1, match_score=0, insert_penalty=-1)
    print("Generating {} optimal alignment{}: took {} seconds".format('single' if single_res else 'all',
                                                                      '' if single_res else 's',
                                                                      (time.clock() - start_time) / 1000))
    for index, alignment_str in zip(range(0, len(alignment_list)), alignment_list):
        print("alignment {}:\n{}\n{}".format(index, alignment_str[0], alignment_str[1]))
    # Test 3
    start_time = time.clock()
    my_score_matrix = align_iupac_dna_sequence(seq_a, seq_b, delete_penalty=0, mismatch_penalty=0, match_score_N=1,
                                               panelty_del_N=0, match_score=100, insert_penalty=0)
    print("Score matrix took {} seconds:\n{}".format((time.clock() - start_time) / 1000, my_score_matrix))
    start_time = time.clock()
    alignment_list = generate_optimal_alignments(seq_a, seq_b, score_matrix=my_score_matrix, return_single=single_res,
                                                 delete_penalty=0, mismatch_penalty=0, match_score_N=1,
                                                 panelty_del_N=0, match_score=100, insert_penalty=0)
    print("Generating {} optimal alignment{}: took {} seconds".format('single' if single_res else 'all',
                                                                      '' if single_res else 's',
                                                                      (time.clock() - start_time) / 1000))
    for index, alignment_str in zip(range(0, len(alignment_list)), alignment_list):
        print("alignment {}:\n{}\n{}".format(index, alignment_str[0], alignment_str[1]))

    # Test 4
    start_time = time.clock()
    my_score_matrix = align_iupac_dna_sequence(seq_a, seq_b, delete_penalty=100, mismatch_penalty=100, match_score_N=0,
                                               panelty_del_N=1, match_score=0, insert_penalty=1, minmax_func=min)
    print("Score matrix took {} seconds:\n{}".format((time.clock() - start_time) / 1000, my_score_matrix))
    start_time = time.clock()
    alignment_list = generate_optimal_alignments(seq_a, seq_b, score_matrix=my_score_matrix, return_single=single_res,
                                                 delete_penalty=100, mismatch_penalty=100, match_score_N=0,
                                                 panelty_del_N=1, match_score=0, insert_penalty=1, minmax_func=min)
    print("Generating {} optimal alignment{}: took {} seconds".format('single' if single_res else 'all',
                                                                      '' if single_res else 's',
                                                                      (time.clock() - start_time) / 1000))
    for index, alignment_str in zip(range(0, len(alignment_list)), alignment_list):
        print("alignment {}:\n{}\n{}".format(index, alignment_str[0], alignment_str[1]))
