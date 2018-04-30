import logging
import sys
import time

MIN_VALUE = -sys.maxsize - 1


class Tree:
    def __init__(self, value, children):
        self.value = value
        self.children = children

    def add_child(self, child):
        self.children.append(child)

    def add_children(self, children):
        self.children.extend(children)

    def __str__(self):
        res = "{} [".format(self.value)
        for child in self.children:
            res += str(child)
            if child != self.children[-1]:
                res += ", "
        res += "]"
        return res
        # return "{} {}".format(self.value, [str(child) for child in self.children])


def def_cmp(x, y):
    if x == y:
        return 1
    return 0


def def_merge(x, y):
    return x


def def_insert_func(x):
    return -0.5


def def_delete_func(x):
    return -0.5


# object with alignment rules.
# cmp_func a function that receives 2 values and returns the bonus for a match (0 if non)
class AlignmentRules:
    def __init__(self, insert_func=def_insert_func, delete_func=def_delete_func, cmp_func=def_cmp,
                 merge_func=def_merge):
        self.cmp_func = cmp_func
        self.merge_func = merge_func
        self.insert_func = insert_func
        self.delete_func = delete_func


# try comparing all children combination WITH ORDER (no cross compare possible) - dynamic programming
def compare_child_combinations(tree_one_children, tree_two_children, alignment_object):
    def get_mex_in_submatrix(max_index_row, max_index_col, matrix):
        res = 0
        if max_index_col > 0:
            for index_row in range(0, max_index_row):
                if res < matrix[index_row][max_index_col - 1]:
                    res = matrix[index_row][max_index_col - 1]
        return res

    if len(tree_one_children) == 0 or len(tree_two_children) == 0:
        return [], 0
    res_child_list = []
    # Generate alignment for all children combinations
    score_matrix = [[0] * len(tree_two_children) for i in range(0, len(tree_one_children))]
    subtree_matrix = [[None] * len(tree_two_children) for i in range(0, len(tree_one_children))]
    sum_matrix = [[0] * len(tree_two_children) for i in range(0, len(tree_one_children))]  # for regression
    for one_index in range(0, len(tree_one_children)):
        for two_index in range(0, len(tree_two_children)):
            subtree_matrix[one_index][two_index], score_matrix[one_index][two_index] = align_trees(
                tree_one_children[one_index], tree_two_children[two_index], alignment_object)
            sum_matrix[one_index][two_index] = score_matrix[one_index][two_index] + get_mex_in_submatrix(
                one_index, two_index, sum_matrix)
    # select largest combination (regression on sum matrix)
    col_index = len(tree_two_children) - 1
    temp_score = res_score = get_mex_in_submatrix(len(tree_one_children), len(tree_two_children), sum_matrix)
    last_taken = len(tree_one_children)
    while col_index >= 0:
        for one_index in range(0, last_taken):
            if temp_score == sum_matrix[one_index][col_index]:
                if score_matrix[one_index][col_index] > 0:
                    last_taken = one_index
                    temp_score -= score_matrix[one_index][col_index]
                    res_child_list.insert(0, subtree_matrix[one_index][col_index])
                break
        col_index -= 1
    return res_child_list, res_score


def best_of_children(tree_two, tree_one_children, alignment_object):
    res_subtree = None
    res_score = 0
    if tree_one_children is not None:
        for child in tree_one_children:
            temp_subtree, temp_score = align_trees(child, tree_two, alignment_object)
            if temp_score > res_score:
                res_subtree = temp_subtree
                res_score = temp_score
    return res_subtree, res_score


# recursive tree alignment
def align_trees(tree_one, tree_two, alignment_object):
    # match
    cmp_res = alignment_object.cmp_func(tree_one.value, tree_two.value)
    if cmp_res > 0:
        match_child_list, match_score = compare_child_combinations(tree_one.children, tree_two.children,
                                                                   alignment_object)
        match_score += cmp_res
        match_subtree = Tree(alignment_object.merge_func(tree_one.value, tree_two.value), [])
        match_subtree.add_children(match_child_list)
    else:
        match_score = MIN_VALUE
    # insert
    insert_subtree, insert_score = compare_child_combinations(tree_one.children, [tree_two],
                                                              alignment_object)
    insert_subtree = Tree(tree_one.value, insert_subtree)
    insert_score += alignment_object.insert_func(tree_one.value)
    # delete
    delete_subtree, delete_score = best_of_children(tree_two, tree_one.children, alignment_object)
    delete_score += alignment_object.delete_func(tree_one.value)
    # Take optimal result
    score = max([match_score, delete_score, insert_score])
    if score == match_score:
        subtree = match_subtree
    elif score == delete_score:
        subtree = delete_subtree
    else:
        subtree = insert_subtree
    return subtree, score


if __name__ == "__main__":
    ar = AlignmentRules()
    test_tree_one = Tree(1, [Tree(2, [Tree(3, [])]), Tree(2, [])])
    test_tree_two = Tree(1, [Tree(2, []), Tree(2, [Tree(3, [])])])
    print("Tree one: {}".format(test_tree_one))
    print("Tree two: {}".format(test_tree_two))
    align_tree, align_score = align_trees(test_tree_one, test_tree_two, ar)
    print("Joined Tree ({}): {}".format(align_score, align_tree))

    # questions:
    # 1. how much does the order matter?
    # match - values are equal get maximum children comparison
    # insert - compare tree_two children to my value take maximum
    # delete - compare my children to tree_two value. take maximum of children
