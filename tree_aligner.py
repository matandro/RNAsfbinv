'''
Dynamic programing recursive alignment of two trees. Tree node defined with value and list of children.
requires an AlignmentRules object which includes functions for:
    1) comparison of two values (return Palely / reward)
    2) merging two tree values (for the merge tree)
    3) Penalty for insertion from tree two
    4) Penalty for deletion from tree one
'''

import logging
import sys
import copy

MIN_VALUE = -sys.maxsize - 1


class Tree:
    # 3 mode, M for match, S for source, T for target. exist for aligned trees.
    def __init__(self, value, children, mode='M'):
        self.value = value
        self.children = children
        self.mode = mode

    def add_child(self, child):
        self.children.append(child)

    def add_children(self, children):
        self.children.extend(children)

    def print_full(self):
        res = "{} M{} [".format(self.value, self.mode)
        for child in self.children:
            res += child.print_matched()
            if child != self.children[-1]:
                res += ", "
        res += "]"
        return res
        # return "{} {}".format(self.value, [str(child) for child in self.children])

    def __str__(self):
        res = ''
        if self.mode == 'M':
            res = '{} ['.format(self.value)
            for child in self.children:
                child_str = str(child)
                if child_str != '':
                    res += child_str + ", "
            if res[-2:] == ', ':
                res = res[:-2]
            res += ']'
        return res

    def __repr__(self):
        return str(self)


def def_cmp(x, y):
    if int(x) == int(y):
        return 10
    return None


def def_merge(x, y):
    return x


def def_insert_func(x):
    return -1


def def_delete_func(x, is_target=False):
    return -1 * (2 if is_target else 1)


# TODO: Add fucntion to say if we are looking for minimum or maximum (score comparison)
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
    def delete_brench(tree, is_target):
        del_score = 0
        cpy_tree = copy.deepcopy(tree)
        tree_stack = [cpy_tree]
        while len(tree_stack) > 0:
            current_tree = tree_stack.pop()
            for child in current_tree.children:
                tree_stack.append(child)
            current_tree.mode = 'T' if is_target else 'S'
            del_score += alignment_object.delete_func(current_tree.value, is_target)
        # logging.debug('delete brench {} with score {} marked tree {}'.format(tree, del_score, cpy_tree))
        return del_score, cpy_tree

    def retrace(score_matrix, transition_matrix):
        recall_subtree = []
        index_one = len(score_matrix) - 1
        index_two = len(score_matrix[0]) - 1
        while index_one > 0 or index_two > 0:
            transition_pair = transition_matrix[index_one][index_two]
            # match
            if index_one > 0 and index_two > 0 and score_matrix[index_one][index_two] == \
                    score_matrix[index_one - 1][index_two - 1] + transition_pair[1]:
                index_one -= 1
                index_two -= 1
            # insert
            elif index_one > 0 and score_matrix[index_one][index_two] == \
                    score_matrix[index_one - 1][index_two] + transition_pair[1]:
                index_one -= 1
            # delete
            elif index_two > 0 and score_matrix[index_one][index_two] == \
                    score_matrix[index_one][index_two - 1] + transition_pair[1]:
                index_two -= 1
            else:
                raise Exception("Score matrix, no match!\n{} @ [{}][{}] {}\n{}"
                                .format(score_matrix, one_index, two_index, transition_pair, transition_matrix))
            recall_subtree.append(transition_pair[0])
        return recall_subtree

    score_matrix = [([0] * (len(tree_two_children) + 1)) for i in range(0, (len(tree_one_children) + 1))]
    transition_matrix = [([0] * (len(tree_two_children) + 1)) for i in range(0, (len(tree_one_children) + 1))]
    for one_index in range(0, len(tree_one_children) + 1):
        for two_index in range(0, len(tree_two_children) + 1):
            alignment_list = []
            transition_list = []
            if one_index == 0 and two_index == 0:
                continue
            # match
            if one_index > 0 and two_index > 0:
                subtree, score = align_trees(tree_one_children[one_index - 1], tree_two_children[two_index - 1],
                                             alignment_object)
                transition_list.append((subtree, score))
                alignment_list.append(score_matrix[one_index - 1][two_index - 1] + score)
            # insert (not using this node from tree one)
            if one_index > 0:
                score, sub_tree = delete_brench(tree_one_children[one_index - 1], is_target=False)
                transition_list.append((sub_tree, score))
                alignment_list.append(score_matrix[one_index - 1][two_index] + score)
            # delete (not using this node from tree two)
            if two_index > 0:
                score, sub_tree = delete_brench(tree_two_children[two_index - 1], is_target=True)
                transition_list.append((sub_tree, score))
                alignment_list.append(score_matrix[one_index][two_index - 1] + score)
            max_score = max(alignment_list)
            transition_matrix[one_index][two_index] = transition_list[alignment_list.index(max_score)]
            score_matrix[one_index][two_index] = max_score
    final_subtree = retrace(score_matrix, transition_matrix)
    final_score = score_matrix[one_index][two_index]
    return final_subtree, final_score


# recursive tree alignment
def align_trees(tree_one, tree_two, alignment_object):
    # match tree one to tree two (None = uncomparable)
    cmp_res = alignment_object.cmp_func(tree_one.value, tree_two.value)
    if cmp_res is not None:
        match_child_list, match_score = compare_child_combinations(tree_one.children, tree_two.children,
                                                                   alignment_object)
        match_score += cmp_res
        match_subtree = Tree(alignment_object.merge_func(tree_one.value, tree_two.value), [], mode='M')
        match_subtree.add_children(match_child_list)
    else:
        match_score = -sys.maxsize - 1

        # No match for current node in tree one.
    temp_subtree, tree_one_score = compare_child_combinations(tree_one.children, [tree_two],
                                                              alignment_object)
    tree_one_subtree = Tree(tree_one.value, temp_subtree, mode='S')
    tree_one_score += alignment_object.delete_func(tree_one.value, is_target=False)

    # No match for current node in tree two.
    temp_subtree, tree_two_score = compare_child_combinations([tree_one], tree_two.children,
                                                              alignment_object)
    tree_two_subtree = Tree(tree_two.value, temp_subtree, mode='T')
    tree_two_score += alignment_object.delete_func(tree_two.value, is_target=True)

    # Take optimal result
    score = max([match_score, tree_one_score, tree_two_score])
    if score == match_score:
        subtree = match_subtree
    elif score == tree_one_score:
        subtree = tree_one_subtree
    elif score == tree_two_score:
        subtree = tree_two_subtree
    else:
        raise Exception('Must match a case {} {} {}'.format(match_score, tree_one_score, tree_two_score))
    return subtree, score


if __name__ == "__main__":
    ar = AlignmentRules()
    test_tree_one = Tree(1.1, [Tree(2.1, [Tree(3.1, [Tree(4.1, [])])]), Tree(2.1, [])])
    test_tree_two = Tree(1.2, [Tree(2.2, []), Tree(2.2, [Tree(3.2, [Tree(3.2, [])])])])
    print("Tree one: {}".format(test_tree_one))
    print("Tree two: {}".format(test_tree_two))
    align_tree, align_score = align_trees(test_tree_one, test_tree_two, ar)
    print("Matched only Joined Tree ({}): {}".format(align_score, align_tree))
    print("Full Joined Tree ({}): {}".format(align_score, align_tree.print_matched()))
