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
from typing import List, Tuple
import itertools

MIN_VALUE = -sys.maxsize - 1


def get_align_tree_distance(aligned_tree):
    node_stack = list()
    node_stack.append(aligned_tree)
    unmatched_counter = 0
    while len(node_stack) > 0:
        node = node_stack.pop()
        if node.mode != 'M':
            unmatched_counter += 1
        for child in node.children:
            node_stack.append(child)
    return unmatched_counter


class Tree:
    # 3 mode, M for match, S for source, T for target. exist for aligned trees.
    def __init__(self, value, children, mode='M', index=None, alignment=None):
        self.value = value
        self.children = children
        self.mode = mode
        self.index = index
        self.value_align = alignment

    def add_child(self, child):
        self.children.append(child)

    def add_children(self, children):
        self.children.extend(children)

    def full_str(self):
        res = "{}{} {}{} [".format('{})'.format(self.index) if self.index is not None else '', self.value, self.mode,
                                   ':{}'.format(self.value_align) if self.value_align is not None else '')
        for child in self.children:
            res += child.full_str()
            if child != self.children[-1]:
                res += ", "
        res += "]"
        return res
        # return "{} {}".format(self.value, [str(child) for child in self.children])

    def __str__(self):
        res = ''
        if self.mode == 'M':
            res = '{}{} ['.format(self.value, ':{}'.format(self.value_align) if self.value_align is not None else '')
            for child in self.children:
                child_str = str(child)
                if child_str != '':
                    res += child_str + ", "
            if res[-2:] == ', ':
                res = res[:-2]
            res += ']'
        else:
            for child in self.children:
                child_str = str(child)
                if child_str != '':
                    res += child_str + ", "
            if res[-2:] == ', ':
                res = res[:-2]
        return res

    def __repr__(self):
        return str(self)


def def_cmp(x, y):
    if int(x) == int(y):
        return 10, None
    return None, None


def def_merge(x, y):
    return x


def def_delete_func(x, is_target=False):
    return -1 * (2 if is_target else 1), None


# object with alignment rules.
# cmp_func : a function that receives 2 tree values and returns the comparison score
# merge_func : a function that receives 2 tree values and returns a merged consensus
# delete_func : a function that receives a single tree value and a boolean marking false as source tree and True as target
class AlignmentRules:
    def __init__(self, minmax_func=max, delete_func=def_delete_func, cmp_func=def_cmp,
                 merge_func=def_merge):
        self.cmp_func = cmp_func
        self.merge_func = merge_func
        self.minmax_func = minmax_func
        self.delete_func = delete_func


'''
# iterative tree alignment ERROR IN VERSION!!! not taking into account tree hierarchy!
def align_trees(tree_one, tree_two, alignment_object):
    def index_tree(tree):
        node_stack = [tree]
        node_list = []
        while len(node_stack) > 0:
            node = node_stack.pop()
            node_list.append(node)
            for child in node.children[::-1]:
                node_stack.append(child)
        index = 0
        node_list = node_list[::-1]
        for node in node_list:
            node.index = index
            index += 1
        return node_list

    def merge_trees(node_list_one, node_list_two):
        def find_children_index(index_list, node_index_list):
            children_indexs = set()
            for index in index_list:
                found = False
                for i in range(0, len(node_index_list)):
                    x = node_index_list[i]
                    if index == x:
                        children_indexs.add(i)
                        break
            return children_indexs

        # backtrack node list
        index_one = len(node_list_one)
        index_two = len(node_list_two)
        merged_node_list = []
        pair_list = []
        while index_one > 0 or index_two > 0:
            merged_node_list.append(transition_matrix[index_one][index_two][1])
            if index_one > 0 and index_two > 0 and transition_matrix[index_one][index_two][0] + \
                    score_matrix[index_one - 1][index_two - 1] == score_matrix[index_one][index_two]:
                index_one -= 1
                index_two -= 1
            elif index_one > 0 and transition_matrix[index_one][index_two][0] + \
                    score_matrix[index_one - 1][index_two] == score_matrix[index_one][index_two]:
                index_one -= 1
            elif index_two > 0 and transition_matrix[index_one][index_two][0] + \
                    score_matrix[index_one][index_two - 1] == score_matrix[index_one][index_two]:
                index_two -= 1
            else:
                raise Exception('Merge trees backtrack error while constructing list. no legal transition\n'
                                'Indexs[{}][{}]\nscore matrix: {}\n transition matrix: {}'
                                .format(index_one, index_two, score_matrix, transition_matrix))
            pair_list.append((index_one, index_two))
        # reconstruct as tree
        seen_children = set()
        for node, (index_one, index_two) in zip(merged_node_list[::-1], pair_list[::-1]):
            if node.mode == 'M':
                children = find_children_index([child.index for child in node_list_one[index_one].children],
                                               [x for x, y in pair_list])
                children = children.union(
                    find_children_index([child.index for child in node_list_two[index_two].children],
                                        [y for x, y in pair_list]))
            elif node.mode == 'S':
                children = find_children_index([child.index for child in node_list_one[index_one].children],
                                               [x for x, y in pair_list])
            elif node.mode == 'T':
                children = find_children_index([child.index for child in node_list_two[index_two].children],
                                               [y for x, y in pair_list])

            children = children.difference(seen_children)
            seen_children = seen_children.union(children)
            node.add_children([merged_node_list[child_index] for child_index in children])
        return merged_node_list[0]

    # index trees (lower is further) and create matrix
    tree_one_indexed = index_tree(tree_one)
    tree_two_indexed = index_tree(tree_two)
    score_matrix = [([0] * (len(tree_two_indexed) + 1)) for i in range(0, len(tree_one_indexed) + 1)]
    transition_matrix = [([0] * (len(tree_two_indexed) + 1)) for i in range(0, len(tree_one_indexed) + 1)]
    #
    for index_one in range(0, len(tree_one_indexed) + 1):
        for index_two in range(0, len(tree_two_indexed) + 1):
            if index_one == 0 and index_two == 0:
                continue
            score_list = []
            transition_list = []
            # match
            if index_one > 0 and index_two > 0:
                cmp_value, cmp_align = alignment_object.cmp_func(tree_one_indexed[index_one - 1].value,
                                                                 tree_two_indexed[index_two - 1].value)
                if cmp_value is not None:
                    transition_value = cmp_value
                    score_list.append(transition_value + score_matrix[index_one - 1][index_two - 1])
                    transition_list.append((transition_value,
                                            Tree(alignment_object.merge_func(tree_one_indexed[index_one - 1].value,
                                                                             tree_two_indexed[index_two - 1].value),
                                                 [], mode='M', alignment=cmp_align)))
            # insert from tree one (ignore node)
            if index_one > 0:
                transition_value, transition_align = alignment_object.delete_func(tree_one_indexed[index_one - 1].value,
                                                                                  is_target=False)
                score_list.append(transition_value + score_matrix[index_one - 1][index_two])
                transition_list.append((transition_value,
                                        Tree(tree_one_indexed[index_one - 1].value, [],
                                             mode='S', alignment=transition_align)))
            # insert from tree two (ignore node)
            if index_two > 0:
                transition_value, transition_align = alignment_object.delete_func(tree_two_indexed[index_two - 1].value,
                                                                                  is_target=True)
                score_list.append(transition_value + score_matrix[index_one][index_two - 1])
                transition_list.append((transition_value,
                                        Tree(tree_two_indexed[index_two - 1].value, [],
                                             mode='T', alignment=transition_align)))
            minmax_value = alignment_object.minmax_func(score_list)
            score_matrix[index_one][index_two] = minmax_value
            transition_matrix[index_one][index_two] = transition_list[score_list.index(minmax_value)]
    return merge_trees(tree_one_indexed, tree_two_indexed), score_matrix[index_one][index_two]
'''


# iterative tree alignment
def align_trees(tree_one: Tree, tree_two: Tree, alignment_object: AlignmentRules) -> Tuple[int, Tree]:
    def index_tree(tree: Tree) -> List[Tree]:
        node_stack = [tree]
        node_list = []
        while len(node_stack) > 0:
            node = node_stack.pop()
            node_list.append(node)
            for index_child in node.children[::-1]:
                node_stack.append(index_child)
        index = 0
        node_list = node_list[::-1]
        for node in node_list:
            node.index = index
            index += 1
        return node_list

    def del_subtrees(children_list: List[Tree], non_del_child: Tree, is_target: bool) -> List[Tuple[int, Tree]]:
        del_tree_list = []
        for del_child in children_list:
            if non_del_child is not None and del_child == non_del_child:
                continue
            del_transition_node = None
            tree_stack: List[Tuple[Tree, Tree]] = [(del_child, None)]
            cost = 0
            while tree_stack:
                curr_node, parent_node = tree_stack.pop()
                del_value, del_align = alignment_object.delete_func(curr_node.value, is_target=is_target)
                cost += del_value
                del_tree = Tree(curr_node.value, [], mode='T' if is_target else 'S',
                                alignment=del_align)
                if del_transition_node is None:
                    del_transition_node = del_tree
                else:
                    parent_node.children.append(del_tree)
                for curr_child in curr_node.children[::-1]:
                    tree_stack.append((curr_child, del_tree))
            del_tree_list.append((cost, del_transition_node))
        return del_tree_list

    def all_combination(source_children: List[Tree], target_children: List[Tree]) -> List[Tuple[int, Tree]]:
        del_source = del_subtrees(source_children, None, False)
        del_target = del_subtrees(target_children, None, True)
        # generate a list of all combinations with order
        all_options = [[]]
        short_list = source_children if len(source_children) <= len(target_children) else target_children
        long_list = source_children if len(source_children) > len(target_children) else target_children
        for locations_long in itertools.combinations(range(len(long_list)), len(short_list)):
            for locations_short in itertools.combinations(range(len(short_list)), len(short_list)):
                m_list = [item for i, item in zip(range(len(short_list)), short_list) if i in locations_short]
                for i, j in zip(locations_long, range(len(short_list))):
                    m_list[j] = (m_list[j], long_list[i]) if short_list == source_children else \
                        (long_list[i], m_list[j])
                all_options.append(m_list)
        best_score_option = -alignment_object.minmax_func(sys.maxsize, -sys.maxsize)
        best_child_list = []
        for single_option in all_options:
            last_added_source = 0
            last_added_target = 0
            curr_score = 0
            curr_children_list = []
            for match_pair in single_option:
                for source_index in range(last_added_source, source_children.index(match_pair[0])):
                    curr_score += del_source[source_index][0]
                    curr_children_list.append(del_source[source_index][1])
                for target_index in range(last_added_target, target_children.index(match_pair[1])):
                    curr_score += del_target[target_index][0]
                    curr_children_list.append(del_target[target_index][1])
                curr_score += score_matrix[match_pair[0].index][match_pair[1].index]
                curr_children_list.append(transition_matrix[match_pair[0].index][match_pair[1].index])
                last_added_source = source_children.index(match_pair[0]) + 1
                last_added_target = target_children.index(match_pair[1]) + 1
            for source_index in range(last_added_source, len(del_source)):
                curr_score += del_source[source_index][0]
                curr_children_list.append(del_source[source_index][1])
            for target_index in range(last_added_target, len(del_target)):
                curr_score += del_target[target_index][0]
                curr_children_list.append(del_target[target_index][1])
            if alignment_object.minmax_func(curr_score, best_score_option) == curr_score:
                best_score_option = curr_score
                best_child_list = curr_children_list
        return best_score_option, best_child_list

    # index trees (lower is further) and create matrix
    tree_one_indexed = index_tree(tree_one)
    tree_two_indexed = index_tree(tree_two)
    score_matrix = [([-alignment_object.minmax_func(sys.maxsize, -sys.maxsize)] * (len(tree_two_indexed)))
                    for i in range(len(tree_one_indexed))]
    transition_matrix = [([None] * (len(tree_two_indexed))) for i in range(len(tree_one_indexed))]
    #
    for index_one in range(len(tree_one_indexed)):
        for index_two in range(len(tree_two_indexed)):
            score_list = []
            transition_list = []
            # match, get best children combination from both
            cmp_value, cmp_align = alignment_object.cmp_func(tree_one_indexed[index_one].value,
                                                             tree_two_indexed[index_two].value)
            if cmp_value is not None:
                transition_value = cmp_value
                best_match_score, child_list = all_combination(tree_one_indexed[index_one].children,
                                                               tree_two_indexed[index_two].children)
                score_list.append(transition_value + best_match_score)
                transition_list.append(Tree(alignment_object.merge_func(tree_one_indexed[index_one].value,
                                                                        tree_two_indexed[index_two].value),
                                            child_list, mode='M', alignment=cmp_align))
            # insert from tree one (ignore node), get best child to index two compare (del rest of children)
            transition_value, transition_align = alignment_object.delete_func(tree_one_indexed[index_one].value,
                                                                              is_target=False)
            opt_t1_child_t2, opt_t1_children = all_combination(tree_one_indexed[index_one].children,
                                                               [tree_two_indexed[index_two]])
            score_list.append(transition_value + opt_t1_child_t2)
            transition_list.append(Tree(tree_one_indexed[index_one].value, opt_t1_children,
                                        mode='S', alignment=transition_align))
            # insert from tree two (ignore node), get best child to index one compare (del rest of children)
            transition_value, transition_align = alignment_object.delete_func(tree_two_indexed[index_two].value,
                                                                              is_target=True)
            opt_t1_t2_child, opt_t2_children = all_combination([tree_one_indexed[index_one]],
                                                               tree_two_indexed[index_two].children)
            score_list.append(transition_value + opt_t1_t2_child)
            transition_list.append(Tree(tree_two_indexed[index_two].value, opt_t2_children,
                                        mode='T', alignment=transition_align))
            # Check best action
            minmax_value = alignment_object.minmax_func(score_list)
            score_matrix[index_one][index_two] = minmax_value
            transition_matrix[index_one][index_two] = transition_list[score_list.index(minmax_value)]
    return transition_matrix[index_one][index_two], score_matrix[index_one][index_two]


if __name__ == "__main__":
    ar = AlignmentRules()
    test_tree_one = Tree(1.1, [Tree(2.1, [Tree(3.1, [Tree(4.1, [])]), Tree(5.1, [Tree(6.1, [])])]),
                               Tree(2.1, [Tree(3.1, [])])])
    test_tree_two = Tree(1.2, [Tree(2.2, [Tree(4.2, [])]), Tree(2.2, [Tree(4.2, [Tree(3.2, [])])])])
    print("Tree one: {}".format(test_tree_one))
    print("Tree two: {}".format(test_tree_two))
    align_tree, align_score = align_trees(test_tree_one, test_tree_two, ar)
    print("Matched only Joined Tree ({}): {}".format(align_score, align_tree))
    print("Full Joined Tree ({}): {}".format(align_score, align_tree.full_str()))
