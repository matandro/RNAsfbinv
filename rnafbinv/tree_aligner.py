'''
Dynamic programing recursive alignment of two trees. Tree node defined with value and list of children.
requires an AlignmentRules object which includes functions for:
    1) comparison of two values (return Palely / reward)
    2) merging two tree values (for the merge tree)
    3) Penalty for insertion from tree two
    4) Penalty for deletion from tree one
'''

import sys
from typing import List, Tuple, Callable, TypeVar, Generic


MIN_VALUE = -sys.maxsize - 1
TreeValue = TypeVar('TreeValue')
AlignmentResult = TypeVar('AlignmentResult')


class Tree(Generic[TreeValue]):
    # 3 mode, M for match, S for source, T for target. exist for aligned trees.
    def __init__(self, value: TreeValue, children: List['Tree'], mode='M', index=None, alignment=None):
        self.value = value
        self.children = children
        self.mode = mode
        self.index = index
        self.value_align = alignment

    def add_child(self, child: 'Tree'):
        self.children.append(child)

    def add_children(self, children: List['Tree']):
        self.children.extend(children)

    def full_str(self) -> str:
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


def get_align_tree_distance(aligned_tree: Tree) -> str:
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


def def_cmp(x: TreeValue, y: TreeValue) -> Tuple[float, AlignmentResult]:
    if int(x) == int(y):
        return 10, None
    return None, None


def def_merge(x: TreeValue, y: TreeValue) -> TreeValue:
    return x


def def_delete_func(x: TreeValue, is_target: bool = False) -> Tuple[float, AlignmentResult]:
    return -1 * (2 if is_target else 1), None


# object with alignment rules.
# T - tree value type
# G - alignment result type
# cmp_func : a function that receives 2 tree values and returns the comparison score
# merge_func : a function that receives 2 tree values and returns a merged consensus
# delete_func : a function that receives a single tree value and a boolean marking false as source tree and True as target
class AlignmentRules(Generic[TreeValue]):
    def __init__(self, minmax_func: Callable[[float, float], float] = max,
                 delete_func: Callable[[TreeValue, bool], Tuple[float, AlignmentResult]] = def_delete_func,
                 cmp_func: Callable[[TreeValue, TreeValue], Tuple[float, AlignmentResult]] = def_cmp,
                 merge_func: Callable[[TreeValue, TreeValue], TreeValue] = def_merge):
        self.cmp_func = cmp_func
        self.merge_func = merge_func
        self.minmax_func = minmax_func
        self.delete_func = delete_func


# iterative tree alignment
def align_trees(tree_one: Tree, tree_two: Tree, alignment_object: AlignmentRules) -> Tuple[Tree, float]:
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
            # tree_stack: List[Tuple[Tree, Tree]] = [(del_child, None)]
            tree_stack = [(del_child, None)]
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

    def all_combination(source_children: List[Tree], target_children: List[Tree]) -> Tuple[int, List[Tree]]:
        del_source = del_subtrees(source_children, None, False)
        del_target = del_subtrees(target_children, None, True)
        # handle an empty tree list
        if not source_children or not target_children:
            if source_children:
                return sum([single_del[0] for single_del in del_source]), [single_del[1] for single_del in del_source]
            elif target_children:
                return sum([single_del[0] for single_del in del_target]), [single_del[1] for single_del in del_target]
            else:
                return 0, []
        # both have at least one child
        children_scores = {}
        children_trees = {}
        for i in range(len(source_children), -1, -1):
            for k in range(len(target_children), -1, -1):
                for j in range(len(source_children), -1, -1):
                    for l in range(len(target_children), -1, -1):
                        # handle tree list vs empty
                        if i >= j or k >= l:
                            if i < j:
                                best_score = sum([del_source[m][0] for m in range(i, j)])
                                best_tree = [del_source[m][1] for m in range(i, j)]
                            elif k < l:
                                best_score = sum([del_target[m][0] for m in range(k, l)])
                                best_tree = [del_target[m][1] for m in range(k, l)]
                            else:
                                best_score = 0
                                best_tree = []
                            children_scores[(i, j, k, l)] = best_score
                            children_trees[(i, j, k, l)] = best_tree
                            continue
                        # handle tree list vs tree list
                        score_options = []
                        tree_options = []
                        # Align first two
                        temp_score = score_matrix[source_children[i].index][target_children[k].index] + \
                                     children_scores[(i + 1, j, k + 1, l)]
                        temp_subtree = [transition_matrix[source_children[i].index][target_children[k].index]] + \
                                       children_trees[(i + 1, j, k + 1, l)]
                        score_options.append(temp_score)
                        tree_options.append(temp_subtree)
                        # Insert first from source
                        temp_score = -alignment_object.minmax_func(sys.maxsize, -sys.maxsize)
                        temp_subtree = []
                        found = False
                        for m in range(k - 1, l):
                            s, tl = all_combination(source_children[i].children, target_children[k:m])
                            m_score = s + children_scores[(i + 1, j, m + 1, l)]
                            if alignment_object.minmax_func(m_score, temp_score) == m_score:
                                found = True
                                temp_score = m_score
                                temp_subtree = tl + children_trees[(i + 1, j, m + 1, l)]
                        if not found:
                            temp_score = 0
                        temp_score += del_source[i][0]
                        temp_subtree = [del_source[i][1]] + temp_subtree
                        score_options.append(temp_score)
                        tree_options.append(temp_subtree)
                        # Del first from target
                        temp_score = -alignment_object.minmax_func(sys.maxsize, -sys.maxsize)
                        temp_subtree = []
                        found = False
                        for m in range(i - 1, j):
                            s, tl = all_combination(source_children[i:m], target_children[k].children)
                            m_score = s + children_scores[(m + 1, j, k + 1, l)]
                            if alignment_object.minmax_func(m_score, temp_score) == m_score:
                                found = False
                                temp_score = m_score
                                temp_subtree = tl + children_trees[(m + 1, j, k + 1, l)]
                        if not found:
                            temp_score = 0
                        temp_score += del_target[k][0]
                        temp_subtree = [del_target[k][1]] + temp_subtree
                        score_options.append(temp_score)
                        tree_options.append(temp_subtree)
                        # find best
                        if score_options:
                            best_score = alignment_object.minmax_func(score_options)
                            best_tree = tree_options[score_options.index(best_score)]
                        else:
                            best_score, best_tree = 0, []
                        children_scores[(i, j, k, l)] = best_score
                        children_trees[(i, j, k, l)] = best_tree

        return (children_scores[(0, len(source_children), 0, len(target_children))],
                children_trees[(0, len(source_children), 0, len(target_children))])

    def min_children(single: Tree, children: List[Tree], is_child_first: bool) -> Tuple[int, List[Tree]]:
        child_scores = [score_matrix[one_child.index][single.index] for one_child in children] if is_child_first else \
            [score_matrix[single.index][one_child.index] for one_child in children]
        if child_scores:
            child_trees = [transition_matrix[one_child.index][single.index] for one_child in children] if \
                is_child_first else [transition_matrix[single.index][one_child.index] for one_child in children]
            score_res = alignment_object.minmax_func(child_scores)
            tree_res = [child_trees[child_scores.index(score_res)]]
        else:
            score_res = 0
            tree_res = []
        return score_res, tree_res

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
            opt_t1_child_t2, opt_t1_children = min_children(tree_two_indexed[index_two],
                                                            tree_one_indexed[index_one].children, True)
            score_list.append(transition_value + opt_t1_child_t2)
            transition_list.append(Tree(tree_one_indexed[index_one].value, opt_t1_children,
                                        mode='S', alignment=transition_align))
            # insert from tree two (ignore node), get best child to index one compare (del rest of children)
            transition_value, transition_align = alignment_object.delete_func(tree_two_indexed[index_two].value,
                                                                              is_target=True)
            opt_t1_t2_child, opt_t2_children = min_children(tree_one_indexed[index_one],
                                                            tree_two_indexed[index_two].children, False)
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
    '''
    test_tree_one = Tree(1.1, [Tree(2.1, [Tree(3.1, [Tree(4.1, [])]), Tree(5.1, [Tree(6.1, [])])]),
                               Tree(2.1, [Tree(3.1, [])])])
    test_tree_two = Tree(1.2, [Tree(2.2, [Tree(4.2, [])]), Tree(2.2, [Tree(4.2, [Tree(3.2, [])])])])
    print("Tree one: {}".format(test_tree_one))
    print("Tree two: {}".format(test_tree_two))
    align_tree, align_score = align_trees(test_tree_one, test_tree_two, ar)
    print("Matched only Joined Tree ({}): {}".format(align_score, align_tree))
    print("Full Joined Tree ({}): {}".format(align_score, align_tree.full_str()))
    '''
    test_tree_one = Tree(1.1, [Tree(2.1, [Tree(3.1, []),
                                          Tree(4.1, [])]),
                               Tree(5.1, [])])
    test_tree_two = Tree(1.2, [Tree(3.2, []),
                               Tree(4.2, []),
                               Tree(5.2, []), ])
    print("Tree one: {}".format(test_tree_one))
    print("Tree two: {}".format(test_tree_two))
    align_tree, align_score = align_trees(test_tree_one, test_tree_two, ar)
    print("Matched only Joined Tree ({}): {}".format(align_score, align_tree))
    print("Full Joined Tree ({}): {}".format(align_score, align_tree.full_str()))
