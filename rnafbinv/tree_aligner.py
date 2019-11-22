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

    def init_deletions():
        init_deletion_tree(tree_one, False)
        init_deletion_tree(tree_two, True)

    def init_deletion_tree(top: Tree, is_target: bool):
        del_dict = del_target if is_target else del_source
        tree_stack = [top]
        tree_list = []
        while tree_stack:
            curr_node = tree_stack.pop()
            del_value, del_align = alignment_object.delete_func(curr_node.value, is_target)
            tree_list.insert(0, (curr_node, del_value, del_align))
            for curr_child in curr_node.children[::-1]:
                tree_stack.append(curr_child)
        for node, del_value, del_align in tree_list:
            del_tree = Tree(node.value, [], mode='T' if is_target else 'S',
                            alignment=del_align)
            child_cost = sum(del_dict[child.index][0] for child in node.children) if node.children else 0
            child_tree = [del_dict[child.index][1] for child in node.children] if node.children else []
            del_tree.add_children(child_tree)
            del_dict[node.index] = (child_cost + del_value, del_tree)

    def get_child_dp_info(source_tree: Tree, i: int, j: int, target_tree: Tree, k: int, l: int) -> Tuple[int, List[Tree]]:
        score = children_scores.get((source_tree.index, i, j, target_tree.index, k, l))
        if score is None:
            if i >= j or k >= l:
                if i < j and source_tree.children[i:j]:
                    score = sum([del_source[child.index][0] for child in source_tree.children[i:j]])
                    tree_list = [del_source[child.index][1] for child in source_tree.children[i:j]]
                elif k < l and target_tree.children[k:l]:
                    score = sum([del_target[child.index][0] for child in target_tree.children[k:l]])
                    tree_list = [del_target[child.index][1] for child in target_tree.children[k:l]]
                else:
                    score = 0
                    tree_list = []
            else:
                raise ValueError(
                    f"{source_tree.index}[{i}:{j}] - {target_tree.index}[{k}:{l}] should have been initialized")
            children_scores[(source_tree.index, i, j, target_tree.index, k, l)] = score
            children_trees[(source_tree.index, i, j, target_tree.index, k, l)] = tree_list
        else:
            tree_list = children_trees[(source_tree.index, i, j, target_tree.index, k, l)]
        return score, tree_list

    def all_combination(source_parent: Tree, target_parent: Tree) -> Tuple[int, List[Tree]]:
        source_children = source_parent.children
        target_children = target_parent.children
        for i in range(len(source_children), -1, -1):
            for j in range(len(source_children), i, -1):
                for k in range(len(target_children), -1, -1):
                    for l in range(len(target_children), k, -1):
                        # handle tree list vs tree list
                        score_options = []
                        tree_options = []
                        # Align first two
                        add_score, add_tree = get_child_dp_info(source_parent, i + 1, j, target_parent, k + 1, l)
                        temp_score = score_matrix[source_children[i].index][target_children[k].index] + add_score
                        temp_subtree = [transition_matrix[source_children[i].index][target_children[k].index]] + \
                                       add_tree
                        score_options.append(temp_score)
                        tree_options.append(temp_subtree)
                        # Insert first from source
                        temp_score = -alignment_object.minmax_func(sys.maxsize, -sys.maxsize)
                        temp_subtree = []
                        temp_contree = None
                        found = False
                        for m in range(k, l):
                            s, tl = get_child_dp_info(source_children[i], 0, len(source_children[i].children),
                                                      target_parent, k, m)
                            add_score, add_tree = get_child_dp_info(source_parent, i + 1, j, target_parent, m, l)
                            m_score = s + add_score
                            if alignment_object.minmax_func(m_score, temp_score) == m_score:
                                found = True
                                temp_score = m_score
                                temp_subtree = tl
                                temp_contree = add_tree
                        if not found:
                            temp_score = del_source[source_children[i].index][0]
                            temp_subtree = [del_source[source_children[i].index][1]]
                        else:
                            only_main_score, temp_align = alignment_object.delete_func(source_children[i].value, False)
                            temp_score += only_main_score
                            temp_subtree = [Tree(source_children[i].value, temp_subtree, mode='S'
                                                 , alignment=temp_align)] + temp_contree
                        score_options.append(temp_score)
                        tree_options.append(temp_subtree)
                        # Del first from target
                        temp_score = -alignment_object.minmax_func(sys.maxsize, -sys.maxsize)
                        temp_subtree = []
                        temp_contree = None
                        found = False
                        for m in range(i, j):
                            s, tl = get_child_dp_info(source_parent, i, m,
                                                      target_children[k], 0, len(target_children[k].children))
                            add_score, add_tree = get_child_dp_info(source_parent, m, j, target_parent, k + 1, l)
                            m_score = s + add_score
                            if alignment_object.minmax_func(m_score, temp_score) == m_score:
                                found = False
                                temp_score = m_score
                                temp_subtree = tl
                                temp_contree = add_tree
                        if not found:
                            temp_score = del_target[target_children[k].index][0]
                            temp_subtree = [del_target[target_children[k].index][1]]
                        else:
                            only_main_score, temp_align = alignment_object.delete_func(target_children[k].value, True)
                            temp_score += only_main_score
                            temp_subtree = [Tree(target_children[k].value, temp_subtree, mode='T',
                                                 alignment=temp_align)] + temp_contree
                        score_options.append(temp_score)
                        tree_options.append(temp_subtree)
                        # find best
                        best_score = alignment_object.minmax_func(score_options)
                        best_tree = tree_options[score_options.index(best_score)]
                        children_scores[(source_parent.index, i, j, target_parent.index, k, l)] = best_score
                        children_trees[(source_parent.index, i, j, target_parent.index, k, l)] = best_tree
        return get_child_dp_info(source_parent, 0, len(source_children),
                                 target_parent, 0, len(target_children))

    def best_of_children(single: Tree, children: List[Tree], is_child_first: bool) -> Tuple[int, List[Tree]]:
        if not children:
            return del_target[single.index] if is_child_first else del_source[single.index]
        best_score = -alignment_object.minmax_func(sys.maxsize, -sys.maxsize)
        best_tree = []
        for i in range(len(children)):
            child_score = score_matrix[children[i].index][single.index] if is_child_first else \
                score_matrix[single.index][children[i].index]
            del_info = [del_source[child.index] for child in children] if is_child_first else \
                [del_target[child.index] for child in children]
            del_score = sum([del_info[j][0] for j in range(len(del_info)) if j != i])
            total_score = child_score + del_score
            if alignment_object.minmax_func(best_score, total_score) == total_score:
                child_tree = transition_matrix[children[i].index][single.index] if is_child_first else \
                    transition_matrix[single.index][children[i].index]
                best_score = total_score
                best_tree = [del_info[j][1] for j in range(i)] + [child_tree] + \
                            [del_info[j][1] for j in range(i + 1, len(del_info))]
        return best_score, best_tree

    # index trees (lower is further) and create matrix
    tree_one_indexed = index_tree(tree_one)
    tree_two_indexed = index_tree(tree_two)
    # pre calculate deletion cost for each subtree
    del_source = {}
    del_target = {}
    init_deletions()
    # setup map for sub-children function
    children_scores = {}
    children_trees = {}
    # init main pass datasets
    score_matrix = [([-alignment_object.minmax_func(sys.maxsize, -sys.maxsize)] * (len(tree_two_indexed)))
                    for i in range(len(tree_one_indexed))]
    transition_matrix = [([None] * (len(tree_two_indexed))) for i in range(len(tree_one_indexed))]
    #
    for index_one in range(len(tree_one_indexed)):
        for index_two in range(len(tree_two_indexed)):
            score_list = []
            transition_list = []
            best_match_score, child_list = all_combination(tree_one_indexed[index_one], tree_two_indexed[index_two])
            # match, get best children combination from both
            cmp_value, cmp_align = alignment_object.cmp_func(tree_one_indexed[index_one].value,
                                                             tree_two_indexed[index_two].value)
            if cmp_value is not None:
                transition_value = cmp_value
                score_list.append(transition_value + best_match_score)
                transition_list.append(Tree(alignment_object.merge_func(tree_one_indexed[index_one].value,
                                                                        tree_two_indexed[index_two].value),
                                            child_list, mode='M', alignment=cmp_align))
            # insert from tree one (ignore node), get best child to index two compare (del rest of children)
            transition_value, transition_align = alignment_object.delete_func(tree_one_indexed[index_one].value,
                                                                              is_target=False)
            opt_t1_child_t2, opt_t1_children = best_of_children(tree_two_indexed[index_two],
                                                                tree_one_indexed[index_one].children, True)
            score_list.append(transition_value + opt_t1_child_t2)
            transition_list.append(Tree(tree_one_indexed[index_one].value, opt_t1_children,
                                        mode='S', alignment=transition_align))
            # insert from tree two (ignore node), get best child to index one compare (del rest of children)
            transition_value, transition_align = alignment_object.delete_func(tree_two_indexed[index_two].value,
                                                                              is_target=True)
            opt_t1_t2_child, opt_t2_children = best_of_children(tree_one_indexed[index_one],
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
    test_tree_one = Tree(1.1, [Tree(2.1, [Tree(3.1, [Tree(4.1, [])]), Tree(5.1, [Tree(6.1, [])])]),
                               Tree(2.1, [Tree(3.1, [])])])
    test_tree_two = Tree(1.2, [Tree(2.2, [Tree(4.2, [])]), Tree(2.2, [Tree(4.2, [Tree(3.2, [])])])])
    print("Tree one: {}".format(test_tree_one))
    print("Tree two: {}".format(test_tree_two))
    align_tree, align_score = align_trees(test_tree_one, test_tree_two, ar)
    print("Matched only Joined Tree ({}): {}".format(align_score, align_tree))
    print("Full Joined Tree ({}): {}".format(align_score, align_tree.full_str()))
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
