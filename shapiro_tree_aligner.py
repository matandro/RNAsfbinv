import re

import shapiro_generator
import tree_aligner
import logging
import IUPAC


def align_sequences(value_one, value_two):
    score_matrix = IUPAC.align_iupac_dna_sequence(value_one.sequence, value_two.sequence, mismatch_penalty=-5)
    return IUPAC.get_best_score(score_matrix)


def cmp_shapiro_tree_values(value_one, value_two):
    result = 0
    value_one_name = value_one.name
    value_two_name = value_two.name
    if value_one_name[0] == value_two_name[0] or \
            ((value_one_name[0] == 'M' or value_one_name[0] == 'I' or value_one_name[0] == 'B' or value_one_name[
                0] == 'H') and
             (value_two_name[0] == 'M' or value_two_name[0] == 'I' or value_two_name[0] == 'B' or value_two_name[
                 0] == 'H')):
        #result += 1 + (align_sequences(value_one, value_two) / max(len(value_one.sequence), len(value_two.sequence)) if value_one_name != 'R' else 0)
        result += 1 + (min(align_sequences(value_one, value_two), 1) if value_one_name != 'R' else 0)
    return result


def merge_shapiro_tree_values(value_one, value_two):
    return value_one


#actually merge string and length etc...
def smart_merge_shapiro_tree_values(value_one, value_two):
    raise NotImplementedError


def insert_shapiro_func(value):
    return -0.125


def delete_shapiro_func(value):
    sequence_false = (sum([(5 if char != 'N' else 0) for char in value.sequence]) / len(value.sequence) if len(value.sequence) > 0 else 0)
    return -5.0 - sequence_false


# recursively turn a shapiro string into a tree
def shapiro_to_tree(shapiro_str, shapiro_index, sequence):
    bracket_stack = []
    str_index = len(shapiro_str) - 1
    bp_index = len(shapiro_index) - 1
    # empty string or not proper string (must start with bracket
    if str_index < 0 or shapiro_str[str_index] != ')':
        if str_index >= 0:
            logging.error("Shapiro does not end with close bracket: {}".format(shapiro_str))
        return None
    # get value of top tree
    bracket_count = 1
    str_index -= 1
    bp_index -= 1
    current_value = ""
    current_bp_list = ""
    while shapiro_str[str_index] != ')' and shapiro_str[str_index] != '(':
        current_value = shapiro_str[str_index] + current_value
        str_index -= 1
    while shapiro_index[bp_index] != ')' and shapiro_index[bp_index] != '(':
        current_bp_list = shapiro_index[bp_index] + current_bp_list
        bp_index -= 1
    res_tree = tree_aligner.Tree(ShapiroTreeValue(current_value, current_bp_list, sequence), [])
    if shapiro_str[str_index] == '(':
        bracket_count -= 1
    # split sub groups and recursively send them forward
    str_start_index = str_index
    bp_start_index = bp_index
    while bracket_count > 0:
        if shapiro_str[str_index] == ')':
            bracket_count += 1
            while shapiro_index[bp_index] != ')':
                bp_index -= 1
        elif shapiro_str[str_index] == '(':
            bracket_count -= 1
            while shapiro_index[bp_index] != '(':
                bp_index -= 1
        if bracket_count == 1:
            res_tree.add_child(shapiro_to_tree(shapiro_str[str_index:str_start_index + 1],
                                               shapiro_index[bp_index:bp_start_index + 1], sequence))
            str_start_index = str_index - 1
            bp_start_index = bp_index - 1
        str_index -= 1
        bp_index -= 1
    return res_tree


class ShapiroTreeValue:
    def __init__(self, shapiro_str_name, shapiro_index, sequence):
        self.name = shapiro_str_name[0]
        try:
            self.size = (int(shapiro_str_name[1:]) if shapiro_str_name[1:] != '' else 0)
        except ValueError:
            self.size = 0
            logging.error("ShapiroTreeValue.__init__ error parsing ShapiroTreeValue({}, {}, {})".format(shapiro_str_name, shapiro_index, sequence))
        if self.size == 0:
            self.index_list = []
            self.sequence = ''
        else:
            self.index_list = []
            for index in shapiro_index[1:len(shapiro_index) - 1].split(','):
                index = index.strip()
                if index.isdigit():
                    self.index_list.append(int(index))
                else:
                    logging.error("Error decyphering shapiro: name={}, index={}, sequence={}".format(shapiro_str_name, shapiro_index, sequence))
            self.index_list.sort()
            indexes = []
            last = None
            current = []
            for index in self.index_list:
                if last is None or index == last + 1:
                    current.append(index)
                else:
                    indexes.append(current)
                    current = [index]
                last = index
            indexes.append(current)
            self.sequence = ".".join(["".join([sequence[index] for index in group]) for group in indexes])

    def __str__(self):
        return "{}{}({})".format(self.name, self.size, self.sequence)

    '''
    def __init__(self, name, size, index_list):
        self.name = name
        self.size = size
        self.index_list = index_list
    '''


def get_tree(structure, sequence):
    shapiro = shapiro_generator.get_shapiro(structure)
    return shapiro_to_tree(shapiro.shapiro, shapiro.shapiro_indexes, sequence)


def align_trees(tree_source, tree_target,
                alignment_rules=tree_aligner.AlignmentRules(insert_func=insert_shapiro_func,
                                                            delete_func=delete_shapiro_func,
                                                            cmp_func=cmp_shapiro_tree_values,
                                                            merge_func=merge_shapiro_tree_values)):
    return tree_aligner.align_trees(tree_source, tree_target, alignment_rules)


def align_shapiro(shapiro_source, sequence_source, shapiro_target, sequence_target,
                  alignment_rules=tree_aligner.AlignmentRules(insert_func=insert_shapiro_func,
                                                              delete_func=delete_shapiro_func,
                                                              cmp_func=cmp_shapiro_tree_values,
                                                              merge_func=merge_shapiro_tree_values)):
    tree_source = shapiro_to_tree(sequence_source.shapiro, shapiro_source.shapiro_indexesm, sequence_source)
    tree_target = shapiro_to_tree(shapiro_target.shapiro, shapiro_target.shapiro_indexesm, sequence_target)
    return tree_aligner.align_trees(tree_source, tree_target,alignment_rules)


if __name__ == "__main__":
    # shapiro_one = shapiro_generator.get_shapiro("(((...(((...)))...(((...)))...(((...(((...)))...(((...)))...)))...)))")
    #:shapiro_two = shapiro_generator.get_shapiro("((((((((.....(((((.......)))))..........(((((((.....)))))))..))))))))")
    #"(((...(((...)))...(((...)))...(((...(((...)))...(((...)))...)))...)))")
    shapiro_one = shapiro_generator.get_shapiro("(.(.).(.).(.(..(..(...)..).)).)")
    shapiro_two = shapiro_generator.get_shapiro("(((((((((.(.((.((((.......)))).)).).)))(((((((((.....)))))))))))))))")
    tree_one = shapiro_to_tree(shapiro_one.shapiro, shapiro_one.shapiro_indexes, 'N' * len(shapiro_one.structure))
    tree_two = shapiro_to_tree(shapiro_two.shapiro, shapiro_two.shapiro_indexes, 'N' * len(shapiro_two.structure))
    print("Shapiro one: {}\nTree one: {}".format(shapiro_one.shapiro, tree_one))
    print("Shapiro two: {}\nTree two: {}".format(shapiro_two.shapiro, tree_two))
    aligned_tree, aligned_score = tree_aligner.align_trees(
        tree_one, tree_two,
        tree_aligner.AlignmentRules(insert_func=insert_shapiro_func, delete_func=delete_shapiro_func,
                                    cmp_func=cmp_shapiro_tree_values, merge_func=merge_shapiro_tree_values))
    print("Aligned tree ({}): {}".format(aligned_score, aligned_tree))
