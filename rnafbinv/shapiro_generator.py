'''
File has the function to build a shapiro representation for a dot bracket structure
Also builds an index list for each motif
'''


import logging


# Object including structure, shapiro and list of sequence indexes for each shapiro motif
class ShapiroObject(object):
    def __init__(self, structure, aux, shapiro, shapiro_indexes):
        self.structure = structure
        self.aux = aux
        self.shapiro = shapiro
        self.shapiro_indexes = shapiro_indexes

    def __str__(self):
        return self.structure + '\n' + self.aux + '\n' + self.shapiro + '\n' + self.shapiro_indexes


# Increase counter by one in index map
def _map_plus_one(orig_map, index):
    value = orig_map.get(index)
    if value is None:
        orig_map[index] = 1
    else:
        orig_map[index] = value + 1


# Add str_index to the list in arr_index on orig_map_list
def _add_to_indexes(orig_map_list, arr_index, str_index):
    value = orig_map_list.get(arr_index, [])
    value.append(str_index)
    orig_map_list[arr_index] = value


def _get_indexes(orig_map_list, arr_index):
    res = "[]"
    value = orig_map_list.get(arr_index)
    if value is not None:
        res = str(value)
    return res


# generate map of closing index to bracket
def _get_closure_map(structure):
    closure_map = {}
    tracking_stack = []
    for i in range(0, len(structure)):
        if structure[i] == '(':
            tracking_stack.append(i)
        elif structure[i] == ')':
            closure_map[i] = tracking_stack.pop()
    return closure_map


# Create a string representation of dot bracket where outer parenthesis are marked with square bracket as a list
def _get_aux_list(structure):
    aux_array = list(structure)
    match_paren = [0] * len(structure)  # list of matching parenthesis
    index = 0  # current index in structure
    o = 0  # depth of open parenthesis
    while index < len(aux_array):
        if aux_array[index] == '(':
            o += 1
            match_paren[o] = index
        elif aux_array[index] == ')':
            p = index  # tracking inner closers
            while p + 1 < len(aux_array) and aux_array[p + 1] == ')' and match_paren[o - 1] == match_paren[o] - 1:
                p += 1
                o -= 1
            aux_array[p] = ']'
            index = p
            aux_array[match_paren[o]] = '['
            o -= 1
        elif aux_array[index] != '.':
            logging.error("Unknown element in dot bracket element: '{}'".format(aux_array[index]))
        index += 1
    return aux_array


# Create a string representation of dot bracket where outer parenthesis are marked with square bracket
def get_aux(structure):
    return "".join(_get_aux_list(structure))


# generate shapiro representation for given structure
def get_shapiro(structure: str) -> ShapiroObject:
    loop_size = {}  # loop sizes of a structure
    helix_size = {}  # helix sizes of a structure
    loop_degree = {0: 0}  # loop degree of a structure
    loops = 0  # number of motif start
    unpaired = 0  # track number on unpaired in motif break
    pairs = 0  # track number of paired in motif
    loop_indexes = {}  # string indexes for a loop
    stem_indexes = {}  # string indexes of stems
    bulge = [0] * int(len(structure) / 3 + 1)  # Track is an interior section is bulge (1) or internal loop (0)
    loop = [0] * int(len(structure) / 3 + 1)  # list of loop depths
    lp = 0  # Track depth of motif
    p = 0  # Track helix size

    aux_list = _get_aux_list(structure)
    closure_map = _get_closure_map(structure)
    temp = '('  # result collector with first bracket for Root
    temp_indexes = '('  # collects indexes for nucleic acids inside the motif

    for i in range(0, len(aux_list)):
        if aux_list[i] == '.':
            unpaired += 1
            _map_plus_one(loop_size, loop[lp])
            _add_to_indexes(loop_indexes, loop[lp], i)
        elif aux_list[i] == '[':
            temp_indexes += '(('
            temp += '(('
            if i > 0 and (aux_list[i - 1] == '(' or aux_list[i - 1] == '['):
                bulge[lp] = 1
            lp += 1
            loops += 1
            loop_degree[loops] = 1
            loop[lp] = loops
            bulge[lp] = 0
        elif aux_list[i] == ')':
            if aux_list[i - 1] == ']':
                bulge[lp] = 1
            p += 1
            _add_to_indexes(stem_indexes, loop[lp], i)
            _add_to_indexes(stem_indexes, loop[lp], closure_map[i])
        elif aux_list[i] == ']':
            if aux_list[i - 1] == ']':
                bulge[lp] = 1
            value = loop_degree.get(loop[lp])
            if value == 1:
                temp += 'H'  # hairpin
            elif value == 2:
                if bulge[lp] == 1:
                    temp += 'B'  # bulge
                else:
                    temp += 'I'  # internal loop
            else:
                temp += 'M'  # multi loop
            helix_size[loop[lp]] = p + 1
            _add_to_indexes(stem_indexes, loop[lp], i)
            _add_to_indexes(stem_indexes, loop[lp], closure_map[i])
            if loop_size.get(loop[lp]) is not None:
                temp += str(loop_size.get(loop[lp]))
            temp += ')'
            temp_indexes += _get_indexes(loop_indexes, loop[lp]) + ')'
            temp += 'S' + str(helix_size[loop[lp]]) + ')'
            temp_indexes += _get_indexes(stem_indexes, loop[lp]) + ')'
            pairs += p + 1
            p = 0
            lp -= 1
            _map_plus_one(loop_degree, loop[lp])
    shapiro = ''
    initial_loop_size = loop_size.get(0)
    if initial_loop_size is not None and initial_loop_size != 0:
        temp += 'E' + str(initial_loop_size) + ')'
        temp_indexes += _get_indexes(loop_indexes, 0) + ')'
        shapiro += '('
        temp_indexes = '(' + temp_indexes
    shapiro += temp + 'R)'
    shapiro_indexes = temp_indexes + '[])'
    return ShapiroObject(structure, "".join(aux_list), shapiro, shapiro_indexes)


# for testing purposes
if __name__ == "__main__":
    '''
    print("1){}\n".format(get_shapiro(
        "(((...)))")))
    print("2){}\n".format(get_shapiro(
        "(((...(((...)))...)))")))
    print("3){}\n".format(get_shapiro(
        "(((...(((...)))...(((...)))...(((...(((...)))...(((...)))...)))...)))")))
    print("4){}\n".format(get_shapiro("((((((((...(.(((((.......)))))(((........((((....)))..).((((....))))..)..))....((((....))..))...)))))))))")))
    '''
    print("5){}\n".format(get_shapiro("(((((((((.(.((.((((.......)))).)).).)))(((((((((.....)))))))))))))))")))

