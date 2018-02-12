'''
Read a tree file and generate distance matrix
'''

def find_path(graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return path
        if start not in graph.keys():
            return None
        for node in graph[start]:
            if node not in path:
                newpath = find_path(graph, node, end, path)
                if newpath: return newpath
        return None

def read_path_to_dictionary(fname):
    item = 0
    dictionary ={}
    weighed_dictionary = {}
    with open(fname) as f:
        content = f.readlines()
    for line in content:
        if '->' in line:
            begin = int(line.split('->')[0])
            end = int(line.split('->')[1].split(':')[0])
            distance = int(line.split('->')[1].split(':')[1])
            if (begin not in dictionary.keys()):
                dictionary[begin] = [end]
                weighed_dictionary[begin] = {}
                weighed_dictionary[begin][end] = distance
            else:
                dictionary[begin].append(end)
                weighed_dictionary[begin][end] = distance
        else:
            item = int(line)
    return item,dictionary, weighed_dictionary

def get_weight_of_path(weighed_graph, path):
    i = 0
    sum = 0
    while i < len(path) - 1:
        next_nodes = weighed_graph[path[i]]
        next_node = [j for j in next_nodes if j == path[i+1]]
        weight = next_nodes[next_node[0]]
        sum = sum + weight
        i = i + 1
    return  sum


if __name__ == '__main__':
    count, paths, weighed_path = read_path_to_dictionary('E:\\Share\\Code\\Fan.Biology\\tests\\dataset_10328_12.txt')
    print(count, paths, weighed_path)
    dict_matrix = {}
    matrix =[]
    for i in range(count):
        dict_matrix[i] = []
        matrix_row=[]
        for j in range(count):
            found_path = find_path(paths,i,j)
            path_weight = get_weight_of_path(weighed_path, found_path)
            dict_matrix[i].append({j:path_weight})
            matrix_row.append(path_weight)
        matrix.append(matrix_row)
    print(dict_matrix)
    print(matrix)




