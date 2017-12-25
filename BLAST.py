import ast
import sys
from suffix_tree import SuffixTree
from fm_index import FMIndex
from data_preprocess import *

def get_match_position_with_suffixtree(neighbors, suffix_trees):
    all_match_position = []
    match_position = []
    for i in range(len(neighbors)):
        for j in range(len(suffix_trees)):
            temp = {}
            for n in neighbors[i]:
                if suffix_trees[j].has_substring(n):
                    p = suffix_trees[j].find_substring(n)
                    if n not in temp:
                        temp[n] = [p]
                    else:
                        temp[n].append(p)
            match_position.append(temp)
        all_match_position.append(match_position)
    return all_match_position


def get_match_position_with_fmindex(neighbors, fm_indexs):
    all_match_position = []
    match_position = []
    for i in range(len(neighbors)):
        for j in range(len(fm_indexs)):
            temp = {}
            for n in neighbors[i]:
                p = fm_indexs[j].get_offset1(n)
                if len(p) > 0:
                    if n not in temp:
                        temp[n] = p
                    else:
                        temp[n].extend(p)
            match_position.append(temp)
        all_match_position.append(match_position)
    return all_match_position


def get_match_position_with_fmindex_file(neighbors, file_name, a, b):

    protein_name =[]
    fread = open(file_name, "r")
    fm_index = FMIndex("a",1,1)
    all_match_position = []
    match_position = []
    fm_indexs = []
    name = fread.readline()
    while name != "":
        protein_name.append(name)
        bwt = fread.readline()[:-1]
        fcol = ast.literal_eval(fread.readline()[:-1])
        tally = ast.literal_eval(fread.readline()[:-1])
        sa = ast.literal_eval(fread.readline()[:-1])
        fm_indexs.append([bwt, fcol, tally, sa])
        name = fread.readline()[:-1]

    for i in range(len(neighbors)):
        for j in range(len(fm_indexs)):
            temp = {}
            for n in neighbors[i]:
                p = fm_index.get_offset(n, a, b, fm_indexs[j][0], fm_indexs[j][1], fm_indexs[j][2], fm_indexs[j][3])
                if len(p) > 0:
                    if n not in temp:
                        temp[n] = p
                    else:
                        temp.extend(p)
            match_position.append(temp)
        all_match_position.append(match_position)
    return all_match_position, protein_name

def get_match_position_org(protein_sequences,query_neighbors, k):

    all_positions = []
    positions = []
    for i in range(len(query_neighbors)):
        for sequece in protein_sequences:
            temp = {}
            for j in range(len(sequece) - k + 1):
                sub_str = sequece[j:j+k]
                if sub_str in query_neighbors[i]:
                    if sub_str not in temp:
                        temp[sub_str] = [j]
                    else:
                        temp[sub_str].append(j)
            positions.append(temp)
        all_positions.append(positions)
    return all_positions

def main():

    a = 3
    b = 3
    read = open("query_output.txt", "r")
    hit_inputs = ["-l", "-s", "-f"]
    extend_inputs =["-d", "-s", "-g"]
    print "Please type your command:",
    inputs = sys.stdin.readline()
    for i in range(len(inputs)):
        if inputs[i] != " ":
            hit = inputs[i:i+2]
            break
    extend = inputs[-3:-1]
    while hit not in hit_inputs or extend not in extend_inputs:
        print "Wrong command, please type again: ",
        input = sys.stdin.readline()
        for i in range(len(input)):
            if input[i] != " ":
                hit = input[i:2]
                break
        extend = input[-3:-1]
    query_sequences = []
    query_sequence = ""
    while True:
        s = read.readline()
        if s == "":
            if query_sequence != "":
                query_sequences.append(query_sequence)
                break
        if s[0] == ">":
            if query_sequence != "":
                query_sequences.append(query_sequence)
                query_sequence = ""
        else:
            if s[-1] == "\n":
                s = s[:-1]
            query_sequence += s
    k_size = 3
    threshold = 11
    query_neighbors = generate_all_neighbors(query_sequences, k_size, threshold)


    protein_name = []
    protein_sequences = get_protein_data("ecoliProtSeq.fasta", protein_name)

    if hit == "-l":
        all_match_positions = get_match_position_org(protein_sequences, query_neighbors, k_size)


    if hit=="-s":
        proteins_suffix_trees = []
        for i in range(len(protein_sequences)):
             proteins_suffix_trees.append(SuffixTree(protein_sequences[i]))

        all_match_positions = get_match_position_with_suffixtree(query_neighbors, proteins_suffix_trees)


    if hit == "-f":
        protein_fmindexs =[]
        for i in range(len(protein_sequences)):
             protein_fmindexs.append(FMIndex(protein_sequences[i]+"$", a, b))

        all_match_positions = get_match_position_with_fmindex(query_neighbors, protein_fmindexs)

    # print len(all_match_positions)

    result = []
    for k in range(len(all_match_positions)):
        for i in range(len(protein_sequences)):
            result.append({})
            for l in range(len(query_neighbors)):
                for k_mer in query_neighbors[l]:
                    for j in range(len(query_neighbors[l][k_mer])):
                        if k_mer in all_match_positions[k][i]:
                            for pos in all_match_positions[k][i][k_mer]:
                                try:
                                    if extend == "-g":
                                        score, q_offset, s_offset, comp_seg, q_seg, s_seg = extending_with_gap(query_sequences[l], query_neighbors[l][k_mer][j], protein_sequences[i], pos, k_size, -5, 10)
                                    if extend == "-s":
                                        score, q_offset, s_offset, comp_seg, q_seg, s_seg = extending_with_score(
                                            query_sequences[l], query_neighbors[l][k_mer][j], protein_sequences[i], pos,
                                            k_size, 10)
                                    if extend == "-d":
                                        score, q_offset, s_offset, comp_seg, q_seg, s_seg = extending_with_drop(
                                            query_sequences[l], query_neighbors[l][k_mer][j], protein_sequences[i], pos,
                                            k_size, 3)
                                except:
                                    continue
                                list = [q_offset, s_offset, comp_seg, q_seg, s_seg]
                                if (score in result):
                                    result[i][score].append(list)
                                else:
                                    result[i][score] = list
        #     print "i", i
        # print "k", k
    sys.stdout = open("output_result.txt", "w")
    for i in range(len(protein_sequences)):
        print(protein_name[i])
        display = result[i]
        sort_order = sorted(display, reverse= True)
        m = 0
        for j in sort_order:
            k = 0
            m += 1
            while k < len(display[j]):
                score = j
                q_offset = display[j][k]
                k+=1
                s_offset = display[j][k]
                k+=1
                comp_seg = display[j][k]
                k+=1
                q_seg = display[j][k]
                k+=1
                s_seg = display[j][k]
                k+=1
                print "Score: ", score
                print "Query:  ", "{0:4d}".format(q_offset), "  ", q_seg, "  ", q_offset+len(q_seg)
                print "                ", comp_seg
                print "Sbjct:  ", "{0:4d}".format(s_offset), "  ", s_seg, "  ", s_offset + len(s_seg)
                print ""
        print("-----------------------------------------------------------------------------------------")

if __name__ == "__main__":
    main()