from extending  import *
import copy

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

def get_scores(s, t):
    score = 0
    for i in range(len(s)):
        score += get_single_score(s[i], t[i])
    return score


def get_k_mers(query_seq, k):
    n = len(query_seq)
    if k > n:
        return []
    result = []
    for i in range(k, n+1):
        result.append(query_seq[i-k:i])
    return result


def generate_neighbors(word, threshold):
    k = len(word)
    collection = list(amino_acids)
    result = []
    for i in range(k-1):
        temp = copy.deepcopy(collection)
        collection = []
        for s in temp:
            for a in amino_acids:
                collection.append(s+a)
    for k_mer in collection:
        if get_scores(k_mer, word) >= threshold:
            result.append(k_mer)
    return result


def generate_all_neighbors(query_sequences, k, threshold):

    k_mers = []
    for query_sequence in query_sequences:
        k_mers.append(get_k_mers(query_sequence, k))
    all_neighbors = []
    for i in range (len (query_sequences)):
        neighbors = {}
        for j in range(len(k_mers[i])):
            neighbor = generate_neighbors(k_mers[i][j], threshold)
            for n in neighbor:
                if n not in neighbors:
                    neighbors[n] = [j,]
                else:
                    neighbors[n].append(j)
        all_neighbors.append(neighbors)
    return all_neighbors


def get_protein_data(file_name, protein_name):
    read_sequence = open(file_name, "r")
    result = []

    protein_name.append(read_sequence.readline()[:-1])
    while True:
        seq = read_sequence.readline()
        times = 0
        protein = ""
        while len(seq) > 0 and seq[0] != ">":
            times += 1
            protein += seq[:-1]
            seq = read_sequence.readline()
        protein_name.append(seq[:-1])

        if protein != "":
            result.append(protein)
        if times == 0:
            break
    return result