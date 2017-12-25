import singlescore

smatrix=singlescore.createScoreM(singlescore.scores)


def get_single_score(x, y):

    ii = singlescore.alp.index(x)
    jj = singlescore.alp.index(y)
    return int(smatrix[ii][jj])

#extact relative info from new query and new sequence
def extract_result(pairscore, new_query, new_sequence):
    Score = pairscore
    comparison_seg = ''
    Identities = 0
    Positives = 0
    Gaps = 0
    for (q_c, s_c) in zip(new_query, new_sequence):
        if q_c == '-' or s_c == '-' :
            comparison_seg += ' '
            Gaps +=1
        elif q_c == s_c:
            comparison_seg += q_c
            Identities += 1
            singlescore = get_single_score(q_c, s_c)
            if singlescore > 0:
                Positives += 1
        else:
            singlescore = get_single_score(q_c, s_c)
            if singlescore > 0:
                comparison_seg += '+'
                Positives += 1
            else:
                comparison_seg += ' '
    results = {}
    results['Score'] = Score
    results['Identities'] = 1.0 * Identities/len(new_query)
    results['Positives'] = 1.0 * Positives/len(new_query)
    results['Gaps'] = 1.0 * Gaps/len(new_query)
    return results, comparison_seg

#extending method one: drop without gap
def extending_with_drop(query, q, sequence, s, k, drop):
    query_seg = query[q:q+k]
    sequence_seg = sequence[s:s+k]
    # computing the score for the word
    pairscore = 0
    for i in range(0,k):
        score = get_single_score(query[q+i], sequence[s+i])
        pairscore += score
    #extending to the left
    left = -1
    left_extend = 0
    while q+left > -1 and s+left > -1:
        mid_pair_score=get_single_score(query[q+left], sequence[s+left])
        if mid_pair_score < -drop+1:
            left_extend = left
            break
        else:
            pairscore+=mid_pair_score
            left -= 1
    #extending to the right
    right = 1
    right_extend = 0
    while q+k+right<len(query)       and s+k+right<len(sequence):
        mid_pair_score = get_single_score(query[q+k+right], sequence[s+k+right])
        if mid_pair_score < -drop+1:
            right_extend = right
            break
        else:
            pairscore+=mid_pair_score
            right+=1

    new_query = query[q+left_extend:q] + query_seg + query[q+k:q+k+right_extend]
    new_sequence = sequence[s+left_extend:s] + sequence_seg + sequence[s+k:s+k+right_extend]
    query_offset = q + left_extend
    sequence_offset = s + left_extend
    result = extract_result(pairscore, new_query, new_sequence)
    return result[0]['Score'], query_offset, sequence_offset, result[1], new_query, new_sequence

#extending method 2: score distance without gap
def extending_with_score(query, q, sequence, s, k, score_distance):
    query_seg = query[q:q+k]
    sequence_seg=sequence[s:s+k]
    pairscore = 0
    for i in range(0,k):
        score = get_single_score(query[q+i], sequence[s+i])
        pairscore += score
    #extending to the left
    left = -1
    left_max = (0, float('-inf'))
    left_score = 0
    while q+left > -1 and s+left > -1:
        mid_pair_score=get_single_score(query[q+left], sequence[s+left])
        left_score +=mid_pair_score
        if left_score > left_max[1]:
            left_max = (left, mid_pair_score)
        if abs(left_max[1] - left_score) >= score_distance:
            break
        left -= 1
    #extending to the right
    right = 1
    right_max = (0, float('-inf'))
    right_score = 0
    while q+right < len(query) and s+right < len(sequence):
        mid_pair_score=get_single_score(query[q+right], sequence[s+right])
        right_score +=mid_pair_score
        if right_score > right_max[1]:
            right_max = (right, mid_pair_score)
        if abs(right_max[1] - right_score) >= score_distance:
            break
        right += 1

    left_query = query[q + left_max[0]:q]
    right_query = query[q + k:q + k + right_max[0]]
    left_sequence = sequence[q + left_max[0]:q]
    right_sequence = sequence[q + k:q + k + right_max[0]]
    new_query = left_query + query_seg + right_query
    new_sequence = left_sequence + sequence_seg + right_sequence
    query_offset = q + left_max[0]
    sequence_offset = s + left_max[0]

    Score = left_max[1] + pairscore + right_max[1]
    result = extract_result(Score, new_query, new_sequence)
    return result[0]["Score"], query_offset, sequence_offset, result[1], new_query, new_sequence

#extending method 3: score distance with gap
def extending_with_gap(query, q, sequence, s, k, gap_score, score_distance):
    query_seg = query[q:q+k]
    sequence_seg=sequence[s:s+k]
    #compute the score for the word
    pairscore = 0
    for i in range(0,k):
        score = get_single_score(query[q+i], sequence[s+i])
        pairscore += score

    #extending to the left
    left = -1
    left_score = []
    left_row_max = []
    left_max = 0
    directions = []
    while q+left > -1 and s+left > -1:
        mid_pair_score=get_single_score(query[q+left], sequence[s+left])
        #fill the top-left
        if left == -1:
            three_directions = [gap_score+gap_score, gap_score+gap_score, mid_pair_score]
            mid_fill = max(three_directions)
            left_score.append([mid_fill])
            directions.append([three_directions.index(mid_fill)])
            left_max = mid_fill
        else:
            #col first
            three_directions = [left_score[0][abs(left)-2]+gap_score,gap_score * abs(left)+gap_score,gap_score * (abs(left)-1)+get_single_score(query[q - 1], sequence[s+left])]
            mid_fill = max(three_directions)
            left_score[0].append(mid_fill)
            directions[0].append(three_directions.index(mid_fill))

            #row first
            three_directions = [gap_score * abs(left)+gap_score,left_score[abs(left)-2][0]+gap_score,gap_score * (abs(left)-1)+get_single_score(query[q + left], sequence[s - 1])]
            mid_fill = max(three_directions)
            left_score.append([mid_fill])
            directions.append([three_directions.index(mid_fill)])

            #fill row
            for i in range(1,abs(left)-1):
                com_fill = get_single_score(query[q+left], sequence[s-i-1])
                three_directions = [left_score[abs(left)-1][i-1]+gap_score,left_score[abs(left)-2][i]+gap_score,left_score[abs(left)-2][i-1]+com_fill]
                mid_fill = max(three_directions)
                left_score[abs(left)-1].append(mid_fill)
                directions[abs(left)-1].append(three_directions.index(mid_fill))

            #fill col
            for i in range(1,abs(left)-1):
                com_fill = get_single_score(query[q-i-1], sequence[s+left])
                three_directions = [left_score[i][abs(left)-2]+gap_score, left_score[i-1][abs(left)-1]+gap_score,left_score[i-1][abs(left)-2]+com_fill]
                mid_fill = max(three_directions)
                left_score[i].append(mid_fill)
                directions[i].append(three_directions.index(mid_fill))

            #fill bottom-right
            three_directions = [left_score[abs(left)-1][abs(left)-2]+gap_score, left_score[abs(left)-2][abs(left)-1]+gap_score,left_score[abs(left)-2][abs(left)-2]+mid_pair_score]
            mid_fill = max(three_directions)
            left_score[abs(left)-1].append(mid_fill)
            directions[abs(left)-1].append(three_directions.index(mid_fill))

        #record row max
        mid_fill_row_max = max(left_score[abs(left)-1])
        if left_max < mid_fill_row_max:
            left_max = mid_fill_row_max
        #when to stop extending
        if len(left_row_max)>0 and abs(left_max-left_row_max[-1]) >= score_distance:
            break
        left_row_max.append(mid_fill_row_max)
        left -= 1


    #left trace back
    new_query = ''
    new_sequence = ''
    left_extend_len = left_row_max.index(left_max) + 1

    query_index = left_extend_len
    left_total_score = left_max
    sequence_index = left_score[query_index - 1][:query_index].index(left_total_score)+1

    while query_index - 1 > -1 and sequence_index - 1 > -1:
        if directions[query_index - 1][sequence_index - 1] == 0:
            new_sequence += sequence[s-sequence_index]
            new_query += '-'
            sequence_index -= 1
        elif directions[query_index - 1][sequence_index - 1] == 1:
            new_sequence += '-'
            new_query += query[q-query_index]
            query_index -= 1
        elif directions[query_index - 1][sequence_index - 1] == 2:
            new_sequence += sequence[s-sequence_index]
            new_query += query[q-query_index]
            sequence_index -= 1
            query_index -= 1
    if query_index > 0:
        new_query += query[q-query_index:q]
        new_sequence += ''.join(['-' for i in range(query_index)])
    if sequence_index > 0:
        new_query += ''.join(['-' for i in range(sequence_index)])
        new_sequence += sequence[s-sequence_index:s]
    left_new_query = new_query
    left_new_sequence = new_sequence


    #fill right matrix
    right = 1
    right_score = []
    right_row_max = []
    right_max = 0
    directions = []
    q2 = q + k - 1
    s2 = s + k - 1
    while q2+right < len(query) and s2+right < len(sequence):
        mid_pair_score=get_single_score(query[q2+right], sequence[s2+right])
        #fill the top-left
        if right == 1:
            three_directions = [gap_score+gap_score, gap_score+gap_score, mid_pair_score]
            mid_fill = max(three_directions)
            right_score.append([mid_fill])
            directions.append([three_directions.index(mid_fill)])
            right_max = mid_fill
        else:
            #col first
            three_directions = [right_score[0][abs(right)-2]+gap_score,gap_score * abs(right)+gap_score,gap_score * (abs(right)-1)+get_single_score(query[q2 + 1], sequence[s2+right])]
            mid_fill = max(three_directions)
            right_score[0].append(mid_fill)
            directions[0].append(three_directions.index(mid_fill))

            #row first
            three_directions = [gap_score * abs(right)+gap_score,right_score[abs(right)-2][0]+gap_score,gap_score * (abs(right)-1)+get_single_score(query[q2 + right], sequence[s2 + 1])]
            mid_fill = max(three_directions)
            right_score.append([mid_fill])
            directions.append([three_directions.index(mid_fill)])

            #fill row
            for i in range(1,abs(right)-1):
                com_fill = get_single_score(query[q2+right], sequence[s2+i+1])
                three_directions = [right_score[abs(right)-1][i-1]+gap_score,right_score[abs(right)-2][i]+gap_score,right_score[abs(right)-2][i-1]+com_fill]
                mid_fill = max(three_directions)
                right_score[abs(right)-1].append(mid_fill)
                directions[abs(right)-1].append(three_directions.index(mid_fill))

            #fill col
            for i in range(1,abs(right)-1):
                com_fill = get_single_score(query[q2+i+1], sequence[s2+right])
                three_directions = [right_score[i][abs(right)-2]+gap_score, right_score[i-1][abs(right)-1]+gap_score, right_score[i-1][abs(right)-2]+com_fill]
                mid_fill = max(three_directions)
                right_score[i].append(mid_fill)
                directions[i].append(three_directions.index(mid_fill))

            #fill bottom-right
            three_directions = [right_score[abs(right)-1][abs(right)-2]+gap_score, right_score[abs(right)-2][abs(right)-1]+gap_score,right_score[abs(right)-2][abs(right)-2]+mid_pair_score]
            mid_fill = max(three_directions)
            right_score[abs(right)-1].append(mid_fill)
            directions[abs(right)-1].append(three_directions.index(mid_fill))

        #record row max
        mid_fill_row_max = max(right_score[abs(right)-1])
        if right_max < mid_fill_row_max:
            right_max = mid_fill_row_max

        #when to stop extending
        if len(right_row_max)>0 and abs(right_max-right_row_max[-1]) >= score_distance:
            break
        right_row_max.append(mid_fill_row_max)
        right += 1

    #right trace back
    new_query = ''
    new_sequence = ''
    right_extend_len = len(right_row_max)

    query_index = right_row_max.index(right_max) + 1
    right_total_score = right_max
    sequence_index = right_score[query_index - 1][:query_index].index(right_total_score)+1
    while query_index - 1 > -1 and sequence_index - 1 > -1:
        if directions[query_index - 1][sequence_index - 1] == 0:
            new_sequence = sequence[s2+sequence_index] + new_sequence
            new_query = '-' + new_query
            sequence_index -= 1
        elif directions[query_index - 1][sequence_index - 1] == 1:
            new_sequence = '-' + new_sequence
            new_query = query[q2+query_index] + new_query
            query_index -= 1
        elif directions[query_index - 1][sequence_index - 1] == 2:
            new_sequence = sequence[s2+sequence_index] + new_sequence
            new_query = query[q2+query_index] + new_query
            sequence_index -= 1
            query_index -= 1
    if query_index > 0:
        new_query = query[q2 + 1:q2 + query_index + 1] + new_query
        new_sequence = ''.join(['-' for i in range(query_index)]) + new_sequence
    if sequence_index > 0:
        new_query = ''.join(['-' for i in range(sequence_index)]) + new_query
        new_sequence = sequence[s2 + 1:s2 + sequence_index+1] + new_sequence
    right_new_query = new_query
    right_new_sequence = new_sequence


    new_query = left_new_query + query_seg + right_new_query
    new_sequence = left_new_sequence + sequence_seg + right_new_sequence
    query_offset = q-left_extend_len
    sequence_offset = s-left_extend_len

    Score = left_total_score + pairscore + right_total_score
    result = extract_result(Score, new_query, new_sequence)
    return result[0]['Score'], query_offset, sequence_offset, result[1], new_query, new_sequence

def extending_test():

    query = 'MRFNRHYNYCTPMQPATRPFTFH'
    query += 'DNFVFGQ'
    q = 18
    sequence = 'QTISGEHGLDGSGVYNGSSDLQLERMNVYFNEASNNKYVPRAVLVDLEPGTMDAVRAGPFGQLFRPDNFVFGQSGAGNNWAKGHYTEG'
    s = 58
    k = 3
    drop = 3
    distance_score = 10
    gap_score = -5
    # print '-------------------------------------'
    #re = extending(query, q, sequence, s, k, drop)
    #re = extending2(query, q, sequence, s, k, distance_score)
    re = extending3(query, q, sequence, s, k, gap_score, distance_score)
    print re[4]
    print re[3]
    print re[5]
    print re
    # score, q_offset, s_offset, comp_seg, q_seg, s_seg = extending2(query, q, sequence, s, k, drop, 5)
    # print score
    # print q_offset, s_offset
    # print q_seg
    # print comp_seg
    # print s_seg

if __name__ == "__main__":
    #main()
    extending_test()
