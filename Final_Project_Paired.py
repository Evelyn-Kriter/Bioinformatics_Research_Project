from collections import Counter

#import xlsxwriter





### TEST RUN DEFAULTS AND SAMPLES

DEFAULT_N = 50



#length of gene1 = 774

#length of gene2 = 766

gene_1 = "TTAACAATTTTTCCATAGTTTTAAAAAATGAAAGAATTTTTTATAACCTAAATTAAAAACACTAGATACCATAAAATAAGGTGTTATTATTAAATTTATTTTTTTTATTTTACTAAAAAAAAAAACATTTTCAATTCCATCTATTATAATCCAATAATTTATATTATTTAAATTTTTAATATATAAAAAATATAAAAAAAATTTTATTAATTTAAAAAATATTTTAAAAAATAATAAATTTATATAAAAAATATAAATTTTAGAAATAAAATTATTTTTATAATTAAAATTGTTTAATAAAAAATTTTCAGAAATATTAATATAACATTTATAATAAAAATATTTGAAATTAATTTTTATAAAAATAATTTTACCTTGTTTTCCAAAATTAGCAATTTTTGAATTATAAAAATTTTTATAACAAAAATTTCCATTACCTCCTTTTCCTCCTTTTAATATTTTAATAAAATAATTGTTTTTAACAATATAAATTTTTTTTTTATTTAATGTTATTATTACACCTAAAGATATTTTTATTATTAAATCTTTACCTTTTTTTCCTTTTTTGTTATTTTTTTCACCATTTTTACCATTTTTTGAAATAAATAAATTATTATTTGGTAATTTAAAATTATTATTAGATAATAAATATAAATCACCACCATCACCACCATCACCTCCATTAGCATAAAATTTTCCATTTAAAAATTTATAACTAATACAACCATTTCCTCCTTTACCACTTTTTAAAAGTAAATATTTTTCAATAAACAT"

gene_2 = "TATTAACTTCTATAAAAGTTTTATTTTTTTTTTTTTTAATTAATACTATTCCATTTATTAATGACTGTATATAATAATTTTTTGAATAAAAAGTATTTTTTCCTGAAATTATTTTAGAACCATTTTGTTTAATTATTATAGATCCTTTATTAACATAATTATTATTATAAATTTTTAATCCTAATCTTTTAGAAATAGAATCTCTTCCATTTCTCGTACTACCTCCAGCTTTTTTTTGTGCCATATTAAATTATATTTTGTAAAAATAATAAACTTTTTTTTTTATAATTTATATTAGTTTTTAAAAAATTTTTTCTTCTTTTTTTTTTTATAGATATAGATTTTATATAAAAATGTTTTAAAATTATAAAACATGCTTTTAAATTTAAACCTTTTTTATCTAAAATTATTTTTTTATTAAAAAAAAAAAATATTTTTTTAAAAAAAATTTTTTTACCAATATTTTTATTTATATAATCAACAATTAAAAAATTATTTAATTTTGAAAAATAAATTTTTTTACCTATACAAAAAGATAACAAATTATTTAATTTTATCATTTTTTAAAATATTATTATTATAATAATTAAAAACACCAAGTTCACCAATTAAAGAACAATCACCTTTTCTATATCCTATTTTTATTAAATAAAAAAAAATTTTTTTTTTAAAATAATTATTTAAAAGAATAAAAATTTTAATTTTATTTTTTTTTAATTTAAAATAATTACTTTTAATTTTACCAAATTTTATTAAATTAAAAAAA"

gene_3 = "TTTTTTAAAGCATTTGATATTTCTAAAAATACAAATTTATTCATATTTTTTAATTTAAATAAATCAAATTCTGATAATTTAATCAAGTCTCCTATTAAAAAAATGTTGTTTCTTTTTAATATATTTGAAGTTTTTATACTTAATTCTAAATTATCTACTGATCTAATAAATATAGGATTAATACTTAAAAAATTTTTTTTTTTAATTTTTTTTGTTTTTTTTTTTCCTAAAACTGAAAAAAAAACATCAAAATATTTTTTTATATAAAATATACAATTATTAAAACATTCTTCAGGTTTAATAATTCCGTCAGTTTCAATATCAAAAAATAAATTTTTAATTTTTTTATTAAAAAATTTTTTATGTATAAAATAATTTATGTTTTTTAAAGAAGACTTTAAAAAATTAATTTTAATAACTTTTGAATTAAAAATTTTATTAATAATATTATTATTAGTAAAATTATAACTAGTATTAATACATTTCATTAAAACAAAAAAAATTATATTATCTGTTATATTTGCTATTATTATATCTGGATTATAAATAATAATATGTTTATCAGAAAAAATATCTTTTGCTTTAACAATACAAGGTCCTTTTTTTTTTATAATTAAAAAAGCAGTTATTGAATTTTTAATTTTTATTAAAATATTATTTATATTATTAATAATTTTTAAAGTATTTTCTTTTATACCTTTTATATTAGAAAATTCAGAA"

gene = gene_1 + gene_2 + gene_3

partial_gene1 = "TTATCCATGTGTGT" ##For testing



#paired_k_values = [5, 10, 20, 40, 80, 160, 320, 640, 720, 774]

#paired_k_values = [2, 5, 6, 7, 8, 9, 10 , 20, 40, 80, 160, 200, 250, 300, 320, 350, 386]
#paired_k_values = [12, 14, 16, 18]
K = 8

D = 2

D_values = [0, 1, 2, 4, 8, 16, 32, 64, 100, 110, 120, 130, 140, 150, 160, 174, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200]
#D_values = [0, 1, 2, 4, 8, 16, 32, 64, 100, 110, 120, 130, 140, 150, 160, 174, 300, 400, 500, 600, 700, 800, 900, 1000, 1100]

#D_values = [0, 1, 2, 4, 8, 32, 100, 120, 140, 160, 600, 800, 1000, 1200]





### CODE START###



### Read a gene

def readGN(filename):

    with open(filename, "r") as file:

        sequence = [line.strip() for line in file.readlines()]

        return sequence

    

    

def paired_kmer_composition(sequence, k, d):

    unique_paired_kmers = []

    for i in range(k+d, len(sequence)-k+1):

        unique_paired_kmers.append(tuple((sequence[i-k-d:i-d],sequence[i:i+k])))

    return unique_paired_kmers





def construct_paired_debruijn(paired_kmers):

    pn_kmers = paired_kmers

    pn_kmers.sort()

    

    graph = {}

    for i in range(0, len(pn_kmers)):

        

        #turn every prefix into a key in the dictionary and add suffix

        prefix1 = pn_kmers[i][0][:-1]

        suffix1 = pn_kmers[i][0][1:]

        

        prefix2 = pn_kmers[i][1][:-1]

        suffix2 = pn_kmers[i][1][1:]

        

        prefix = (prefix1, prefix2)

        suffix = (suffix1, suffix2)

        

        if prefix in graph:

            graph[prefix].append(suffix)

        else:

            graph[prefix] = [suffix]

            

    return graph





def generate_contigs(graph):

    num_incoming = Counter()

    num_outgoing = Counter()

    

    #calculate in-edges

    for node, edges in graph.items():

        num_outgoing[node] += len(edges)

        contigs = []

        for edge in edges:

            num_incoming[edge] += 1

            

        for v in graph: #range(0, len(graph[node])): #length of KEYS in graph

            if num_incoming[v] != 1 or num_outgoing[v] != 1: #v is not a 1-in-1-out node

                if num_outgoing[v] > 0:

                    for w in graph[v]: #range(0, len(num_outgoing)): #for each outgoing edge (v, w) from v

                        #nonbranchingpath = the path consisting of the single edge (v, w)

                        paths = [v, w]

                        while num_incoming[w] == 1 and num_outgoing[w] == 1: #w is a 1-in-1-out node

                            #extend NonBranchingPath by the outgoing edge (w, u) from w

                            w = graph[w][0]

                            paths.append(w)  

                        contigs.append(paths)

                        

        for v in graph:#range(0, len(graph[node])): #length of KEYS in graph

            for w in graph[v]: #range(0, len(num_outgoing)): #for each outgoing edge (v, w) from v

                if num_incoming[v] == 1 or num_outgoing[v] == 1:   

                    cycle = [v]

                while (num_incoming[v] == 1 and num_outgoing[v] == 1 and cycle[0] != cycle[-1]): #w is a 1-in-1-out node

                    #extend NonBranchingPath by the outgoing edge (w, u) from w

                    num_outgoing[v] += -1

                    v = graph[v].pop(0)

                    cycle.append(v)

                    #print("I got here!")

                    #print(cycle)

                    if cycle[0] == cycle[-1]:

                        contigs.append(cycle)

                        break



    #print("contigs :", contigs)

    return contigs





# def calculate_assembly_quality(N, sequences):

#     best_score = 999999999999

#     for k in range(0, len(sequences)):

#         sorted_sequences = sorted(sequences[k], key=len)

#         sorted_sequences.reverse()

#         totallength = 0

#         #print(sorted_sequences)

#         #print(len(sorted_sequences))

#         #calculate total coverage

#         for j in range(0, len(sorted_sequences)):

#             #print("sequence :",sequence)

#             #print(sorted_sequences[j])

#             totallength += len(sorted_sequences[j][0])

            

        

#         #calculate bases in N% coverage

#         percentage = N/100

#         bp_coverage = totallength*percentage

        

#         base_count = 0

#         i = -1

#         while base_count < bp_coverage:

#             i += 1

#             base_count += len(sorted_sequences[i][0])

            

#         #print(sorted_sequences[i][0])

#         nscore = len(sorted_sequences[i][0])

#         print(nscore)

#         if nscore > best_score:

#             best_score = nscore

#         #print(nscore)

#         return best_score



def contig_string_unpaired(path):

    contig_string = path[0]

    for j in range(1, len(path)):

        contig_string += path[j][-1]

    return contig_string



def contig_string(path, k, d):

    prefix = contig_string_unpaired([node[0] for node in path])

    suffix = contig_string_unpaired([node[1] for node in path])

    return prefix + suffix[-(k+d):]



def calculate_assembly_quality(N, sequences, k, d):

    sorted_sequences = sorted([contig_string(path, k, d) for path in sequences], key=len)

    sorted_sequences.reverse()

    totallength = 0

    

    #calculate total coverage

    for sequence in sorted_sequences:

        totallength += len(sequence)



    #calculate bases in N% coverage

    percentage = N/100

    bp_coverage = totallength*percentage

    

    base_count = 0

    i = -1

    while base_count < bp_coverage:

        i += 1

        base_count += len(sorted_sequences[i])

        

    nscore = len(sorted_sequences[i])

    return nscore







def paired_kmer_score(paired_kmers, k, d):

    #paired_kmers = paired_kmer_composition(sequence, k, d)

    paired_graph = construct_paired_debruijn(paired_kmers)

    paired_k_contigs = generate_contigs(paired_graph)

    score = calculate_assembly_quality(DEFAULT_N, paired_k_contigs, k, d)

    return score





def main_paired(sequence, vary_d):

    if vary_d == False:

        row = 0

        for i in range(0, len(paired_k_values)):

            #sheet1.write(row, 0, paired_k_values[i])

            paired_kmers = paired_kmer_composition(sequence, paired_k_values[i], D)

            k_score = paired_kmer_score(paired_kmers, paired_k_values[i], D)

            print(paired_k_values[i],D, k_score)

            #sheet1.write(row, 1, k_score)

            #sheet1.write(row, 2, D)

            row += 1

    else:

        row = 0

        for i in range(0, len(D_values)):

            #sheet1.write(row, 0, D_values[i])

            paired_kmers = paired_kmer_composition(sequence, K, D_values[i])

            #paired_kmers = paired_kmer_composition(sequence, K, D)

            k_score = paired_kmer_score(paired_kmers, K, D_values[i])

            print(K, D_values[i], k_score)

            #sheet1.write(row, 1, k_score)

            #sheet1.write(row, 2, K)

            row += 1





### OUTPUT

#WORKBOOK = xlsxwriter.Workbook('paired_N50comparison__1.xlsx')

#sheet1 = WORKBOOK.add_worksheet()

#main_paired(gene, False)

main_paired((gene), True)

#WORKBOOK.close()



### OUTPUT--READ FILE

# gene = readGN("gene.txt")

# WORKBOOK = xlsxwriter.Workbook('paired_N50comparison__0.xlsx')

# sheet1 = WORKBOOK.add_worksheet()

# sample_output = main_paired(gene)

# WORKBOOK.close()