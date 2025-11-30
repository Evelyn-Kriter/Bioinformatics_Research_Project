from collections import Counter

#import xlsxwriter





### TEST RUN DEFAULTS AND SAMPLES

DEFAULT_N = 50



#k_values = [2, 5, 10, 20, 40, 80, 160, 200, 250, 300, 320, 350, 400, 450, 500, 550, 600, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760, 774]

#k_values = [500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000]

k_values = [2, 4, 6, 12, 16, 17, 18, 19, 20, 160, 400, 600, 800, 1100]



#Carsonella Ruddii Gene--has length 774

sample_gene = "TTAACAATTTTTCCATAGTTTTAAAAAATGAAAGAATTTTTTATAACCTAAATTAAAAACACTAGATACCATAAAATAAGGTGTTATTATTAAATTTATTTTTTTTATTTTACTAAAAAAAAAAACATTTTCAATTCCATCTATTATAATCCAATAATTTATATTATTTAAATTTTTAATATATAAAAAATATAAAAAAAATTTTATTAATTTAAAAAATATTTTAAAAAATAATAAATTTATATAAAAAATATAAATTTTAGAAATAAAATTATTTTTATAATTAAAATTGTTTAATAAAAAATTTTCAGAAATATTAATATAACATTTATAATAAAAATATTTGAAATTAATTTTTATAAAAATAATTTTACCTTGTTTTCCAAAATTAGCAATTTTTGAATTATAAAAATTTTTATAACAAAAATTTCCATTACCTCCTTTTCCTCCTTTTAATATTTTAATAAAATAATTGTTTTTAACAATATAAATTTTTTTTTTATTTAATGTTATTATTACACCTAAAGATATTTTTATTATTAAATCTTTACCTTTTTTTCCTTTTTTGTTATTTTTTTCACCATTTTTACCATTTTTTGAAATAAATAAATTATTATTTGGTAATTTAAAATTATTATTAGATAATAAATATAAATCACCACCATCACCACCATCACCTCCATTAGCATAAAATTTTCCATTTAAAAATTTATAACTAATACAACCATTTCCTCCTTTACCACTTTTTAAAAGTAAATATTTTTCAATAAACAT"



gene_1 = "TTAACAATTTTTCCATAGTTTTAAAAAATGAAAGAATTTTTTATAACCTAAATTAAAAACACTAGATACCATAAAATAAGGTGTTATTATTAAATTTATTTTTTTTATTTTACTAAAAAAAAAAACATTTTCAATTCCATCTATTATAATCCAATAATTTATATTATTTAAATTTTTAATATATAAAAAATATAAAAAAAATTTTATTAATTTAAAAAATATTTTAAAAAATAATAAATTTATATAAAAAATATAAATTTTAGAAATAAAATTATTTTTATAATTAAAATTGTTTAATAAAAAATTTTCAGAAATATTAATATAACATTTATAATAAAAATATTTGAAATTAATTTTTATAAAAATAATTTTACCTTGTTTTCCAAAATTAGCAATTTTTGAATTATAAAAATTTTTATAACAAAAATTTCCATTACCTCCTTTTCCTCCTTTTAATATTTTAATAAAATAATTGTTTTTAACAATATAAATTTTTTTTTTATTTAATGTTATTATTACACCTAAAGATATTTTTATTATTAAATCTTTACCTTTTTTTCCTTTTTTGTTATTTTTTTCACCATTTTTACCATTTTTTGAAATAAATAAATTATTATTTGGTAATTTAAAATTATTATTAGATAATAAATATAAATCACCACCATCACCACCATCACCTCCATTAGCATAAAATTTTCCATTTAAAAATTTATAACTAATACAACCATTTCCTCCTTTACCACTTTTTAAAAGTAAATATTTTTCAATAAACAT"

gene_2 = "TATTAACTTCTATAAAAGTTTTATTTTTTTTTTTTTTAATTAATACTATTCCATTTATTAATGACTGTATATAATAATTTTTTGAATAAAAAGTATTTTTTCCTGAAATTATTTTAGAACCATTTTGTTTAATTATTATAGATCCTTTATTAACATAATTATTATTATAAATTTTTAATCCTAATCTTTTAGAAATAGAATCTCTTCCATTTCTCGTACTACCTCCAGCTTTTTTTTGTGCCATATTAAATTATATTTTGTAAAAATAATAAACTTTTTTTTTTATAATTTATATTAGTTTTTAAAAAATTTTTTCTTCTTTTTTTTTTTATAGATATAGATTTTATATAAAAATGTTTTAAAATTATAAAACATGCTTTTAAATTTAAACCTTTTTTATCTAAAATTATTTTTTTATTAAAAAAAAAAAATATTTTTTTAAAAAAAATTTTTTTACCAATATTTTTATTTATATAATCAACAATTAAAAAATTATTTAATTTTGAAAAATAAATTTTTTTACCTATACAAAAAGATAACAAATTATTTAATTTTATCATTTTTTAAAATATTATTATTATAATAATTAAAAACACCAAGTTCACCAATTAAAGAACAATCACCTTTTCTATATCCTATTTTTATTAAATAAAAAAAAATTTTTTTTTTAAAATAATTATTTAAAAGAATAAAAATTTTAATTTTATTTTTTTTTAATTTAAAATAATTACTTTTAATTTTACCAAATTTTATTAAATTAAAAAAA"

gene_3 = "TTTTTTAAAGCATTTGATATTTCTAAAAATACAAATTTATTCATATTTTTTAATTTAAATAAATCAAATTCTGATAATTTAATCAAGTCTCCTATTAAAAAAATGTTGTTTCTTTTTAATATATTTGAAGTTTTTATACTTAATTCTAAATTATCTACTGATCTAATAAATATAGGATTAATACTTAAAAAATTTTTTTTTTTAATTTTTTTTGTTTTTTTTTTTCCTAAAACTGAAAAAAAAACATCAAAATATTTTTTTATATAAAATATACAATTATTAAAACATTCTTCAGGTTTAATAATTCCGTCAGTTTCAATATCAAAAAATAAATTTTTAATTTTTTTATTAAAAAATTTTTTATGTATAAAATAATTTATGTTTTTTAAAGAAGACTTTAAAAAATTAATTTTAATAACTTTTGAATTAAAAATTTTATTAATAATATTATTATTAGTAAAATTATAACTAGTATTAATACATTTCATTAAAACAAAAAAAATTATATTATCTGTTATATTTGCTATTATTATATCTGGATTATAAATAATAATATGTTTATCAGAAAAAATATCTTTTGCTTTAACAATACAAGGTCCTTTTTTTTTTATAATTAAAAAAGCAGTTATTGAATTTTTAATTTTTATTAAAATATTATTTATATTATTAATAATTTTTAAAGTATTTTCTTTTATACCTTTTATATTAGAAAATTCAGAA"

gene = gene_1+gene_2+gene_3



#Test genes

sample_gene1 = "TTATCCATGTGTGT"

sample_gene2 = "TTAACAATTTTTCCATAGTTTTAAAAAATGAAAGAATTTTTTATAACCTA" #sample of genome[:30]









### CODE START###



### Read a genome

def readGN(filename):

    with open(filename, "r") as file:

        sequence = [line.strip() for line in file.readlines()]

        return sequence[0]



def kmer_composition(sequence, k):

    #Return set of unique k-mers in sequence

    #sequence is a string

    #k is an int

    unique_kmers = []

    for i in range(0, len(sequence)-k+1):

        unique_kmers.append(sequence[i:i+k])

    #unique_kmers.sort()

    return unique_kmers





def construct_debruijn(kmers):

    pn_kmers = kmers

    pn_kmers.sort()

    

    graph = {}

    for i in range(0, len(pn_kmers)):

        

        #turn every prefix into a key in the dictionary and add suffix

        prefix = pn_kmers[i][:-1]

        suffix = pn_kmers[i][1:]

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

                if num_incoming[v] == 1 and num_outgoing[v] == 1:   

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



    return contigs



def contig_string(path):

    contig_string = path[0]

    for j in range(1, len(path)):

        contig_string += path[j][-1]

    return contig_string





def calculate_assembly_quality(N, sequences):

    sorted_sequences = sorted([contig_string(path) for path in sequences], key=len)

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



def kmer_score(sequence, k):

    kmers = kmer_composition(sequence, k)

    graph = construct_debruijn(kmers)

    k_contigs = generate_contigs(graph)

    score = calculate_assembly_quality(DEFAULT_N, k_contigs)

    return score





def main_unpaired(sequence):

    #score_list = []

    row = 0

    for i in range(0, len(k_values)):

        #sheet1.write(row, 0, k_values[i])

        k_score = kmer_score(sequence, k_values[i])

        print(k_values[i], k_score)

        #sheet1.write(row, 1, k_score)

        row += 1

    



### OUTPUT

#WORKBOOK = xlsxwriter.Workbook('unpaired_N50comparison__0.xlsx')

#sheet1 = WORKBOOK.add_worksheet()

main_unpaired(gene)

#WORKBOOK.close()



### OUTPUT--READ FILE

# genome = readGN("genome.txt")

# WORKBOOK = xlsxwriter.Workbook('N50comparison_genome.xlsx')

# sheet1 = WORKBOOK.add_worksheet()

# sample_output = main_unpaired(genome)

# WORKBOOK.close()



