# Bioinformatics-Research-Project
# Abstract
### Independent research N50 scores for De Bruijn Graph assembly. I generated paired and unpaired De Bruijn graphs from sequencing reads with varying k and d values
### for three genes from Carsonella Ruddii. Using a contig generator, I extended non-branching paths in 1-in-1-out nodes and calculated N50 scores across parameter 
### ranges. Results interestingly showed spikes in assembly quality before linear decline, which may suggest a possible mathematical formula that could predict 
### lowest optimal k and d values. 
# 
# Background
### Unpaired Debruijn graphs provide a better approach to assembling reads compared to algorithms that try every combination, though doing so sacrifices ability to 
### solve repeats well. Runtime can also be an issue as well. The Paired Debruijn graph was introduced as an alternative, using kmer pairs of length k separated by 
### bases of length d. For Debruijn graphs in general, higher k values improve ability to resolve repeats, but N50 score is reduced the likelihood an error in the 
### read in increased. Lower k values increase N50 score and decreases likelihood of errors in reads but ability to solve repeats is decreased. Lower k values also 
### take longer to run because more kmers are generated. But with Paired Debruijn graphs, contig sizes are improved. (Medvedev et. al., 2011). 
### N scores can be used to measure contig coverage. For N50, half of the genome is covered by contigs larger than or equal to the N50 score.
# 
# Goal
### The goal of my project is to obtain N50 scores for a range of k values for Unpaired Debruijn Graphs and a range of k and d values for Paired Debruijn graphs.
# 
# Algorithm
### Paired and Unpaired Debruijn graphs were created from reads with varying d and k values. A list of contigs was generated using the code on the left from each
### Debruijn graph type. The contig generator finds a non branching path and extends it for as long as each node is 1-in-1-out. The N50 score was calculated for a 
### series of k and d values. These scores went into a list that was exported to excel to be visualized. 
#
# Experiment 
### Error free k-mers were generated with perfect coverage from three different subsets of a reference genome: one gene at 774bp, two genes with total 1540bp, or 
### three genes with total 2260bp (NIH, 2017).
### For Unpaired Debruijn graphs, Figure 1 confirms that N50 scores increase as k-values increase and then plateaus after a certain point, as expected from error-### free reads. For Paired Debruijn graphs, Figures 2, 3, and 4 show that increasing k by one yields either the same or better N50 score before reaching a plateau, 
### modeling a similar trend to varying k for Unpaired Graphs in Figure 1. Increasing d while keeping k constant shows a general increase in N50 score with spikes of
### high assembly quality before decreasing linearly as d approaches the sequence length. There are many small fluctuations in the raw data when k is constant and d
### is varied that could be elucidated with more time and processing power. Future research might extend this experiment by predicting the cause of any correlations
### in the data. 
