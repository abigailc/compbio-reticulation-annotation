I had downloaded all the sequences for a specific gene available on NCBI, built a massive phylogenetic tree, broke it into subtrees using a taxonomy-based algorithm, and then I wanted to analyse what horizontal gene transfer events may have occured within each sub-tree.

For each sub-tree dataset, I generated a gene tree based on the SOD sequences using RAxML with 100 bootstraps. I also generated a species-tree using 30 concatenated ribosomal proteins, also run through RAxML. 

Here I have created a program that takes in a folder containing my species trees and a folder containing my gene trees, and then 

1) creates input files suitable for rangerDTL 
	100 input files per dataset, to capture all information across bootstraps
	this involved changing the tip names in each tree so the format of the species tree and of the gene tree matched, as well as removing underscores, and writing each to a new file.

2) runs each input file through OptRoot 
	OptRoot is part of RangerDTL package, used to optimally root the gene tree. 
	my species trees were all rooted at an appropriate outgroup species.

3) runs each input file through RangerDTL
	saves the result, and replaces the underscores where they should be in the tip names

4) parses the resulting output files, and scores for each node of each bootstrapped gene tree / best species tree whether that node was found to be a transfer event, duplication event, or speciation event (losses are not really inferrable from these datasets, but I left in support so add them in the future)

5) uses #4 to write to file four ways of interpreting the results

	a) the bootstraps file is annotated at each node with what RangerDTL classified it as (transfer/duplication/speciation)
		also nodes downsteam of a transfer event are annotated with either an R(recipient) or a D(donor) if polarity can be inferred by comparison to the species tree. data_ANNboots.newick

	b) the best species tree file is annotated at each node with the number of times RangerDTL classified it as a (transfer donor/transfer recipient/duplication/speciation). nodes in this sense are each defined by their downstream tips, regardless of the topology of internal branches. That way slightly different topologies across bootstraps are counted correctly at deeper nodes. data_ANNspecies.newick

	c) for the best species tree, a tab-separated file is generated containing every sub-clade and the number of times it was classified in each way by RangerDTL. data_master_species_spreadsheet

	d) for all gene trees combined, a tab-separated file is generated containing every sub-clade and the number of times it was classified in each way by RangerDTL. data_master_gene_spreadsheet

Anyways, here's how to get it to work locally

1. clone the repository
2. make sure you have OptRoot.linux and RangerDTL.linux downloaded and in the same folder as the code (or executable in $PATH). Downloadable at http://compbio.engr.uconn.edu/software/RANGER-DTL/
3. run $  python annotate_nodes_main.py test_data/ -p SOD_new -b gene_boots/ -s species_best/
4. watch as it takes some time to run (the slow step is RangerDTL actually computing the reticulations, at least on my computer) - I've uploaded some of the smaller datasets so it shouldn't take too long.
5. view the output in the results folder. the trees should open nicely in FigTree (available at http://tree.bio.ed.ac.uk/software/figtree/)
