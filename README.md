# PhylogeneticTrees
Written in Java, this program builds and manipulates phylogenetic trees by comparing the difference between sets of aligned amino acid sequences, with one sequence for each species. By comparing the differences between these sequences this program can infer a hierarchial structuring of the species with respect to the evolutionary relationships between them. Species live on the leaves of the phylogenetic tree.
Each edge from parent to child has an associated weight that indicates to what extent the child deviates from the proposed parent.
Once the tree is generated, the user can query the trees to find common ancestors, measure evolutionary distance between species and compute basic properties of the tree and its nodes.
This program analyizes amino acid sequences which are in a FASTA format (see https://en.wikipedia.org/wiki/FASTA_format for more information).

