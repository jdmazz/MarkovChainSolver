# MarkovChainSolver
@author J.D. Mazz

This program solves for the steady state of any arbitrary probability graph.  This program can handle cyclic Digraphs with/without terminal nodes and  determine the probability distributions of all nodes within the network. Takes a matrix of counts, the desired level of convergence (necessary for cycles),  and the max precision of the denominator of the fractional solution and returns The counts indicate the ratios of the transtions and the column sums will become the denominators of the transition matrix. Each row of the transition matrix will sum to one. Once the steady state is achieved, the row vector solution  will be normalized and indicates the relative frequency of each node.

Ex: [terminal_node_numerator1, ..., common_denominator]
