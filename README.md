# SAF_optimization
Implements the optimization of synchronization for networks of heterogeneous phase oscillators 
with natural frequencies drawn from a normal distribution. For any usage of this code, please cite:

D Taylor, PS Skardal, J Sun. (2016) ''Synchronization of heterogeneous oscillators under network modifications: Perturbation and optimization of the synchrony alignment function.'' SIAM Journal on Applied Mathematics 76(5), 1984-2008
https://epubs.siam.org/doi/abs/10.1137/16M1075181

The following files are included:

1. Demo1_chain_network.m - a demo of the included scripts for a chain
2. Demo2_SF_network.m - a demo of the included scripts for a power-law network
3. create_chain.m - makes a chain network
4. create_SF.m - makes a power law network network
5. compute_SAF.m - compute the synchrony alignment function for a given w and L	
6. compute_Q_matrix - computes the first order approximation Q_{pq}
7. rank_edges.m - ranks potential edges according to Q_{pq}
8. add_edge.m - adds a newedge to a network
9. change_in_SAF_under_edge_addition.m - computes the change to the SAF
10. algorithm_6_1.m - implements Algorithm 6.1
11. algorithm_6_2.m - implements Algorithm 6.2

If you have any questions about the code, please email Dane Taylor at danet@buffalo.edu.
