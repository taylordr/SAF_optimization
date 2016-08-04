# SAF_optimization
This code implements the optimization of phase synchronization for an undirected chain in which the 
oscillators are drawn from a normal distribution. For any usage of this code, please cite:
D Taylor, PS Skardal, J Sun. (2016) ''Synchronization of heterogeneous oscillators under network 
modifications: Perturbation and optimization of the synchrony alignment function.'' 
arXiv preprint:  http://arxiv.org/abs/1605.04009

This folder includes the following files:

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

If you have any questions about the code, please email Dane Taylor at dane.r.taylor@gmail.com.
