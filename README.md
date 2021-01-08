## PEBBLE GAME 
This code uses pebble game to identify rigid regions on a given graph in 2D.
The graph is assumed to be planar. In this work each node (vertex) in the graph 
has two degrees of freedom (translations in vertical and horizontal directions)
and the graph has three globals degress of freedom ( two translations and one 
rotation). A graph is rigid if we cannot apply displacements to the vertices 
without changing the length of the edges in the graph. 

Please refer to the following paper for details of the algorithm.
 
[Jacobs, Donald J., and Bruce Hendrickson. "An algorithm for two-dimensional rigidity percolation: the pebble game." Journal of Computational Physics 137.2 (1997): 346-365.][1]

[This lecture from MIT Opencourseware by Dr Erik Demaine gives an excellent explanation of the algorithm] [2]

The input file for the code is the edge list of the graph. Note that the code does not require the cartesian co-ordinate of the edges, just the indices of the vertices connected by an 
edge. For vertices listed from 0 to N_v, if there is a edge between i and j, the input file should contain a column with  i j. To run the code use 

`python 2d_pebble.py [edge_list] [N_v]`

[1]: https://people.engr.tamu.edu/ajiang/PebbleGame.pdf
[2]: https://www.youtube.com/watch?v=yvatNaV6Bog
