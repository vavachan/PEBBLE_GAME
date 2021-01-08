## PEBBLE GAME 
This code uses pebble game to identify rigid regions on a given graph in 2D.
The graph is assumed to be planar. In this work each node (vertex) in the graph 
has two degrees of freedom (translations in vertical and horizontal directions)
and the graph has three globals degress of freedom ( two translations and one 
rotation). A graph is rigid if we cannot apply displacements to the vertices 
without changing the length of the edges in the graph. 

Please refer to the following paper for details of the algorithm. 
[Jacobs, Donald J., and Bruce Hendrickson. "An algorithm for two-dimensional rigidity percolation: the pebble game." Journal of Computational Physics 137.2 (1997): 346-365.][1]

[1]: https://people.engr.tamu.edu/ajiang/PebbleGame.pdf
