# CalfemPlus
An extension of finite element toolbox CALFEM.
This package "CalfemPlus" is for adding the extended functions and elements for the toolbox of "Calfem" (Computer Aided Learning of the Finite Element Method).
"Calfem" Websites and downloads:
[https://www.byggmek.lth.se/english/calfem/]
[https://github.com/CALFEM/calfem-matlab]

Now the package "CalfemPlus" is workable with "Calfem" version 3.4. Some usage examples are in the directory "bookex", most of which implementing the examples and problems in book "Chandrupatla TR, Belegundu AD. Introduction to Finite Elements in Engineering, 4th edition".

## A short Functions introduction

### ExtractElInfo module functions

This module mainly contains the following functions: ExtractElInfo, DofRaw and CoordElNd. The core function is ExtractElInfo. The call format is:

`[EDof,GDof,Ex,Ey,Ez]=ExtractElInfo(El,ElTp,Coord)`

The input of the function are the element-node connection relation matrix El, element type ElTp, node coordinate table Coord. And return the element-dof connection matrix(topology matrix) EDof, the node-dof connection matrix GDof, and the element coordinate matrix Ex, Ey, Ez. See the application example ex4_1B.m in the directory "bookex" for details.

The function of DofRaw is to set the connection relation table between nodes and the overall degree of freedom for each element type. If the user defines a new element and calls ExtractElInfo, the corresponding information of the element needs to be added to this function. The function of CoordElNd is to extract coordinate vectors of x, y and z of a certain group of elements. This module requires nodes and cells to be numbered consecutively starting with the number 1.

### extended elements

The list of elements supplemented by CalfemPlus is shown in Table 1. 

*Table 1 Extended elements*

|Element function 	|Internal force or gradient function |Element explanation|
|-------------------|------------------------------------|-------------------|
|quad1de 			|quad1ds 	|one-dimensional quadratic element                   |
|flw2i6e   			|flw2i6s 	|planar six-node triangular heat transfer element    |
|plani6e   			|plani6s 	|planar six-node triangular solid element            |
|tetra4e 			|tetra4s 	|space tetrahedral element with four nodes    |
|flw3tet10e 		|flw3tet10s |space ten-node tetrahedron heat transfer element    |
|tetra10e 			|tetra10s 	|space tetrahedral element with ten nodes     |
|flw3i20e 			|flw3i20s 	|space twenty node hexahedron heat transfer element  |
|soli20e 			|soli20s 	|space twenty node hexahedron solid element         |
|axisyme 			|axisyms 	| axisysymmetrical three-node triangle element|
|tin1e 				|- 			|one-dimensional heat dissipation element            |

### other tools
Additional functions in CalfemPlus are listed in Table 2.
Table 2 Additional functions
|Function Name |Function Description |Adaptable element|
|-------------------|------------------------------------|-------------------|
|solvep 	|penalty method solving MPC constraints 						|-               |
|assemf 	|Assembel fe into the global force vector f 									|-               |
|bar2theta 	|bar element temperature load generating function 								|bar2e           |
|EYp2EZo 	|Calculate the unit vector along local z axis from the reference point |beam3e          |
|plant2T2f 	|Create and assemble traction force										|plante, plani4e |
|planttheta |triangle element temperature load generating function 							|plante          |
|plantefit 	|nodal values from known constant triangle element values	|plante, flw2te  |

