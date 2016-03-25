# calfemPlus
An extension of finite element package CALFEM3.4

Contents:

 function [d,Q]=solvep(K,f,bc) ------ Solve static FE-equations considering boundary conditions with penalty method.This function can handle multipoint Constraints.
 
 function [Ke]=quad1De(ep) ------ Compute element stiffness matrix for quadratic 1D element.
 
 function [es]=quad1Ds(ep,ed,xi) ------ Compute element force quadratic 1D element (quad1De).
