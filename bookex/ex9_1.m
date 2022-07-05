% example ex9_1
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-05-24
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear
%-----------Global coordinate matrix-----------------------------
Coord=[1   0   1   1;
	   2   0   0   1;
	   3   1   0   1;
	   4   0   0   0]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 3 4];
	   
%-----------Element type----------------------------------------
ElTp = 'tetra4e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------

[EDof,GDof,Ex,Ey,Ez]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof))
% K=zeros(nDof,nDof)
 K=sparse(nDof,nDof)
 f=zeros(nDof,1); 
 f(GDof(1,3))=-1000;
 
%----- Element properties  --------------------------------------
 E=30e6; v=0.3; ptype=4;
 D=hooke(ptype,E,v);
 
%----- Assemble Ke into K ---------------------------------------
Ke=tetra4e(Ex,Ey,Ez,D); %Element stiffness matrices
[K]=assem(EDof,K,Ke);
 
%----- Solve the system of equations ----------------------------
 bc= [GDof(2,1) 0; GDof(2,2) 0; GDof(2,3) 0;
      GDof(3,1) 0; GDof(3,2) 0; GDof(3,3) 0;
	  GDof(4,1) 0; GDof(4,2) 0; GDof(4,3) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
ed=extract(EDof,a);
[es,et]=tetra4s(Ex,Ey,Ez,D,ed)
%---------------------------- end -------------------------------

