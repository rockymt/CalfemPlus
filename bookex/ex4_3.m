% example ex4_3
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% 龙驭球 包世华
% 结构力学教程
%----------------------------------------------------------------
clear 
%-----------Global coordinate matrix-----------------------------
Coord=[1  -2   0   0;
	   2   0   0   4;
	   3   2   0   0;
	   4   0  -3   0]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2;
       2  2 3;
	   3  2 4];
	   
%-----------Element type----------------------------------------
ElTp = 'bar3e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------

[EDof,GDof,Ex,Ey,Ez]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof))
% K=zeros(nDof,nDof)
 K=sparse(nDof,nDof)
 f=zeros(nDof,1); 
 f(GDof(2,2))=-3000
 
%----- Element properties  --------------------------------------
 E=200e9; A=900e-6;
 ep=[E A]
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:3
   Ke=bar3e(Ex(i,:),Ey(i,:),Ez(i,:),ep); %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,1) 0; GDof(1,2) 0; GDof(1,3) 0;
      GDof(3,1) 0; GDof(3,2) 0; GDof(3,3) 0;
	  GDof(4,1) 0; GDof(4,2) 0; GDof(4,3) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:3
	 ed(i,:)=extract(EDof(i,:),a);
	 N(i)=bar3s(Ex(i,:),Ey(i,:),Ez(i,:),ep,ed(i,:));
	 sigma(i)=N(i)/A;
 end
%---------------------------- end -------------------------------
