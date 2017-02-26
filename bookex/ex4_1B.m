% example ex4_1B
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
 
%-----------Global coordinate matrix-----------------------------
Coord=[1     0   0;
	   2    40   0;
	   3    40  30;
	   4     0  30]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2;
       2  3 2;
	   3  1 3;
	   4  4 3];
	   
%-----------Element type----------------------------------------
ElTp = 'bar2e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------

[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof))
 K=zeros(nDof,nDof)
 f=zeros(nDof,1); 
 f(GDof(3,2))=-25000; f(GDof(2,1))=20000
 
%----- Element properties  --------------------------------------
 E=29.5e6; A=1.0;
 ep=[E A]
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:4
   Ke=bar2e(Ex(i,:),Ey(i,:),ep) %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,1) 0; GDof(1,2) 0;
      GDof(2,2) 0;
	  GDof(4,1) 0; GDof(4,2) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:4
	 ed(i,:)=extract(EDof(i,:),a);
	 N(i)=bar2s(Ex(i,:),Ey(i,:),ep,ed(i,:));
	 sigma(i)=N(i)/A;
 end
%---------------------------- end -------------------------------
