% example ex6_7
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-20
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
 
%-----------Global coordinate matrix-----------------------------
Coord=[1     3   0;
	   2     3   2;
	   3     0   2;
	   4     0   0]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 4;
       2  3 4 2];
%-----------Element type----------------------------------------
ElTp = 'plante';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 f(GDof(2,2))=-1000
 
%----- Element properties  --------------------------------------
 ptype=1; thick=0.5;
 ep=[ptype thick];
 E=30e6; v=0.25;
 D=hooke(ptype,E,v);
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:2
   Ke=plante(Ex(i,:),Ey(i,:),ep,D) %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(3,1) 0; GDof(3,2) 0;
      GDof(1,2) 0;
	  GDof(4,1) 0; GDof(4,2) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:2
	 ed(i,:)=extract(EDof(i,:),a);
	 [es(i,:) et(i,:)]=plants(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
 end
 r(GDof(1,2))
 r(GDof(3,1:2))
 r(GDof(4,1:2))
%---------------------------- end -------------------------------
