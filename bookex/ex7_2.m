% example ex7_2
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
Coord=[1    40  10;
	   2    40   0;
	   3    60   0;
	   4    60  10]; 
Coord(:,2:3)=Coord(:,2:3)*1e-3;
%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 4;
       2  2 3 4];
%-----------Element type----------------------------------------
ElTp = 'axisyme';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 f(GDof(1,1))=2514;
 f(GDof(2,1))=2514;
 
%----- Element properties  --------------------------------------
 ptype=3
 E=200e9; v=0.3;
 D=hooke(ptype,E,v);
 
%----- Assemble Ke into K ---------------------------------------
 for i = 1:2
   Ke=axisyme(Ex(i,:),Ey(i,:),D);
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------
 bc= [GDof(1,2) 0; GDof(2,2) 0;
      GDof(3,1) 0; GDof(3,2) 0;
	  GDof(4,1) 0; GDof(4,2) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:2
	 ed(i,:)=extract(EDof(i,:),a);
	 [es(i,:) et(i,:)]=axisyms(Ex(i,:),Ey(i,:),D,ed(i,:));
 end
%---------------------------- end -------------------------------
