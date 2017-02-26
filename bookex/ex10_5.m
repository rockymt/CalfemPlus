% example ex10_5
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-06-02
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
 
%-----------Global coordinate matrix-----------------------------
Coord=[1     0   0;
	   2     2 1.5;
	   3     4   0;
	   4     4   3;
	   5     0   3];
Coord(:,2:end)=1e-2*Coord(:,2:end);
%-----------Element connectivity matrix--------------------------	   
El  = [1  1 3 2;
       2  3 4 2;
	   3  4 5 2;
	   4  5 1 2];
%-----------Element type----------------------------------------
ElTp = 'flw2te';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 
%----- Element properties  --------------------------------------
 ep=1;
 D=[1 0;
    0 1];
 eq=2;
%----- Assemble Ke into K ---------------------------------------
 for i = 1:4
   [Ke fe]=flw2te(Ex(i,:),Ey(i,:),ep,D,eq);
   [K f]=assem(EDof(i,:),K,Ke,f,fe);
 end
 
%----- Solve the system of equations ----------------------------
 bc= [GDof(3,1) 0; GDof(4,1) 0; GDof(5,1) 0];   
 [a,r]=solveq(K,f,bc)

%----- angle of twist -------------------------------------------
 alpha = 0;
 for i = 1:4
	 ed=extract(EDof(i,:),a);
	 x13=Ex(i,1)-Ex(i,3); x23=Ex(i,2)-Ex(i,3);
	 y13=Ey(i,1)-Ey(i,3); y23=Ey(i,2)-Ey(i,3);
	 Area=0.5*(x13*y23-x23*y13);
	 alpha = alpha+sum(ed)*Area*8/3;
 end
 alpha=1/alpha
%---------------------------- end -------------------------------
