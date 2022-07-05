% example exs8_11
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-20
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
% Problems
%----------------------------------------------------------------
clear 
%-----------Global coordinate matrix-----------------------------
Coord=[1     0 -0.4;
	   2     2 -0.4;
	   3     2  0.4;
	   4     0  0.4;
	   5     1 -0.4;
	   6     2    0;
	   7     1  0.4;
	   8     0    0;
	   9     4 -0.4;
	   10    4  0.4;
	   11    3 -0.4;
	   12    4    0;
	   13    3  0.4;
	   14    6 -0.4;
	   15    6  0.4;
	   16    5 -0.4;
	   17    6    0;
	   18    5  0.4]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 3 4 5 6 7 8;
       2  2 9 10 3 11 12 13 6;
	   3  9 14 15 10 16 17 18 12];
%-----------Element type----------------------------------------
ElTp = 'plani8e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 f(GDof(15,2))=-10000;
 
%----- Element properties  --------------------------------------
 ptype=1; thick=1;
 ep=[ptype thick 3];
 E=30e6; v=0.3;
 D=hooke(ptype,E,v);
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:3
   Ke=plani8e(Ex(i,:),Ey(i,:),ep,D); %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(8,1) 0; GDof(8,2) 0;
	  GDof(1,1) 0; GDof(4,1) 0];   
 [a,r]=solveq(K,f,bc);

%----- Element forces -------------------------------------------
 for i = 1:3
	 ed(i,:)=extract(EDof(i,:),a);
	 [es{i} et{i} eci{i}]=plani8s(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
 end
 
 sfac=scalfact2(Ex(2,1:4),Ey(2,1:4),ed(2,1:8),0.2);
 for i = 1:3
     eldraw2(Ex(i,1:4),Ey(i,1:4),[2 1 0]);	 
	 eldisp2(Ex(i,1:4),Ey(i,1:4),ed(i,1:8),[1 2 1],sfac);
 end
 %pltscalb2(sfac,[1e-1 60 0]);
 title('displacements')
 a(GDof(17,2))
 
%---------------------------- end -------------------------------
