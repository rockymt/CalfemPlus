% example exs6_37b
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
Coord=[1     0   0;
	   2    50   0;
	   3    50  50;
	   4     0  50;
	   5    30  10]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 5;
       2  2 3 5;
	   3  3 4 5;
	   4  4 1 5];
%-----------Element type----------------------------------------
ElTp = 'plante';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 
%----- Element properties  --------------------------------------
 ptype=2; thick=10;
 ep=[ptype thick];
 E=1e6; v=0.25;
 D=hooke(ptype,E,v);
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:4
   Ke=plante(Ex(i,:),Ey(i,:),ep,D) %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,1) 0; GDof(1,2) 0;
      GDof(2,1) 5e-3; GDof(2,2) 2.5e-3;
      GDof(3,1) 7.5e-3; GDof(3,2) 7.5e-3;
	  GDof(4,1) 2.5e-3; GDof(4,2) 5.0e-3];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:4
	 ed(i,:)=extract(EDof(i,:),a);
	 [es(i,:) et(i,:)]=plants(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
 end
 
 sfac=scalfact2(Ex(2,:),Ey(2,:),ed(2,:),0.2);
 for i = 1:size(El,1)
     eldraw2(Ex(i,:),Ey(i,:),[2 1 0]);	 
	 eldisp2(Ex(i,:),Ey(i,:),ed(i,:),[1 2 1],sfac);
 end
 %pltscalb2(sfac,[1e-1 60 0]);
 title('displacements')
 a(GDof(5,1))
 a(GDof(5,2))
 es
%---------------------------- end -------------------------------
