% example exs8_27
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
Coord=[1     0   0;
	   2    50   0;
	   3   100   0;
	   4   100  50;
	   5   100 100;
	   6    50 100;
	   7     0 100;
	   8     0  50;
	   9    60  20]; 
 Coord(:,2:end)=Coord(:,2:end)*1e-3;
%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 9 8;
       2  2 3 4 9;
	   3  9 4 5 6;
	   4  8 9 6 7];
%-----------Element type----------------------------------------
ElTp = 'plani4e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 f(GDof(3,1))=2500;
 f(GDof(4,1))=5000;
 f(GDof(5,1))=2500;
 
%----- Element properties  --------------------------------------
 ptype=1; thick=10e-3;
 ep=[ptype thick 2];
 E=200e9; v=0.25;
 D=hooke(ptype,E,v);
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:4
   Ke=plani4e(Ex(i,:),Ey(i,:),ep,D); %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,1) 0; GDof(1,2) 0;
	  GDof(8,1) 0; GDof(7,1) 0];   
 [a,r]=solveq(K,f,bc);

%----- Element forces -------------------------------------------
 for i = 1:4
	 ed(i,:)=extract(EDof(i,:),a);
	 [es{i} et{i}]=plani4s(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
 end
 
 sfac=scalfact2(Ex(2,:),Ey(2,:),ed(2,:),0.2);
 for i = 1:size(El,1)
     eldraw2(Ex(i,:),Ey(i,:),[2 1 0]);	 
	 eldisp2(Ex(i,:),Ey(i,:),ed(i,:),[1 2 1],sfac);
 end
 %pltscalb2(sfac,[1e-1 60 0]);
 title('displacements')
 r(GDof(1,1:2))
 r(GDof(7,1))
 r(GDof(8,1))
%---------------------------- end -------------------------------
