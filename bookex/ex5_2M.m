% example ex5_2M: with joint connectivity
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-06
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear;
%-----------Global coordinate matrix-----------------------------
Coord=[1     0   8;
	   2     12  8;
	   3     0   0;
	   4     12  0;
	   5     12  8]; 
Coord(:,2:end)=Coord(:,2:end)*12;

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2;
       2  1 3;
	   3  5 4];
	   
%-----------Element type----------------------------------------
ElTp = 'beam2e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------

[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end));	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=sparse(nDof,nDof);
 f=zeros(nDof,1);
 f(GDof(1,1))=3000;
 eq=[0   -500/12;
     0     0;
	 0     0];
 
%----- Element properties  --------------------------------------
 E=30e6; A=6.8; I=65;
 ep=[E A I];
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:3
   %Element stiffness matrices
   [Ke,fe]=beam2e(Ex(i,:),Ey(i,:),ep,eq(i,:)); 
   [K,f]=assem(EDof(i,:),K,Ke,f,fe);
 end
 
%----- Solve the system of equations ----------------------------
 bc= [GDof(3,1)  0; GDof(3,2)  0; GDof(3,6) 0;
	  GDof(4,1)  0; GDof(4,2)  0; GDof(4,6) 0;
	  GDof(2,1)  0;-GDof(5,1) -1;
	  GDof(2,2)  0;-GDof(5,2) -1];   
 [a,r]=solvep(K,f,bc)

%----- Section forces -------------------------------------------
 for i = 1:3
	Ed(i,:)=extract(EDof(i,:),a);
	es{i}=beam2s(Ex(i,:),Ey(i,:),ep,Ed(i,:),eq(i,:),21); 
 end

 %----- Draw deformed frame ---------------------------------------
 
 figure(1)
 hold on;
 sfac=scalfact2(Ex(1,:),Ey(1,:),Ed(3,:),0.1);
 for i = 1:3
     eldraw2(Ex(i,:),Ey(i,:),[2 1 0]);	 
	 eldisp2(Ex(i,:),Ey(i,:),Ed(i,:),[1 2 1],sfac);
 end
 axis([-25 175 -20 120]);
 pltscalb2(sfac,[1e-1 60 0]);
 title('displacements')
 
%----- Draw normal force diagram --------------------------------
 
 figure(2)
 sfac=scalfact2(Ex(3,:),Ey(3,:),es{3}(:,1),0.2);
 for i = 1:3
     eldia2(Ex(i,:),Ey(i,:),-es{i}(:,1),[2 1],sfac);
 end
 axis([-25 175 -20 120]);
 pltscalb2(sfac,[5e3 60 0]);
 title('normal force')

%----- Draw shear force diagram ---------------------------------
 
 figure(3)
 sfac=scalfact2(Ex(1,:),Ey(1,:),es{1}(:,2),0.2);
 for i = 1:3
     eldia2(Ex(i,:),Ey(i,:),es{i}(:,2),[2 1],sfac);
 end
 axis([-25 175 -20 140]);
 pltscalb2(sfac,[5e3 60 0]);
 title('shear force') 

%----- Draw moment diagram --------------------------------------
 
 figure(4)
 sfac=scalfact2(Ex(1,:),Ey(1,:),es{1}(:,3),0.2);
 for i = 1:3
     eldia2(Ex(i,:),Ey(i,:),es{i}(:,3),[2 1],sfac);
 end
 axis([-40 175 -20 140]);
 pltscalb2(sfac,[1e5 60 0]);
 title('moment') 

%------------------------ end -----------------------------------