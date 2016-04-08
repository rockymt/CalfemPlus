% example ex5_1
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-06
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------

%-----------Global coordinate matrix-----------------------------
Coord=[1     0   0;
	   2     1   0;
	   3     2   0]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2;
       2  2 3];
	   
%-----------Element type----------------------------------------
ElTp = 'beam2e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------

[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof))
 K=zeros(nDof,nDof)
 f=zeros(nDof,1); 
 eq=[0 0; 0 -12e3];
 
%----- Element properties  --------------------------------------
 E=200e9; A=1e10; I=4e-6;
 ep=[E A I];
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:2
   [Ke,fe]=beam2e(Ex(i,:),Ey(i,:),ep,eq(i,:)) %Element stiffness matrices
   [K,f]=assem(EDof(i,:),K,Ke,f,fe);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,1) 0; GDof(1,2) 0; GDof(1,6) 0;
      GDof(2,2) 0;
	  GDof(3,2) 0];   
 [a,r]=solveq(K,f,bc)


%----- Section forces -------------------------------------------
 for i = 1:2
	Ed(i,:)=extract(EDof(i,:),a);
	es1=beam2s(Ex(i,:),Ey(i,:),ep,Ed(i,:),eq(i,:),21); 
 end
 %----- Draw deformed frame ---------------------------------------
 
 figure(1)
 plotpar=[2 1 0];
 for i = 1:2
	eldraw2(Ex(i,:),Ey(i,:),plotpar);
 end

 sfac=scalfact2(Ex(1,:),Ey(1,:),Ed(1,:),0.1);
 plotpar=[1 2 1];
 for i = 1:2
	eldisp2(Ex(i,:),Ey(i,:),Ed(i,:),plotpar,sfac);
 end

 axis([-1.5 7.5 -0.5 5.5]); 
 pltscalb2(sfac,[1e-2 0.5 0]);
 axis([-1.5 7.5 -0.5 5.5]);
 title('displacements')
 
%----- Draw normal force diagram --------------------------------
 
 figure(2)
 plotpar=[2 1];
 sfac=scalfact2(ex1,ey1,es1(:,1),0.2);
 eldia2(ex1,ey1,es1(:,1),plotpar,sfac);
 eldia2(ex2,ey2,es2(:,1),plotpar,sfac);
 eldia2(ex3,ey3,es3(:,1),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]);
 pltscalb2(sfac,[3e4 1.5 0]);
 title('normal force')

%----- Draw shear force diagram ---------------------------------
 
 figure(3)
 plotpar=[2 1];
 sfac=scalfact2(ex3,ey3,es3(:,2),0.2);
 eldia2(ex1,ey1,es1(:,2),plotpar,sfac);
 eldia2(ex2,ey2,es2(:,2),plotpar,sfac);
 eldia2(ex3,ey3,es3(:,2),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]);
 pltscalb2(sfac,[3e4 0.5 0]);
 title('shear force') 

%----- Draw moment diagram --------------------------------------
 
 figure(4)
 plotpar=[2 1];
 sfac=scalfact2(ex3,ey3,es3(:,3),0.2);
 eldia2(ex1,ey1,es1(:,3),plotpar,sfac);
 eldia2(ex2,ey2,es2(:,3),plotpar,sfac);
 eldia2(ex3,ey3,es3(:,3),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]);
 pltscalb2(sfac,[3e4 0.5 0]);
 title('moment') 

%------------------------ end -----------------------------------
 echo off