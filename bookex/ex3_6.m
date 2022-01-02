 % example ex3_6
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
 
%----- Topology matrix Edof -------------------------------------
 Edof=[1   3  1;
       2   4  2];
%----- Stiffness matrix K and load vector f ---------------------
 K=zeros(5,5); 
 f=zeros(5,1);	f(5)=30e3;
%----- Element stiffness matrices  ------------------------------
 E1=2.0e11;  A1=1.2e-3;  l1=4.5;
 E2=0.7e11;  A2=0.9e-3;  l2=3.0;
 ep1 = E1*A1/l1;  ep2 = E2*A2/l2;
 Ke1=spring1e(ep1)
 Ke2=spring1e(ep2)
%----- Assemble Ke into K ---------------------------------------
 K=assem(Edof(1,:),K,Ke1)
 K=assem(Edof(2,:),K,Ke2) 
%----- Solve the system of equations ----------------------------
 bc = [3 0; 
	   4 0; 
	   1 0; -5 -0.333333; 
	   2 0; -5 -0.833333];
 [a,r]=solvep(K,f,bc)
 %----- Element forces -------------------------------------------
 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)
 es1=spring1s(ep1,ed1)/A1
 es2=spring1s(ep2,ed2)/A2
 %---------------------------- end -------------------------------
