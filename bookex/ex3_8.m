% example ex3_8
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

 Edof=[1 1 2;
       2 2 3];
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(3,3) 
 f=zeros(3,1);  f(2)=300e3
 
%----- Element stiffness matrices  ------------------------------
 E1=70e9; A1=900e-6; L1=0.2; a1 = 23e-6
 E2=200e9; A2=1200e-6; L2=0.3; a2 = 11.7e-6
 dT = 40;
 ep1=E1*A1/L1; ep2=E2*A2/L2; 
 Ke1=spring1e(ep1); Ke2=spring1e(ep2); 
 fe1 = E1*A1*a1*dT*[-1 1]';
 fe2 = E2*A2*a2*dT*[-1 1]';
%----- Assemble Ke into K ---------------------------------------

 [K,f]=assem(Edof(1,:),K,Ke1,f,fe1)
 [K,f]=assem(Edof(2,:),K,Ke2,f,fe2) 
 
%----- Solve the system of equations ----------------------------

 bc= [1 0; 3 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)
 sigma1=spring1s(ep1,ed1)/A1-E1*a1*dT
 sigma2=spring1s(ep2,ed2)/A2-E2*a2*dT
%---------------------------- end -------------------------------
