% example ex3_4
%----------------------------------------------------------------
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
 f=zeros(3,1);  f(2)=200e3
 
%----- Element stiffness matrices  ------------------------------
 E1=70e9; E2=200e9; A1=2400e-6; A2=600e-6; L1=0.3; L2=0.4
 ep1=E1*A1/L1;  ep2=E2*A2/L2;
 Ke1=spring1e(ep1)
 Ke2=spring1e(ep2)
%----- Assemble Ke into K ---------------------------------------

 [K]=assem(Edof(1,:),K,Ke1)
 [K]=assem(Edof(2,:),K,Ke2) 
 
%----- Solve the system of equations ----------------------------

 bc= [1 0; 3 0];   
 [a,r]=solvep(K,f,bc)

%----- Element forces -------------------------------------------

 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)

 sigma1=spring1s(ep1,ed1)/A1
 sigma2=spring1s(ep2,ed2)/A2


%---------------------------- end -------------------------------
