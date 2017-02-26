% example ex3_5
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
 f=zeros(3,1);  f(2)=60e3
 
%----- Element stiffness matrices  ------------------------------
 E=20e9; A=250e-6; L=0.15
 ep=E*A/L;
 Ke=spring1e(ep)
%----- Assemble Ke into K ---------------------------------------

 [K]=assem(Edof(1,:),K,Ke)
 [K]=assem(Edof(2,:),K,Ke) 
 
%----- Solve the system of equations ----------------------------

 bc= [1 0; 3 1.2e-3];   
 [a,r]=solvep(K,f,bc)

%----- Element forces -------------------------------------------

 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)

 sigma1=spring1s(ep,ed1)/A
 sigma2=spring1s(ep,ed2)/A


%---------------------------- end -------------------------------
