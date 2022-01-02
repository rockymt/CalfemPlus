% example ex3_3
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
 clear all
%----- Topology matrix Edof -------------------------------------

 Edof=[1 1 2;
       2 2 3];
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(3,3) 
 f=zeros(3,1);  f(2)=100
 
%----- Element stiffness matrices  ------------------------------
 E=30e6; A1=5.25; A2=3.75; L=12; rhog=0.2836 
 ep1=E*A1/L;  ep2=E*A2/L;
 Ke1=spring1e(ep1)
 Ke2=spring1e(ep2)
 fe1=rhog*A1*L/2*[1 1]'
 fe2=rhog*A2*L/2*[1 1]'
%----- Assemble Ke into K ---------------------------------------

 [K,f]=assem(Edof(1,:),K,Ke1,f,fe1)
 [K,f]=assem(Edof(2,:),K,Ke2,f,fe2) 
 det(K)
 
%----- Solve the system of equations ----------------------------

 bc= [1 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------

 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)

 sigma1=spring1s(ep1,ed1)/A1
 sigma2=spring1s(ep2,ed2)/A2


%---------------------------- end -------------------------------
