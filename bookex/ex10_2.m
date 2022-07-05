% example ex10_2
%----------------------------------------------------------------
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-05-24
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear
%----- Topology matrix Edof -------------------------------------
 Edof=[1 1 2;
       2 2 3];
 
%----- Stiffness matrix K and load vector f ---------------------
 K=zeros(3,3);
 f=zeros(3,1);
 
%----- Element stiffness matrices  ------------------------------
 k=0.8; l=6.25e-2;
 Tif=30; h=20; Q=4000;
 ep=k/l;
 eq=Q*l/2;
%----- Assemble Ke into K ---------------------------------------
 for i=1:2
    Ke=spring1e(ep);
    [K f]=assem(Edof(i,:),K,Ke,f,eq*[1 1]');
 end
 
%----- Solve the system of equations ----------------------------
 K(3,3)=K(3,3)+h;
 f(3)=f(3)+h*Tif;  
 [a,r]=solveq(K,f)

%----- Element forces -------------------------------------------
 ed=[];
 for i=1:2
    ed(i,:)=extract(Edof(i,:),a);
    qflux(i)=-spring1s(ep,ed(i,:));
 end
 
%---------------------------- end -------------------------------
