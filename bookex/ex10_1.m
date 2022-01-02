% example ex10_1
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
       2 2 3;
	   3 3 4];
 
%----- Stiffness matrix K and load vector f ---------------------
 K=zeros(4,4);
 f=zeros(4,1);
 
%----- Element stiffness matrices  ------------------------------
 k1=20; k2=30; k3=50;
 l1=0.3; l2=0.15; l3=0.15;
 Tif=800; h=25; 
 ep=[];
 ep(1)=k1/l1; ep(2)=k2/l2; ep(3)=k3/l3;
 
%----- Assemble Ke into K --------------------------------------- 
 for i=1:3
    Ke=spring1e(ep(i));
    [K]=assem(Edof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------
 K(1,1)=K(1,1)+h;
 f(1)=f(1)+h*Tif;
 bc= [4 20];   
 [a,r]=solvep(K,f,bc)

%----- Element forces -------------------------------------------
 ed=[];
 for i=1:3
    ed(i,:)=extract(Edof(i,:),a);
    qflux(i)=-spring1s(ep(i),ed(i,:));
 end
 
%---------------------------- end -------------------------------
