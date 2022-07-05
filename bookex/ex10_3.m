% example ex10_3
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
 k=360; l=0.1/3.0; t=1e-3;
 Tif=20; h=9;
 ep=[k/l h*l/3.0/t];
 eq=h*Tif*l/t;
%----- Assemble Ke into K ---------------------------------------
 for i=1:3
    [Ke fe]=tin1e(ep,eq);
    [K f]=assem(Edof(i,:),K,Ke,f,fe);
 end

%----- Solve the system of equations ----------------------------
 bc= [1 235];
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 ed=[];
 for i=1:3
    ed(i,:)=extract(Edof(i,:),a);
    qflux(i)=-spring1s(ep(1),ed(i,:));
 end
 H=r(1)*t*1
%---------------------------- end -------------------------------
