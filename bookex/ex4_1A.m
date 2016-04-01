% example ex4_1A
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

 EDof=[ 1   1   2   3   4;
        2   5   6   3   4;
        3   1   2   5   6;
        4   7   8   5   6];
	   
%----- Element coordinates --------------------------------------

 Ex=[ 0   40;
     40   40;
      0   40;
      0   40];
 Ey=[ 0    0;
     30    0;
      0   30;
     30   30];      
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(8,8)
 f=zeros(8,1); 
 f(3)=20000; f(6)=-25000
 
%----- Element properties  --------------------------------------
 E=29.5e6; A=1.0;
 ep=[E A]
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:4
   Ke=bar2e(Ex(i,:),Ey(i,:),ep) %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [1 0; 2 0;
      4 0;
	  7 0; 8 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:4
	 ed(i,:)=extract(EDof(i,:),a);
	 N(i)=bar2s(Ex(i,:),Ey(i,:),ep,ed(i,:));
	 sigma(i)=N(i)/A;
 end
%---------------------------- end -------------------------------
