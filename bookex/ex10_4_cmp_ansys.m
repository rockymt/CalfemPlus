% example ex10_4
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-06-02
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear
%-----------Global coordinate matrix-----------------------------
Coord=[1     0   0;
	   2   0.4   0;
	   3   0.4   0.15;
	   4   0.4   0.3;
	   5     0   0.3]; 
%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2 3;
       2  5 1 3;
	   3  5 3 4];
%-----------Element type----------------------------------------
ElTp = 'flw2te';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 
%----- Element properties  --------------------------------------
 ep=1;
 D=1.5*[1 0;
        0 1];
%  h{1}=[0 0 0;0 2.5 1.25; 0 1.25 2.5];
%  h{2}=h{1}*0
%  h{3}=[0 0 0;0 2.5 1.25; 0 1.25 2.5];
 h{1}=[0 0 0;0 3.75 0; 0 0 3.75];
 h{2}=h{1}*0;
 h{3}=[0 0 0;0 3.75 0; 0 0 3.75];
 fe{1}=50*25*0.15/2*[0 1 1]';
 fe{2}=fe{1}*0;
 fe{3}=50*25*0.15/2*[0 1 1]';
%----- Assemble Ke into K ---------------------------------------
 K1=K
 for i = 1:3
   Ke=flw2te(Ex(i,:),Ey(i,:),ep,D)
   [K f]=assem(EDof(i,:),K,Ke+h{i},f,fe{i});
 end
 
%----- Solve the system of equations ----------------------------
 bc= [GDof(4,1) 180; GDof(5,1) 180];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:3
	 ed(i,:)=extract(EDof(i,:),a);
	 [es(i,:) et(i,:)]=flw2ts(Ex(i,:),Ey(i,:),D,ed(i,:));
 end
%---------------------------- end -------------------------------
