% example ex5_3
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-06
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear
%-----------Global coordinate matrix-----------------------------
Coord=[1     0   0  0;
	   2     0   3  0;
	   3     3   3  0;
	   4     6   3  0;
	   5     9   0  3]; 
Ref6=[6 6 0];
Ref7=[-3 0 0];
%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2;
       2  2 3;
	   3  3 4;
	   4  4 5];
	   
%-----------Element type----------------------------------------
ElTp = 'beam3e';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------
[EDof,GDof,Ex,Ey,Ez]=ExtractElInfo(El,ElTp,Coord(:,2:end))
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 %K=zeros(nDof,nDof);
 K=sparse(nDof,nDof);
 f=zeros(nDof,1);
 f(GDof(3,3))=240e3;
 f(GDof(4,2))=-60e3; f(GDof(4,6))=-180e3;
 eq=[0 -40e3   0   0;
     0     0   0   0;
	 0     0   0   0;
	 0     0   0   0];
 
%----- Element properties  --------------------------------------
 E=200E9; G=80E9; A=0.01; Iy=0.001; Iz=0.001; J=0.002; %nu=0.25
 ep=[E G A Iy Iz J];
 ERefp=[Ref7;
        Ref6;
        Ref6;
        Ref6];
 
%----- Assemble Ke into K ---------------------------------------

 for i = 1:4
   eo(i,:)=EYp2EZo(ERefp(i,:),Ex(i,:),Ey(i,:),Ez(i,:));
   [Ke,fe]=beam3e(Ex(i,:),Ey(i,:),Ez(i,:),eo(i,:),ep,eq(i,:)); 
   [K,f]=assem(EDof(i,:),K,Ke,f,fe);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,1) 0; GDof(1,2) 0; GDof(1,3) 0;
      GDof(1,4) 0; GDof(1,5) 0; GDof(1,6) 0;
	  GDof(5,1) 0; GDof(5,2) 0; GDof(5,3) 0;
      GDof(5,4) 0; GDof(5,5) 0; GDof(5,6) 0];   
 [a,r]=solveq(K,f,bc)

%----- Section forces -------------------------------------------
 for i = 1:4
	Ed(i,:)=extract(EDof(i,:),a);
	es{i}=beam3s(Ex(i,:),Ey(i,:),Ez(i,:),...
	eo(i,:),ep,Ed(i,:),eq(i,:),21); 	
 end
 
%-----Drawing Section forces -------------------------------------------
 plotflag = ['r','g','b','k'];
 xi = linspace(-1,1,21);
 figure(1)
 for i = 1:4 
     plot(xi,es{i}(:,5),plotflag(i));
	 MyMaxe(i) = max(abs(es{i}(:,5)));	 
	 hold on;
 end
 legend('MyE1','MyE2','MyE3','MyE4')
 Mymax = max(MyMaxe)
 
 figure(2)
 for i = 1:4 
     plot(xi,es{i}(:,6),plotflag(i));
	 MzMaxe(i) = max(abs(es{i}(:,6)));	 
	 hold on;
 end
 legend('MzE1','MzE2','MzE3','MzE4')
 Mzmax = max(MzMaxe)