% example ex7_4
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-20
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear
%-----------Global coordinate matrix-----------------------------
Coord=[1     3    0;
	   2     3  0.5;
	   3   7.5    0;
	   4   7.5  0.5;
	   5    12    0;
	   6    12  0.5]; 

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 3 2;
       2  2 3 4;
	   3  4 3 5;
	   4  4 5 6];
%-----------Element type----------------------------------------
ElTp = 'axisyme';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Er,Ez]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof);
 f=zeros(nDof,1); 
 
%----- Element properties  --------------------------------------
 ptype=3
 E=30e6; v=0.3;
 D=hooke(ptype,E,v);
 rhog=0.283; w=3000*2*pi/60;g=32.2*12;
 
%----- Assemble Ke into K ---------------------------------------
 for i = 1:4
   rb(i,:)=sum(Er(i,:))/3.0
   [Ke,fe]=axisyme(Er(i,:),Ez(i,:),D,[rhog/g*rb(i,:)*w^2; 0]);
   [K,f]=assem(EDof(i,:),K,Ke,f,fe);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(1,2) 0;
      GDof(3,2) 0;
	  GDof(5,2) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:4
	 ed(i,:)=extract(EDof(i,:),a);
	 [es(i,:) et(i,:)]=axisyms(Er(i,:),Ez(i,:),D,ed(i,:));
 end
 rb1=[3; rb];
 es1=[interp1(rb,es(:,3),3,'spline','extrap'); es(:,3)];
 plot(rb1,es1,'-o')
%---------------------------- end -------------------------------
