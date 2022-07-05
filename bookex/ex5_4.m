% example ex5_4: bar2d&beam2e
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-06
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% Calfem3.4 /examples/exs7.m 
%----------------------------------------------------------------
clear
%-----------Global coordinate matrix-----------------------------
 Coord=[1 0 0;
 	   2 1 0;
 	   3 0 1;
 	   4 1 1;
 	   5 0 2;
 	   6 1 2];
 
%-----------Element connectivity matrix--------------------------	   
 El1 = [1 1 3;
        2 3 5;
 	   3 2 4;
 	   4 4 6;
 	   5 3 4;
 	   6 5 6];
 El2 = [7 1 4
        8 3 6
 	   9 3 2
 	   10 5 4];
 El = {El1 El2};
 	   
%-----------Element type----------------------------------------
 ElTp = {'beam2e' 'bar2e'};

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------

 [EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end))	   
 
%----- Draw the fe-mesh as a check of the model -----------------
 figure(1)
 hold on;
 eldraw2(Ex{1},Ey{1},[1 3 1]);     
 eldraw2(Ex{2},Ey{2},[1 2 1]);
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=sparse(nDof,nDof);
 f=zeros(nDof,1);
 f(GDof(5,1))=1;
 
%----- Element properties  --------------------------------------
 E1=1; A1=1; I=1;
 E2=1; A2=1;
 ep={[E1 A1 I] [E2 A2]};
 
%----- Assemble Ke into K ---------------------------------------
 for i = 1:size(EDof{1},1)
   Ke=beam2e(Ex{1}(i,:),Ey{1}(i,:),ep{1}); 
   K=assem(EDof{1}(i,:),K,Ke);
 end
 for i = 1:size(EDof{2},1)
   Ke=bar2e(Ex{2}(i,:),Ey{2}(i,:),ep{2}); 
   K=assem(EDof{2}(i,:),K,Ke);
 end
% spy(K)
%----- Solve the system of equations ----------------------------
 bc= [GDof(1,1)  0; GDof(1,2)  0; GDof(1,6) 0;
	  GDof(2,1)  0; GDof(2,2)  0; GDof(2,6) 0];   
 [a,r]=solveq(K,f,bc)

 %----- Draw deformed frame ---------------------------------------
 Ed{1}=extract(EDof{1},a);  
 Ed{2}=extract(EDof{2},a);
 sfac=scalfact2(Ex{1},Ey{1},Ed{1},0.1);
 eldisp2(Ex{1},Ey{1},Ed{1},[2 1 1],sfac);
 eldisp2(Ex{2},Ey{2},Ed{2},[2 1 1],sfac);
 axis([-0.5 1.5 -0.5 2.5]);
 pltscalb2(sfac,[1 60 0]);
 title('displacements')
 
 %----- Section forces -------------------------------------------
  
 for i = 1:size(EDof{1},1)
	es1{i}=beam2s(Ex{1}(i,:),Ey{1}(i,:),ep{1},Ed{1}(i,:)); 
 end
 for i = 1:size(EDof{2},1)
	es2{i}=bar2s(Ex{2}(i,:),Ey{2}(i,:),ep{2},Ed{2}(i,:)); 
 end
%----- Draw normal force diagram --------------------------------
 
 figure(2)
 sfac=scalfact2(Ex{1}(1,:),Ey{1}(1,:),es1{1}(:,1),0.2);
 for i = 1:size(EDof{1},1)
     eldia2(Ex{1}(i,:),Ey{1}(i,:),-es1{i}(:,1),[2 1],sfac);
 end
 for i = 1:size(EDof{2},1)
     eldia2(Ex{2}(i,:),Ey{2}(i,:),-es2{i}*[1; 1],[2 1],sfac);
 end
 axis([-0.5 1.5 -0.5 2.5]);
 pltscalb2(sfac,[0.1 0.5 0]);
 title('normal force')

%----- Draw shear force diagram ---------------------------------
 
 figure(3)
 sfac=scalfact2(Ex{1}(1,:),Ey{1}(1,:),es1{1}(:,2),0.2);
 for i = 1:size(EDof{1},1)
     eldia2(Ex{1}(i,:),Ey{1}(i,:),es1{i}(:,2),[2 1],sfac);
 end
 eldraw2(Ex{2},Ey{2},[1 1 0]);
 axis([-0.5 1.5 -0.5 2.5]);
 pltscalb2(sfac,[0.5 0.5 0]);
 title('shear force') 

%----- Draw moment diagram --------------------------------------
 
 figure(4)
 sfac=scalfact2(Ex{1}(1,:),Ey{1}(1,:),es1{1}(:,3),0.2);
 for i = 1:size(EDof{1},1)
     eldia2(Ex{1}(i,:),Ey{1}(i,:),es1{i}(:,3),[2 1],sfac);
 end
 eldraw2(Ex{2},Ey{2},[1 1 0]);
 axis([-0.5 1.5 -0.5 2.5]);
 pltscalb2(sfac,[0.5 0.5 0]);
 title('moment') 

%------------------------ end -----------------------------------