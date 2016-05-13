% example ex6_9free
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-23
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% 
%----------------------------------------------------------------

%-----------Global coordinate matrix-----------------------------
 msh=load_gmsh2('ex6_9free.msh',[3,2])
 Coord=[(1:msh.nbNod)',msh.POS(:,1:2)];
 
%-----------Element connectivity matrix--------------------------	   
 El = [(1:msh.nbTriangles)',msh.TRIANGLES(:,1:3)];

%-----------Element type----------------------------------------
 ElTp = 'plante';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------
 [EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end));	   
 
%----- Draw the fe-mesh as a check of the model -----------------
 figure(1)
 hold on;
 eldraw2(Ex,Ey,[1 3 1]);     
 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=sparse(nDof,nDof);
 f=zeros(nDof,1);

%----- Element properties  --------------------------------------
 ptype=1; thick=0.5;
 ep=[ptype thick];
 E=30e6; v=0.25;
 D=hooke(ptype,E,v);
 
%-------- Create and assemble traction force -------------------
 TBd=[2 1  10  -400 0 -400 0;
      2 10 2   -400 0 -400 0];
 f=plant2T2f(f,TBd,El,EDof,Ex,Ey,ep)

%----- Assemble Ke into K ---------------------------------------

 for i = 1:size(El,1)
   Ke=plante(Ex(i,:),Ey(i,:),ep,D); %Element stiffness matrices
   [K]=assem(EDof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [GDof(4,1) 0; GDof(4,2) 0;
      GDof(3,1) 0; GDof(3,2) 0;
	  GDof(14,1) 0; GDof(14,2) 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
 for i = 1:size(El,1)
	 ed(i,:)=extract(EDof(i,:),a);
	 [es(i,:) et(i,:)]=plants(Ex(i,:),Ey(i,:),ep,D,ed(i,:));
 end
 
 esy=plantefit(El,Coord,es(:,2));

 fill(Ex',Ey',esy');
 axis equal;
 colorbar;
 %----- Draw deformed frame ---------------------------------------
 figure(2)
 hold on;
 sfac=scalfact2(Ex(3,:),Ey(3,:),ed(3,:),0.2);
 for i = 1:size(El,1)
     eldraw2(Ex(i,:),Ey(i,:),[2 1 0]);	 
	 eldisp2(Ex(i,:),Ey(i,:),ed(i,:),[1 2 1],sfac);
 end
 %pltscalb2(sfac,[1e-1 60 0]);
 title('displacements')
 
 