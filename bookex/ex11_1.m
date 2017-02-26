% example ex11_1
%----------------------------------------------------------------
% PURPOSE 
%    Structural Dynamics, time integration, full system.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------

% ------ Generate the model ------------------------------------------
%-----------Global coordinate matrix-----------------------------
Coord=[1     0   0;
	   2     0 1.5;
	   3     0   3;
	   4     1   3;
	   5     2   3];

%-----------Element connectivity matrix--------------------------	   
El  = [1  1 2;
       2  2 3;
	   3  3 4;
	   4  4 5];
%-----------Element type----------------------------------------
ElTp = 'beam2d';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end));	   
 
%----- Stiffness & Mass matrix K & M and load vector f ---------
 nDof=max(max(GDof));
 K=sparse(nDof,nDof); M=K;
 
% ------ material data ------------------------------------------
E=3e10;                 rho=2500;
Av=0.1030e-2;           Iv=0.0171e-4;             % IPE100
Ah=0.0764e-2;           Ih=0.00801e-4;            % IPE80
epv=[E Av Iv rho*Av];   eph=[E Ah Ih rho*Ah];

% ------ generate element matrices, assemble in global matrices - 
for i=1:2
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),epv);
  K=assem(EDof(i,:),K,k);  M=assem(EDof(i,:),M,m);  
end
for i=3:4
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),eph);
  K=assem(EDof(i,:),K,k);  M=assem(EDof(i,:),M,m);  
end

% ----- Draw a plot of the element mesh --------------------------
clf;     eldraw2(Ex,Ey,[1 2 2],EDof);
grid;    title('2-D Frame Structure') 

%---------time step & total times--------------------------------
dt=0.002;    T=1;
% ------ the load -----------------------------------------------
G=[0 0; 0.15 1; 0.25 0; T 0];   [t,g]=gfunc(G,dt);
f=zeros(15, length(g));         f(4,:)=1000*g;
% ------ boundary condition, initial condition ------------------
bc=[GDof(1,1) 0; GDof(1,2) 0; GDof(1,6) 0; GDof(5,2) 0];
d0=zeros(15,1);                 v0=zeros(15,1);
% ------ output parameters --------------------------------------
ntimes=[0.1:0.1:1];             nhist=[4 11]; 
% ------ time integration parameters ----------------------------
ip=[dt T 0.25 0.5 10 2 ntimes nhist];

% ------ time integration ---------------------------------------
[Dsnap,D,V,A]=step2(K,[],M,d0,v0,ip,f,bc);

% ----- Plot time history for two DOF:s -------------------------

figure(1), plot(t,D(1,:),'-',t,D(2,:),'--')
grid, xlabel('time (sec)'), ylabel('displacement (m)')
title('Displacement(time) at the 4th and 11th degree-of-freedom')
text(0.3,0.017,'solid line = impact point, x-direction')
text(0.3,0.012,'dashed line = center, horizontal beam, y-direction')

% ----- Plot displacements for some time increments -------------

figure(2),clf, axis('equal'), hold on, axis off
sfac=25;  
title('Snapshots (sec), magnification = 25');
for i=1:5;
  Ext=Ex+(i-1)*3;            eldraw2(Ext,Ey,[2 3 0]); 
  Edb=extract(EDof,Dsnap(:,i));
  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
  Time=num2str(ntimes(i));   text(3*(i-1)+.5,1.5,Time);
end;
Eyt=Ey-4; 
for i=6:10;
  Ext=Ex+(i-6)*3;            eldraw2(Ext,Eyt,[2 3 0]); 
  Edb=extract(EDof,Dsnap(:,i));
  eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
  Time=num2str(ntimes(i));   text(3*(i-6)+.5,-2.5,Time);
end

% ----------------------- end -----------------------------------