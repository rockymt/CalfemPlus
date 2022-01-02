% example ex11_2
%----------------------------------------------------------------
% PURPOSE 
%    Structural Dynamics, time integration, reduced system.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------
clear
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

%-------- Extract topology matrix EDof &-----------------------
%-------- element coordinate matrices -------------------------	   
[EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end));	   
 
%----- Stiffness & Mass matrix K & M and load vector f ---------
 nDof=max(max(GDof));
 K=zeros(nDof,nDof); M=K;
 
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

% ----- Eigenvalue analysis --------------------------------------
b=[GDof(1,1) GDof(1,2) GDof(1,6) GDof(5,2)]';
[La,Egv]=eigen(K,M,b);
Freq=sqrt(La)/(2*pi);

% ----- plot one eigenmode ---------------------------------------

figure(1),    clf,     grid,     title('The first eigenmode'), 
eldraw2(Ex,Ey,[2 3 1]); 
Edb=extract(EDof,Egv(:,1));      eldisp2(Ex,Ey,Edb,[1 2 2]);
FreqText=num2str(Freq(1));       text(.5,1.75,FreqText);

% ----- plot eight eigenmodes ------------------------------------

figure(2), clf, axis('equal'), hold on, axis off
sfac=0.5;
title('The first eight eigenmodes (Hz)' )
for i=1:4;
  Edb=extract(EDof,Egv(:,i));
  Ext=Ex+(i-1)*3;                eldraw2(Ext,Ey,[2 3 1]); 
  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
  FreqText=num2str(Freq(i));     text(3*(i-1)+.5,1.5,FreqText);
end;
Eyt=Ey-4; 
for i=5:8;
  Edb=extract(EDof,Egv(:,i));
  Ext=Ex+(i-5)*3;                eldraw2(Ext,Eyt,[2 3 1]); 
  eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
  FreqText=num2str(Freq(i));     text(3*(i-5)+.5,-2.5,FreqText);
end

%---------time step & total times. TWO EIGENVECTORS ARE USED------
dt=0.002;      T=1;      nev=2;
% ----- the load ------------------------------------------------
G=[0 0; 0.15 1; 0.25 0; T 0];        [t,g]=gfunc(G,dt);
f=zeros(nDof, length(g));            f(GDof(2,1),:)=1000*g;
fr=Egv(:,1:nev)'*f;
% ----- reduced system matrices ---------------------------------
kr=Egv(:,1:nev)'*K*Egv(:,1:nev);
mr=Egv(:,1:nev)'*M*Egv(:,1:nev);
% ----- initial condition ---------------------------------------
dr0=zeros(nev,1);                    vr0=zeros(nev,1);
% ----- output parameters ---------------------------------------
ntimes=[0.1:0.1:1];    nhistr=[1:1:nev];   nhist=[4 11];
% ----- time integration parameters -----------------------------
ip=[dt T 0.25 0.5 10 nev ntimes nhistr];
% ----- time integration ----------------------------------------
[Dsnapr,Dr,Vr,Ar]=step2(kr,[],mr,dr0,vr0,ip,fr,[]);
% ----- mapping back to original coordinate system --------------
DsnapR=Egv(:,1:nev)*Dsnapr;
DR=Egv(nhist,1:nev)*Dr;
% ----- plot time history for two DOF:s -------------------------
figure(3), plot(t,DR(1,:),'-',t,DR(2,:),'--')
axis([0    1.0000   -0.0100    0.0200])
grid, xlabel('time (sec)'), ylabel('displacement (m)')
title('Displacement(time) at the 4th and 11th degree-of-freedom')
text(0.3,0.017,'solid line = impact point, x-direction')
text(0.3,0.012,'dashed line = center, horizontal beam, y-direction')
text(0.3,-0.007,'TWO EIGENVECTORS ARE USED')

% ----------------------- end -----------------------------------