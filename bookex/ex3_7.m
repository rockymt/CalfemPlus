 % example ex3_7
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
 Edof=[1   1 3 2;
       2   3 5 4];
%----- Stiffness matrix K and load vector f ---------------------
 K=zeros(5); 
 f=zeros(5,1);
%----- Element stiffness matrices  ------------------------------
 E=1.0e7;  A=0.6;  L=21;
 rho = 0.2836/32.2/12; w=30; 
 r1=10.5; r2= 31.5;
 f1=rho*r1*w^2; f2=rho*r2*w^2;
 ep = E*A/L;
 Ke=quad1De(ep)
 fe1=A*L*f1*[1/6 1/6 2/3]';
 fe2=A*L*f2*[1/6 1/6 2/3]';
%----- Assemble Ke into K ---------------------------------------
 [K f]=assem(Edof(1,:),K,Ke,f,fe1)
 [K f]=assem(Edof(2,:),K,Ke,f,fe2) 
%----- Solve the system of equations ----------------------------
 bc = [1 0];
 [a,r]=solvep(K,f,bc)
 %----- Element forces -------------------------------------------
 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)
 xi = -1:0.2:1;
 sigma1 = quad1Ds(ep,ed1,xi)/A;
 sigma2 = quad1Ds(ep,ed2,xi)/A;
 plot(xi,sigma1,'r'); hold on
 plot(xi,sigma2);
 legend('stress of element 1','stress of element 2')
 xlabel('\xi');ylabel('\sigma (psi)');grid on
 %---------------------------- end -------------------------------
