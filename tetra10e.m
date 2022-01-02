function [Ke,fe]=tetra10e(ex,ey,ez,ep,D,eq)
% Ke=tetra10e(ex,ey,ez,ep,D)
% [Ke,fe]=tetra10e(ex,ey,ez,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a 10 node (tetra) element.
%
% INPUT:   ex = [x1 x2 x3 ... x10]
%          ey = [y1 y2 y3 ... y10]  element coordinates
%          ez = [z1 z2 z3 ... z10]
%
%          ep : = [Ir]               Ir: integration rule
%
%          D                       constitutive matrix
%
%          eq = [bx; by; bz]       bx: body force in x direction
%                                  by: body force in y direction
%                                  bz: body force in z direction
%
% OUTPUT: Ke : element stiffness matrix
%         fe : equivalent nodal forces 
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2020-02-02
% Copyright (c)  School of Civil Engineexing.
%                Ludong Univexsity
%------------------------------------------------------------- 

  ir=ep(1);
  if ir==1  
	ngp=1;
  elseif ir==2
	ngp=4;
  elseif  ir==3
	ngp=5;
  end

  if nargin==5   eq=zeros(3,1);  end
%--------- iron points --------------------------------------
  if ir==1 % 1 pt. quadrature O(h^2)
    w1=1.0/6;
    gp=[0.25 0.25 0.25 0.25];  wp=[ w1 ];
  elseif ir==2 % 4 pt. quadrature O(h^3)
    g1=0.5854101966309658; 
	g2=0.1381966011250105;
	w2=0.25/6;
    gp(:,1)=[g1;g2;g2;g2];  
	gp(:,2)=[g2;g1;g2;g2];
	gp(:,3)=[g2;g2;g1;g2];
	gp(:,4)=[g2;g2;g2;g1];
    wp=[ w2; w2; w2; w2];   
	
  elseif ir==3 % 5 pt. quadrature O(h^4)
    g1=0.5; g2=1.0/6;
    w1=-0.8/6; w2=0.45/6;
    gp(:,1)=[0.25;g1;g2;g2;g2];  
	gp(:,2)=[0.25;g2;g1;g2;g2];
	gp(:,3)=[0.25;g2;g2;g1;g2];
	gp(:,4)=[0.25;g2;g2;g2;g1];
    wp=[ w1; w2; w2; w2; w2]; 
  else
    disp('Used number of integration points not implemented');
    return
  end

  L1=gp(:,1);  L2=gp(:,2);  L3=gp(:,3);   L4=gp(:,4); r2=ngp*3;
%--------- shape functions -----------------------------------
  N(:,1)=2*(L1-0.5).*L1; N(:,5)=4*L1.*L2;
  N(:,2)=2*(L2-0.5).*L2; N(:,6)=4*L1.*L3;
  N(:,3)=2*(L3-0.5).*L3; N(:,7)=4*L1.*L4;
  N(:,4)=2*(L4-0.5).*L4; N(:,8)=4*L2.*L3;
  N(:,9)=4*L3.*L4;
  N(:,10)=4*L2.*L4;

  dNr(1:3:r2,1)=4*L1-1;    dNr(1:3:r2,2)=0;
  dNr(1:3:r2,3)=0;    		dNr(1:3:r2,4)=-4*L4+1;
  dNr(1:3:r2,5)=4.*L2;    	dNr(1:3:r2,6)=4*L3;
  dNr(1:3:r2,7)=4.*(L4-L1); dNr(1:3:r2,8)=0;
  dNr(1:3:r2,9)=-4.*L3;    	dNr(1:3:r2,10)=-4*L2;
  
  dNr(2:3:r2+1,1)=0;  		dNr(2:3:r2+1,2)=4*L2-1;
  dNr(2:3:r2+1,3)=0;  		dNr(2:3:r2+1,4)=-4*L4+1;
  dNr(2:3:r2+1,5)=4*L1;  	dNr(2:3:r2+1,6)=0;
  dNr(2:3:r2+1,7)=-4*L1;  	dNr(2:3:r2+1,8)=4*L3;
  dNr(2:3:r2+1,9)=-4*L3;  	dNr(2:3:r2+1,10)=4*(L4-L2);
  
  dNr(3:3:r2+2,1)=0;  		dNr(3:3:r2+2,2)=0;
  dNr(3:3:r2+2,3)=4*L3-1;  dNr(3:3:r2+2,4)=-4*L4+1;
  dNr(3:3:r2+2,5)=0;  		dNr(3:3:r2+2,6)=4*L1;
  dNr(3:3:r2+2,7)=-4*L1;  	dNr(3:3:r2+2,8)=4*L2;
  dNr(3:3:r2+2,9)=4*(L4-L3); dNr(3:3:r2+2,10)=-4*L2;

%   dN1(1:ngp,1)=4*L1-1;   dN1(1:ngp,2)=0; dN1(1:ngp,3)=0; dN1(1:ngp,4)=0;
%   dN1(1:ngp,5)=4*L2;   dN1(1:ngp,6)=4*L3; dN1(1:ngp,7)=4*L4; 
%   dN1(1:ngp,8)=0;   dN1(1:ngp,9)=0; dN1(1:ngp,10)=0;
%   
%   dN2(1:ngp,1)=0;   dN2(1:ngp,2)=4*L2-1; dN2(1:ngp,3)=0; dN2(1:ngp,4)=0;
%   dN2(1:ngp,5)=4*L1;   dN2(1:ngp,6)=0; dN2(1:ngp,7)=0; 
%   dN2(1:ngp,8)=4*L3;   dN2(1:ngp,9)=0; dN2(1:ngp,10)=4*L4;
%   
%   dN3(1:ngp,1)=0;   dN3(1:ngp,2)=0; dN3(1:ngp,3)=4*L3-1; dN3(1:ngp,4)=0;
%   dN3(1:ngp,5)=0;   dN3(1:ngp,6)=4*L1; dN3(1:ngp,7)=0; 
%   dN3(1:ngp,8)=4*L2;   dN3(1:ngp,9)=4*L4; dN3(1:ngp,10)=0;
%   
%   dN4(1:ngp,1)=0;   dN4(1:ngp,2)=0; dN4(1:ngp,3)=0; dN4(1:ngp,4)=4*L4-1;
%   dN4(1:ngp,5)=0;   dN4(1:ngp,6)=0; dN4(1:ngp,7)=4*L1; 
%   dN4(1:ngp,8)=0;   dN4(1:ngp,9)=4*L3; dN4(1:ngp,10)=4*L2;
%   
%   dNrr(1:3:r2,:)=dN1-dN4; dNrr(2:3:r2+1,:)=dN2-dN4; dNrr(3:3:r2+2,:)=dN3-dN4;

  Ke=zeros(30,30);  
  fe=zeros(30,1);
  JT=dNr*[ex;ey;ez]';

%--------- three dimensional case ----------------------------
  for i=1:ngp
    indx=[ 3*i-2; 3*i-1; 3*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
%     JTinv=inv(JT(indx,:));
%     dNx=JTinv*dNr(indx,:);

    dNx=JT(indx,:)\dNr(indx,:);
    
    B(1,1:3:30-2)=dNx(1,:);
    B(2,2:3:30-1)=dNx(2,:);
    B(3,3:3:30)  =dNx(3,:);
    B(4,1:3:30-2)=dNx(2,:);
    B(4,2:3:30-1)=dNx(1,:);
    B(5,1:3:30-2)=dNx(3,:);
    B(5,3:3:30)  =dNx(1,:);
    B(6,2:3:30-1)=dNx(3,:);
    B(6,3:3:30)  =dNx(2,:);

    N3(1,1:3:30-2)=N(i,:);
    N3(2,2:3:30-1)=N(i,:);
    N3(3,3:3:30)  =N(i,:);

    Ke=Ke+B'*D*B*detJ*wp(i);
    fe=fe+N3'*eq*detJ*wp(i);
  end
%--------------------------end--------------------------------
