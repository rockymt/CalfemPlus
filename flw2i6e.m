function [Ke,fe]=flw2i6e(ex,ey,ep,D,eq)
% Ke=flw2i6e(ex,ey,ep,D)
% [Ke,fe]=flw2i8e(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness (conductivity)
%  matrix for 6 node isoparametric field element
%
% INPUT:  ex = [x1 ... x6]    element coordinates
%         ey = [y1 ... y6]
%                             
%         ep = [t ir]          thickness and 
%                              integration rule
%
%         D  = [kxx kxy;
%               kyx kyy]       constitutive matrix
%
%         eq                   heat supply per unit volume
%
% OUTPUT: Ke :  element 'stiffness' matrix (6 x 6)
%         fe :  element load vector (6 x 1)
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2020-09-30
% Copyright (c)  School of Civil Engineexing.
%                Ludong Univexsity
%------------------------------------------------------------- 
  t=ep(1); ir=ep(2);   
  if ir==1  
	ngp=1;
  elseif ir==2
	ngp=3;
  elseif  ir==3
	ngp=4;
  end
  if nargin==4; eq=0 ; end

%--------- iron points --------------------------------------
  if ir==1 % 1 pt. quadrature O(h^2)
    w1=0.5;
	g0=0.3333333333333333;
    gp=[g0 g0 g0];  wp=[ w1 ];
  elseif ir==2 % 3 pt. quadrature O(h^3)
    g1=0.1666666666666667;
	g2=0.6666666666666667;
	w2=0.1666666666666667;
    gp(:,1)=[g1;g1;g2];  
	gp(:,2)=[g1;g2;g1];
	gp(:,3)=[g2;g1;g1];
    wp=[ w2; w2; w2 ];   
	
  elseif ir==3 % 4 pt. quadrature O(h^4)
    g0=0.3333333333333333; g1=0.2; g2=0.6;
    w1=-0.28125; w2=0.2604166666666667;
    gp(:,1)=[g0;g1;g1;g2];  
	gp(:,2)=[g0;g1;g2;g1];
	gp(:,3)=[g0;g2;g1;g1];
    wp=[ w1; w2; w2; w2 ]; 
  else
    disp('Used number of integration points not implemented');
    return
  end


  L1=gp(:,1);  L2=gp(:,2);  L3=gp(:,3);  r2=ngp*2;
%--------- shape functions -----------------------------------
  N(:,1)=2*(L1-0.5).*L1; N(:,4)=4*L1.*L2;
  N(:,2)=2*(L2-0.5).*L2; N(:,5)=4*L2.*L3;
  N(:,3)=2*(L3-0.5).*L3; N(:,6)=4*L1.*L3;

  dNr(1:2:r2,1)=4*L1-1;     dNr(1:2:r2,2)=0;
  dNr(1:2:r2,3)=-4*L3+1;    dNr(1:2:r2,4)=4*L2;
  dNr(1:2:r2,5)=-4*L2;    	dNr(1:2:r2,6)=4*(L3-L1);

  dNr(2:2:r2+1,1)=0;  		dNr(2:2:r2+1,2)=4*L2-1;
  dNr(2:2:r2+1,3)=-4*L3+1;  dNr(2:2:r2+1,4)=4*L1;
  dNr(2:2:r2+1,5)=4*(L3-L2);  	dNr(2:2:r2+1,6)=-4*L1;


  Ke1=zeros(6,6);  fe1=zeros(6,1);
  JT=dNr*[ex;ey]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminanten lika med noll!')
    end
    % JTinv=inv(JT(indx,:));
    % B=JTinv*dNr(indx,:);

    B=JT(indx,:)\dNr(indx,:);
    Ke1=Ke1+B'*D*B*detJ*wp(i);
    fe1=fe1+N(i,:)'*detJ*wp(i);
  end

  Ke=Ke1*t;  fe=fe1*t*eq;
%--------------------------end--------------------------------
