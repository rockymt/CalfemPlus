function [Ke,fe]=plani6e(ex,ey,ep,D,eq)
% Ke=plani6e(ex,ey,ep,D)
% [Ke,fe]=plani6e(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a 6 node isoparametric
%  element in plane strain or plane stress.
%
% INPUT:  ex = [x1 ...   x6]  element coordinates
%         ey = [y1 ...   y6]
%                             
%         ep =[ptype t ir]    ptype: analysis type
%                             ir: integration rule
%                             t : thickness
%
%         D                   constitutive matrix
%
%         eq = [bx; by]       bx: body force in x direction
%                             by: body force in y direction
%
% OUTPUT: Ke : element stiffness matrix (12 x 12)
%         fe : equivalent nodal forces (12 x 1)
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2020-02-04
% Copyright (c)  School of Civil Engineexing.
%                Ludong Univexsity
%------------------------------------------------------------- 
  ptype=ep(1); t=ep(2);  ir=ep(3);  

  if ir==1  
	ngp=1;
  elseif ir==2
	ngp=3;
  elseif  ir==3
	ngp=4;
  end

  if nargin==4   eq=zeros(2,1);  end
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

  Ke=zeros(12,12);  
  fe=zeros(12,1);
  JT=dNr*[ex;ey]';
%--------- plane stress --------------------------------------
if ptype==1

  colD=size(D,2);
  if colD>3
    Cm=inv(D);
    Dm=inv(Cm([1 2 4],[1 2 4]));
  else
    Dm=D;
  end

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
%    JTinv=inv(JT(indx,:));
%    dNx=JTinv*dNr(indx,:);
	dNx=JT(indx,:)\dNr(indx,:);

    B(1,1:2:12-1)=dNx(1,:);
    B(2,2:2:12)  =dNx(2,:);
    B(3,1:2:12-1)=dNx(2,:);
    B(3,2:2:12)  =dNx(1,:);

    N2(1,1:2:12-1)=N(i,:);
    N2(2,2:2:12)  =N(i,:);

    Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
    fe=fe+N2'*eq*detJ*wp(i)*t;
  end
%--------- plane strain --------------------------------------
elseif ptype==2
  
  colD=size(D,2);
  if colD>3
    Dm=D([1 2 4],[1 2 4]);
  else
    Dm=D;
  end

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
    % JTinv=inv(JT(indx,:));
    % dNx=JTinv*dNr(indx,:);
    dNx=JT(indx,:)\dNr(indx,:);
	
    B(1,1:2:12-1)=dNx(1,:);
    B(2,2:2:12)  =dNx(2,:);
    B(3,1:2:12-1)=dNx(2,:);
    B(3,2:2:12)  =dNx(1,:);

    N2(1,1:2:12-1)=N(i,:);
    N2(2,2:2:12)  =N(i,:);

    Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
    fe=fe+N2'*eq*detJ*wp(i)*t;
  end
  
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------

