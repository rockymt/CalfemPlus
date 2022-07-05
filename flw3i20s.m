function [es,et,eci]=flw3i20s(ex,ey,ez,ep,D,ed)
% [es,et,eci]=flw3i20s(ex,ey,ez,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%   Compute flows or corresponding quantities in the
%   20 node (3-dim) isoparametric field element.
%   Arangement of nodes numbering follows the convention 
%   of Calculix or Ansys
%
%  INPUT:  ex = [x1 x2 x3 ... x20] 
%          ey = [y1 y2 y3 ... y20]       element coordinates
%          ez = [z1 z2 z3 ... z20] 
%  
%          ep = [ir]                    Ir: Integration rule
%
%          D  = [kxx kxy kxz;
%                kyx kyy kyz;
%                kzx kzy kzz]           constitutive matrix
%
%          ed =[ el u1 .... u20]         element  nodal values
%                . ... ..   ..]
%                                       
% OUTPUT: es=[qx qy qz
%             . ...  ..]                 element flow(s)
%
%         et=[gx gy gz
%             . ...  ..]                 element gradient(s)
%  
%         eci=[ix1  iy1  iz1;             location vector
%              ....          ;            nint: number of 
%             ix(nint) iy(nint) iz(nint)]   integration points
%-------------------------------------------------------------
% LAST MODIFIED: Yan Liu   2022-07-03
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
ir=ep(1);  ngp=ir*ir*ir;

%--------- Gauss points --------------------------------------
  if ir==1% 1 pt. quadrature 
    g1=0.0; w1=2.0;
    gp=[ g1 g1 g1];  w=[ w1 w1 w1];
  elseif ir==2 % 8 pt. quadrature 
	g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1;-g1; g1;-g1; g1];  
	gp(:,2)=[-g1;-g1; g1; g1;-g1;-g1; g1; g1]; 
	gp(:,3)=[-g1;-g1;-g1;-g1; g1; g1; g1; g1];
    w(:,1)= w1*ones(8,1);   
	w(:,2)=w(:,1);
	w(:,3)=w(:,1);
  elseif ir==3% 27 pt. quadrature 
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=kron(ones(9,1),[-g1; g2; g1]);
    gp(:,2)=kron(kron(ones(3,1),[-g1; g2; g1]),ones(3,1));
	gp(:,3)=kron([-g1; g2; g1],ones(9,1));
    w(:,1)=kron(ones(9,1),[w1; w2; w1]);
    w(:,2)=kron(kron(ones(3,1),[w1; w2; w1]),ones(3,1));
	w(:,3)=kron([w1; w2; w1],ones(9,1));
  elseif ir==4% 64 pt. quadrature 
    g1=0.861136311594053; g2=0.339981043584856;
    w1=0.339981043584856; w2=0.652145154862546;
    gp(:,1)=kron(ones(16,1),[-g1; -g2; g2; g1]);
    gp(:,2)=kron(kron(ones(4,1),[-g1; -g2; g2; g1]),ones(4,1));
	gp(:,3)=kron([-g1; -g2; g2; g1],ones(16,1));
    w(:,1)=kron(ones(16,1),[w1; w2; w2; w1]);
    w(:,2)=kron(kron(ones(4,1),[w1; w2; w2; w1]),ones(4,1));
	w(:,3)=kron([w1; w2; w2; w1],ones(16,1));
  else
    disp('Used number of integration points not implemented');
    return
  end;

  wp=w(:,1).*w(:,2).*w(:,3);
  xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;

%--------- shape functions -----------------------------------

% N=[N_1 N_2 N_3 N_4 ... N_20    _the first Gauss point
%    N_1 N_2 N_3 N_4 ... N_20    _the second Gauss point
%    ...
%    N_1 N_2 N_3 N_4 ... N_20]   _the last Gauss point

  N(:,1)=(1-xsi).*(1-eta).*(1-zet).*(-xsi-eta-zet-2)/8;  
  N(:,2)=(1+xsi).*(1-eta).*(1-zet).*(+xsi-eta-zet-2)/8;  
  N(:,3)=(1+xsi).*(1+eta).*(1-zet).*(+xsi+eta-zet-2)/8;  
  N(:,4)=(1-xsi).*(1+eta).*(1-zet).*(-xsi+eta-zet-2)/8;  
  N(:,5)=(1-xsi).*(1-eta).*(1+zet).*(-xsi-eta+zet-2)/8;
  N(:,6)=(1+xsi).*(1-eta).*(1+zet).*(+xsi-eta+zet-2)/8;
  N(:,7)=(1+xsi).*(1+eta).*(1+zet).*(+xsi+eta+zet-2)/8;
  N(:,8)=(1-xsi).*(1+eta).*(1+zet).*(-xsi+eta+zet-2)/8;
  
  N(:,9) =(1-xsi.^2).*(1-eta).*(1-zet)/4;
  N(:,11)=(1-xsi.^2).*(1+eta).*(1-zet)/4;
  N(:,13)=(1-xsi.^2).*(1-eta).*(1+zet)/4;
  N(:,15)=(1-xsi.^2).*(1+eta).*(1+zet)/4;
  
  N(:,10)=(1-eta.^2).*(1+xsi).*(1-zet)/4;
  N(:,12)=(1-eta.^2).*(1-xsi).*(1-zet)/4;
  N(:,14)=(1-eta.^2).*(1+xsi).*(1+zet)/4;
  N(:,16)=(1-eta.^2).*(1-xsi).*(1+zet)/4;
  
  N(:,17)=(1-zet.^2).*(1-xsi).*(1-eta)/4;
  N(:,18)=(1-zet.^2).*(1+xsi).*(1-eta)/4;
  N(:,19)=(1-zet.^2).*(1+xsi).*(1+eta)/4;
  N(:,20)=(1-zet.^2).*(1-xsi).*(1+eta)/4;

% dNr=[dN_1/dxsi dN_2/dxsi  ... dN_20/dxsi    _the first Gauss point
%      dN_1/deta dN_2/deta  ... dN_20/deta    _the first Gauss point
%      dN_1/dzet dN_2/dzet  ... dN_20/dzet    _the first Gauss point
%      dN_1/dxsi dN_2/dxsi  ... dN_20/dxsi    _the second Gauss point
%      dN_1/deta dN_2/deta  ... dN_20/deta    _the second Gauss point
%      dN_1/dzet dN_2/dzet  ... dN_20/dzet    _the second Gauss point
%      ...]
%    

  dNr(1:3:r2,1)=-(1-eta).*(1-zet).*(-xsi-eta-zet-1-xsi)/8.;    
  dNr(1:3:r2,2)= (1-eta).*(1-zet).*(+xsi-eta-zet-1+xsi)/8.;
  dNr(1:3:r2,3)= (1+eta).*(1-zet).*(+xsi+eta-zet-1+xsi)/8.;    
  dNr(1:3:r2,4)=-(1+eta).*(1-zet).*(-xsi+eta-zet-1-xsi)/8.;
  dNr(1:3:r2,5)=-(1-eta).*(1+zet).*(-xsi-eta+zet-1-xsi)/8.;    
  dNr(1:3:r2,6)= (1-eta).*(1+zet).*(+xsi-eta+zet-1+xsi)/8.;
  dNr(1:3:r2,7)= (1+eta).*(1+zet).*(+xsi+eta+zet-1+xsi)/8.;    
  dNr(1:3:r2,8)=-(1+eta).*(1+zet).*(-xsi+eta+zet-1-xsi)/8.;
  
  dNr(1:3:r2,9) =-2*xsi.*(1-eta).*(1-zet)/4;
  dNr(1:3:r2,11)=-2*xsi.*(1+eta).*(1-zet)/4;
  dNr(1:3:r2,13)=-2*xsi.*(1-eta).*(1+zet)/4;
  dNr(1:3:r2,15)=-2*xsi.*(1+eta).*(1+zet)/4;
  
  dNr(1:3:r2,10)= (1-eta.^2).*(1-zet)/4;
  dNr(1:3:r2,12)=-(1-eta.^2).*(1-zet)/4;
  dNr(1:3:r2,14)= (1-eta.^2).*(1+zet)/4;
  dNr(1:3:r2,16)=-(1-eta.^2).*(1+zet)/4;
  
  dNr(1:3:r2,17)=-(1-zet.^2).*(1-eta)/4;
  dNr(1:3:r2,18)= (1-zet.^2).*(1-eta)/4;
  dNr(1:3:r2,19)= (1-zet.^2).*(1+eta)/4;
  dNr(1:3:r2,20)=-(1-zet.^2).*(1+eta)/4;
  
  
  dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet).*(-xsi-eta-zet-1-eta)/8.;  
  dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet).*(+xsi-eta-zet-1-eta)/8.;
  dNr(2:3:r2+1,3)= (1+xsi).*(1-zet).*(+xsi+eta-zet-1+eta)/8.;  
  dNr(2:3:r2+1,4)= (1-xsi).*(1-zet).*(-xsi+eta-zet-1+eta)/8.;
  dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet).*(-xsi-eta+zet-1-eta)/8.;  
  dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet).*(+xsi-eta+zet-1-eta)/8.;
  dNr(2:3:r2+1,7)= (1+xsi).*(1+zet).*(+xsi+eta+zet-1+eta)/8.;  
  dNr(2:3:r2+1,8)= (1-xsi).*(1+zet).*(-xsi+eta+zet-1+eta)/8.;
  
  dNr(2:3:r2+1,9) = -(1-xsi.^2).*(1-zet)/4;
  dNr(2:3:r2+1,11)=  (1-xsi.^2).*(1-zet)/4;
  dNr(2:3:r2+1,13)= -(1-xsi.^2).*(1+zet)/4;
  dNr(2:3:r2+1,15)=  (1-xsi.^2).*(1+zet)/4;
  
  dNr(2:3:r2+1,10)= -2*eta.*(1+xsi).*(1-zet)/4;
  dNr(2:3:r2+1,12)= -2*eta.*(1-xsi).*(1-zet)/4;
  dNr(2:3:r2+1,14)= -2*eta.*(1+xsi).*(1+zet)/4;
  dNr(2:3:r2+1,16)= -2*eta.*(1-xsi).*(1+zet)/4;
  
  dNr(2:3:r2+1,17)=-(1-zet.^2).*(1-xsi)/4;
  dNr(2:3:r2+1,18)=-(1-zet.^2).*(1+xsi)/4;
  dNr(2:3:r2+1,19)= (1-zet.^2).*(1+xsi)/4;
  dNr(2:3:r2+1,20)= (1-zet.^2).*(1-xsi)/4;
  
  
  
  dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta).*(-xsi-eta-zet-1-zet)/8.;  
  dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta).*(+xsi-eta-zet-1-zet)/8.;
  dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta).*(+xsi+eta-zet-1-zet)/8.;  
  dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta).*(-xsi+eta-zet-1-zet)/8.;
  dNr(3:3:r2+2,5)= (1-xsi).*(1-eta).*(-xsi-eta+zet-1+zet)/8.;  
  dNr(3:3:r2+2,6)= (1+xsi).*(1-eta).*(+xsi-eta+zet-1+zet)/8.;
  dNr(3:3:r2+2,7)= (1+xsi).*(1+eta).*(+xsi+eta+zet-1+zet)/8.;  
  dNr(3:3:r2+2,8)= (1-xsi).*(1+eta).*(-xsi+eta+zet-1+zet)/8.;
  
  dNr(3:3:r2+2,9) = -(1-xsi.^2).*(1-eta)/4;
  dNr(3:3:r2+2,11)= -(1-xsi.^2).*(1+eta)/4;
  dNr(3:3:r2+2,13)=  (1-xsi.^2).*(1-eta)/4;
  dNr(3:3:r2+2,15)=  (1-xsi.^2).*(1+eta)/4;
  
  dNr(3:3:r2+2,10)=-(1-eta.^2).*(1+xsi)/4;
  dNr(3:3:r2+2,12)=-(1-eta.^2).*(1-xsi)/4;
  dNr(3:3:r2+2,14)= (1-eta.^2).*(1+xsi)/4;
  dNr(3:3:r2+2,16)= (1-eta.^2).*(1-xsi)/4;

  dNr(3:3:r2+2,17)=-2*zet.*(1-xsi).*(1-eta)/4;
  dNr(3:3:r2+2,18)=-2*zet.*(1+xsi).*(1-eta)/4;
  dNr(3:3:r2+2,19)=-2*zet.*(1+xsi).*(1+eta)/4;
  dNr(3:3:r2+2,20)=-2*zet.*(1-xsi).*(1+eta)/4;

  

  eci=N*[ex;ey;ez]';  [red,ced]=size(ed);
  JT=dNr*[ex;ey;ez]';

  for i=1:ngp
    indx=[ 3*i-2; 3*i-1; 3*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
    % JTinv=inv(JT(indx,:));
    % B=JTinv*dNr(indx,:);
	B=JT(indx,:)\dNr(indx,:);
    p1=-D*B*ed';
    p2=B*ed';
    es(i:ngp:ngp*red,:)=p1';
    et(i:ngp:ngp*red,:)=p2';
  end
%--------------------------end--------------------------------
