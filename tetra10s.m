  function [es,et,eci]=tetra10s_ccx(ex,ey,ez,ep,D,ed)
% [es,et,eci]=tetra10s(ex,ey,ez,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a
%  a 10 node (tetra) isoparametric element.
%  Arangement of nodes numbering follows the convention 
%  of Calculix or Ansys
%
% INPUT:  ex = [x1 x2 x3 ... x10]
%         ey = [y1 y2 y3 ... y10]  element coordinates
%         ez = [z1 z2 z3 ... z10]
%
%         ep = [Ir]               Ir: integration rule
%
%         D                       constitutive matrix
%
%         ed = [u1 u2 ..u30;      element displacement vector
%               ...........]      one row for each element
%  
% OUTPUT: es = [ sigx sigy sigz tauxy tauyz tauxz ;  
%                  ......       ...               ] 
%         element stress matrix, one row for each element  
%         et = [ epsx epsy epsz gamxy gamyz gamxz ; 
%                  ......       ...               ]; 
%         element strain matrix， one row for each element
%         eci=[ix1  iy1  iz1;             location vector
%              ....          ;            nint: number of 
%             ix(nint) iy(nint) iz(nint)]   integration points 
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2020-02-04
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  ir=ep(1);    
  if ir==1  
	ngp=1;
  elseif ir==2
	ngp=4;
  elseif  ir==3
	ngp=5;
  end

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
  N(:,2)=2*(L2-0.5).*L2; N(:,6)=4*L2.*L3;
  N(:,3)=2*(L3-0.5).*L3; N(:,7)=4*L1.*L3;
  N(:,4)=2*(L4-0.5).*L4; N(:,8)=4*L1.*L4;
  N(:,9)=4*L2.*L4;
  N(:,10)=4*L3.*L4;

  dNr(1:3:r2,1)=4*L1-1;     dNr(1:3:r2,2)=0;
  dNr(1:3:r2,3)=-4*L3+1;    dNr(1:3:r2,4)=0;
  dNr(1:3:r2,5)=4*L2;    	dNr(1:3:r2,6)=-4*L2;
  dNr(1:3:r2,7)=4*(L3-L1);  dNr(1:3:r2,8)=4*L4;
  dNr(1:3:r2,9)=0;    	    dNr(1:3:r2,10)=-4*L4;
  
  dNr(2:3:r2+1,1)=0;  		dNr(2:3:r2+1,2)=4*L2-1;
  dNr(2:3:r2+1,3)=-4*L3+1;  dNr(2:3:r2+1,4)=0;
  dNr(2:3:r2+1,5)=4*L1;  	dNr(2:3:r2+1,6)=4*(L3-L2);
  dNr(2:3:r2+1,7)=-4*L1;  	dNr(2:3:r2+1,8)=0;
  dNr(2:3:r2+1,9)=4*L4;  	dNr(2:3:r2+1,10)=-4*L4;
  
  dNr(3:3:r2+2,1)=0;  		dNr(3:3:r2+2,2)=0;
  dNr(3:3:r2+2,3)=-4*L3+1;  dNr(3:3:r2+2,4)=4*L4-1;
  dNr(3:3:r2+2,5)=0;  		dNr(3:3:r2+2,6)=-4*L2;
  dNr(3:3:r2+2,7)=-4*L1;  	dNr(3:3:r2+2,8)=4*L1;
  dNr(3:3:r2+2,9)=4*L2;     dNr(3:3:r2+2,10)=4*(L3-L4);

%--------- three dimensional case ----------------------------
  rowed=size(ed,1);
  rowex=size(ex,1);

  if rowex==1 incie=0; else incie=1; end

  es=[]; et=[]; eci=[]; ie=1;
  for ied=1:rowed
    eci=[eci N*[ex(ie,:);ey(ie,:);ez(ie,:)]']; 
    JT=dNr*[ex(ie,:);ey(ie,:);ez(ie,:)]';

    for i=1:ngp
      indx=[ 3*i-2; 3*i-1; 3*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
%       JTinv=inv(JT(indx,:));
%       dNx=JTinv*dNr(indx,:);
      
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

      ee=B*ed(ied,:)';
      et=[et; ee'];
      es=[es; (D*ee)'];
    end
    
    ie=ie+incie;
 end
%--------------------------end--------------------------------
