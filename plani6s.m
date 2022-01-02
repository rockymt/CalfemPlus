 function [es,et,eci]=plani6s(ex,ey,ep,D,ed)
% [es,et,eci]=plani6s(ex,ey,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a 6 node 
%  isoparametric element in plane strain or plane stress.
%
% INPUT:   ex = [x1 ... x6]    element coordinates
%          ey = [y1 ... y6]
%              
%          ep = [ptype t Ir ]  ptype: analysis type
%                              t : thickness 
%                              Ir: integration rule
%                             
%          D                   constitutive matrix
%
%          ed = [u1 u2 ..u12;   element displacement vector
%               ..........]    one row for each element
%
% OUTPUT: es = [ sigx sigy [sigz] tauxy    element stress matrix
%                  ......              ]   one row for each element
%
%         et = [ epsx epsy [epsz] gamxy    element strain matrix
%                  ......              ]   one row for each element
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

%--------- iron points --------------------------------------
  if ir==1 % 1 pt. quadrature O(h)
    w1=0.5;
	g0=0.3333333333333333;
    gp=[g0 g0 g0];  wp=[ w1 ];
  elseif ir==2 % 3 pt. quadrature O(h^2)
    g1=0.1666666666666667;
	g2=0.6666666666666667;
	w2=0.1666666666666667;
    gp(:,1)=[g1;g1;g2];  
	gp(:,2)=[g1;g2;g1];
	gp(:,3)=[g2;g1;g1];
    wp=[ w2; w2; w2 ];   
	
  elseif ir==3 % 4 pt. quadrature O(h^3)
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
  
  rowed=size(ed,1);
  rowex=size(ex,1);
  colD =size(D ,2);
  
  if colD>3
    Cm=inv(D);
    Dm=inv(Cm([1 2 4],[1 2 4]));
  else
    Dm=D;
  end
  
  if rowex==1 incie=0; else incie=1; end
  
  es=[]; et=[]; eci=[]; ie=1;
  for ied=1:rowed
    eci=[eci N*[ex(ie,:);ey(ie,:)]']; 
    JT=dNr*[ex(ie,:);ey(ie,:)]';

    for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:12-1)=dNx(1,:);
      B(2,2:2:12)  =dNx(2,:);
      B(3,1:2:12-1)=dNx(2,:);
      B(3,2:2:12)  =dNx(1,:);

      ee=B*ed(ied,:)';
      if colD>3
         ss=zeros(colD,1);
         ss([1 2 4])=Dm*ee;
         ee=Cm*ss;
      else
         ss=Dm*ee;
      end
      
      et=[et; ee'];
      es=[es; ss'];
    end
    
    ie=ie+incie;
  end
%--------- plane strain --------------------------------------
elseif ptype==2
  
  rowed=size(ed,1);
  rowex=size(ex,1);
  colD =size(D ,2);
    
  if rowex==1 incie=0; else incie=1; end
  
  es=[]; et=[]; eci=[]; ie=1; ee=zeros(colD,1);
  for ied=1:rowed
    eci=[eci N*[ex(ie,:);ey(ie,:)]']; 
    JT=dNr*[ex(ie,:);ey(ie,:)]';

    for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:12-1)=dNx(1,:);
      B(2,2:2:12)  =dNx(2,:);
      B(3,1:2:12-1)=dNx(2,:);
      B(3,2:2:12)  =dNx(1,:);

      e=B*ed(ied,:)';
      if colD>3 ee([1 2 4])=e; else ee=e; end
      
      et=[et; ee'];
      es=[es; (D*ee)'];
    end
    
    ie=ie+incie;
  end
  
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
