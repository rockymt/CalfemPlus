function [elefit,fit]=plantefit(El,Coord,phi)

%-----------Element type----------------------------------------
 ElTp = 'flw2te';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------
 [EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end));	 
%----- Stiffness matrix K and load vector f ---------------------
 nDof=max(max(GDof));
 K=sparse(nDof,nDof);
 f=zeros(nDof,1);
 
 %----- Assemble Ke into K ---------------------------------------
 for i = 1:size(EDof,1)
   [Ke fe]=We(Ex(i,:),Ey(i,:),phi(i)); %Element stiffness matrices
   [K f]=assem(EDof(i,:),K,Ke,f,fe);
 end
 fit=K\f;
 for i = 1:size(EDof,1)
    elefit(i,:)=extract(EDof(i,:),fit);
 end
 
 function [Ke fe]=We(ex,ey,phi)
  x23=ex(2)-ex(3); x31=ex(3)-ex(1);
  y23=ey(2)-ey(3); y31=ey(3)-ey(1);
  A=0.5*(-x31*y23+y31*x23);
  Ke=A*1/12.0*[2 1 1;
               1 2 1;
               1 1 2];
  fe=phi*A/3.0*[1 1 1]';
 