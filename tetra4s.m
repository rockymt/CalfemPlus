  function [es,et]=tetra4s(ex,ey,ez,D,ed)
% [es,et]=tetra4s(ex,ey,ez,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for 4 node (tetra) 
%  element.
%
% INPUT:   ex = [x1 x2 x3 ... x4]
%          ey = [y1 y2 y3 ... y4]  element coordinates
%          ez = [z1 z2 z3 ... z4]
%
%          D                       constitutive matrix
%
%         ed = [u1 u2 ..u12;      element displacement vector
%               ...........]      one row for each element
%  
% OUTPUT: es = [ sigx sigy sigz sigxy sigyz sigxz ;  element stress matrix
%                  ......       ...               ]; one row for each element
%         et = [ epsx epsy epsz gamxy gamyz gamxz ; element strain matrix
%                  ......       ...               ]; one row for each element
%
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-05-24
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

  x14=ex(1)-ex(4); x24=ex(2)-ex(4); x34=ex(3)-ex(4);
  y14=ey(1)-ey(4); y24=ey(2)-ey(4); y34=ey(3)-ey(4);
  z14=ez(1)-ez(4); z24=ez(2)-ez(4); z34=ez(3)-ez(4);

  detJ=x14*(y24*z34-y34*z24)+y14*(z24*x34-z34*x24)+z14*(x24*y34-x34*y24);
  A=[y24*z34-y34*z24 y34*z14-y14*z34 y14*z24-y24*z14;
     z24*x34-z34*x24 z34*x14-z14*x34 z14*x24-z24*x14;
     x24*y34-x34*y24 x34*y14-x14*y34 x14*y24-x24*y14];
  A=A/detJ;
  B11=diag(A(:,1));B12=diag(A(:,2));B13=diag(A(:,3));
  B21=[0 A(3,1) A(2,1);
      A(3,1) 0 A(1,1);
	  A(2,1) A(1,1) 0];
  B22=[0 A(3,2) A(2,2);
      A(3,2) 0 A(1,2);
	  A(2,2) A(1,2) 0]; 
  B23=[0 A(3,3) A(2,3);
      A(3,3) 0 A(1,3);
	  A(2,3) A(1,3) 0];
  B=[B11 B12 B13 -(B11+B12+B13);
     B21 B22 B23 -(B21+B22+B23)];
  e=B*ed';
  et=e';
  es=(D*e)';  
  end