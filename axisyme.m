  function [Ke,fe]=axisyme(er,ez,D,eq)
% [Ke,fe]=axisyme(er,ez,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a triangular 
%  axis-symmetric element.
%
% INPUT:  er = [r1 r2 r3]         element coordinates
%         ez = [z1 z2 z3]
%         D                       constitutive matrix
%         eq = [br;               br: body force r-dir
%               bz]               bz: body force z-dir
%
% OUTPUT: Ke : element stiffness matrix (6 x 6)
%         fe : equivalent nodal forces (6 x 1)
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-04-23
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%------------------------------------------------------------- 
  br=0.; bz=0.; if nargin==4;  br=eq(1); bz=eq(2); end
  A=1/2*det([ones(3,1) er' ez']);
  if size(D,2)>4
    Dm=D([1 2 3 4],[1 2 3 4]);
  else
    Dm=D;
  end
  z12=ez(1)-ez(2); z23=ez(2)-ez(3); z31=ez(3)-ez(1);
  r21=er(2)-er(1); r32=er(3)-er(2); r13=er(1)-er(3);
  r=sum(er)/3;
  B=1/(2*A)*[z23   0 z31   0 z12   0;
               0 r32   0 r13   0 r21;
			 r32 z23 r13 z31 r21 z12];
  B=[B(1,:);B(2,:);
     1.0/3/r 0 1.0/3/r 0 1.0/3/r 0;
	 B(3,:)];
  Ke=2*pi*r*A*B'*Dm*B;
  fe=2*pi*r*A/3*[br bz br bz br bz]';