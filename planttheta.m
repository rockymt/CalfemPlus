function [fe]=planttheta(ex,ey,ep,D,epsilon0)
% [fe]=planttheta(ex,ey,ep,D,epsilon0)
%-----------------------------------------------------------------
% PURPOSE
%  Compute the element load vector induced by initial strain 
%  for "plante" element.
%
% INPUT:  ex = [x1 x2 x3]         element coordinates
%         ey = [y1 y2 y3]
% 
%         ep = [ptype t ]         ptype: analysis type
%                                 t: thickness
% 
%         D                       constitutive matrix
%
%         epsilon0: initial strain dim=3 x 1
%
% OUTPUT: fe   element load vector, dim(Ke)= 6 x 1
%-----------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-04-20
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  
  thick=ep(2);
  x12=ex(1)-ex(2); x23=ex(2)-ex(3); x31=ex(3)-ex(1);
  y12=ey(1)-ey(2); y23=ey(2)-ey(3); y31=ey(3)-ey(1);
  A=0.5*(-x31*y23+y31*x23);
  B=[ y23    0   y31    0  y12    0;
        0 -x23     0 -x31    0 -x12; 
     -x23  y23  -x31  y31 -x12  y12]/2.0/A;
  fe=thick*A*B'*D*epsilon0;
%--------------------------end--------------------------------