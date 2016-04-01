function [fe]=bar2theta(ex,ey,ep,epsilon0)
% [fe]=bar2theta(ex,ey,ep,epsilon0)
%-----------------------------------------------------------------
% PURPOSE
%  Compute the element load vector induced by initial strain 
%  for two dimensional bar element.
%
% INPUT:  ex = [x1 x2];
%         ey = [y1 y2];      element node coordinates
%
%         ep = [E A]         E: Young's modulus
%                            A: Cross section area
%         epsilon0: initial strain
%
% OUTPUT: fe   element load vector, dim(Ke)= 4 x 1
%-----------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  E=ep(1);  A=ep(2); 
  fe=E*A*epsilon0*[-1 1]';
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);

  n=b'/L;  G=[   n      zeros(size(n));  
              zeros(size(n))     n   ];

  fe=G'*fe;
%--------------------------end--------------------------------
