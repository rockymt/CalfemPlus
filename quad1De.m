function [Ke]=quad1De(ep);
% Ke=spring1e(ep)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix for quadratic 1D element.
%
% INPUT:  ep = [k];       stiffness or analog quantity, 
%                         for bar: ep = EA/L
%
% OUTPUT: Ke :            stiffness matrix, dim(Ke)= 3 x 3
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-22
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
k = ep;  
Ke = k/3*[ 7  1 -8;
           1  7 -8;
		  -8 -8  16];
%--------------------------end--------------------------------
