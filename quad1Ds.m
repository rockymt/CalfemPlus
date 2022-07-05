function [es]=quad1Ds(ep,ed,xi)
% es=spring1s(ep,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute element force quadratic 1D element (quad1De).
%
% INPUT:  ep = [k]        spring stiffness or analog quantity
%         ed = [u1 u2]    element displacements
%                         u1, u2: nodal displacements
%         xi = [xi1 xi2 ...]  local element coordinate                 
%
% OUTPUT: es  = [N]       element force on xi
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-22
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  k = ep;
  u = ed';
  es = (k*[2*xi'-1 2*xi'+1 -4*xi']*u)';
%--------------------------end--------------------------------

