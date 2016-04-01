  function [GDof]=ReNumbGDof(GDof)
% [GDof]=ReNumbGDof(GDof)
%-------------------------------------------------------------
% PURPOSE
%  ReNumber global dof.
%
% INPUT: 
%        GDof : raw global dofs vector
%              [x1 y1 z1 Mx1 My1 Mz1 x2 y2 z2 Mx2 ...]
%
% OUTPUT:  
%        GDof : renumbered global dofs vector
%              [x1 y1 z1 Mx1 My1 Mz1 x2 y2 z2 Mx2 ...]
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-03-13
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  
  RealDof = find(GDof==1);
  GDof(RealDof) = 1:length(RealDof);

%--------------------------end--------------------------------