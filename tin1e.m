function [Ke fe]=tin1e(ep,eq);
% Ke=tin1e(ep)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix for 1D tin element.
%  Ref: Chap. 10 of Chandrupatla TR, Belegundu AD. 
%       Introduction to Finite Elements in Engineering, 4th edition. 
% INPUT:  ep = [ep1 ep2];   ep1: conduction stiffness; 
%                           ep2: convection stiffness;
%
% OUTPUT: Ke :            stiffness matrix, dim(Ke)= 2 x 2
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-05-24
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
b=0.; if nargin==2;  b=eq; end
k1 = ep(1);   k2 = ep(2);
Ke = k1*[1  -1;
        -1  1];
Ke = Ke +k2*[2   1;
             1   2];
fe = b*[1 1]';
%--------------------------end--------------------------------
