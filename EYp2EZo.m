  function [Eo]=EYp2EZo(Yp,ex,ey,ez)
% [Eo]=EYp2EZo(Yp,ex,ey,ez)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the global unit vector parallel with
%  the positive local z axis of the beam3e element, 
%  according to the reference point Yp. 

%
% INPUT: Yp=[xp yp zp] 
%             the coordinate of Yp, which is in the
%             first quadrant of local xy plane.
%
%        ex = [x1 x2]        
%        ey = [y1 y2]   
%        ez = [z1 z2]   node coordinates        
%
% OUTPUT:  Eo : the global unit vector
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-04-15
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
 vx =[ ex(2)-ex(1) ey(2)-ey(1) ez(2)-ez(1) ];
 v13=[ Yp(1)-ex(1) Yp(2)-ey(1) Yp(3)-ez(1) ];
 Eo=cross(vx,v13);
 Eo=Eo/sqrt(Eo*Eo');

%--------------------------end--------------------------------