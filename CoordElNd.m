  function [Ex,Ey,Ez]=CoordElNd(El,Coord)
% [Ex,Ey,Ez]=CoordElNd(El,Coord)
%-------------------------------------------------------------
% PURPOSE
%  Extract the coordinates of element nodes
%
% INPUT: El : El : the element connectivity matrix
%             example:elementNumber nodeNumber
%             [1            n1 n2;
%              2            n3 n2;...]
%        Coord : global coordinate matrix
%              in sequence of natural numbers
%              [x1 (y1 z1);
%               x2 (y2 z2);...];
%
% OUTPUT:  Ex,Ey,Ez : element coordinate matrices
%         Ex=[x1 x2 ...xnen;    one row for each element
%             ...     ...  ;
%             nel     ...  ]  
%             dim= nel x nen ;   nel:number of elemnts 
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-03-13
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  Ex = []; Ey = []; Ez = [];
  nsd = size(Coord,2);
  for i = 2:size(El,2)
    Ex = [Ex Coord(El(:,i),1)];
	if nsd>1
	  Ey = [Ey Coord(El(:,i),2)];
	end
	if nsd>2
	  Ez = [Ez Coord(El(:,i),3)];
	end
  end

%--------------------------end--------------------------------