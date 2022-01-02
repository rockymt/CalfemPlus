  function [EDof,GDof,Ex,Ey,Ez]=ExtractElInfo(El,ElTp,Coord)
% [EDof,GDof,Ex,Ey,Ez]=ExtractElInfo(El,ElTp,Coord)
%-------------------------------------------------------------
% PURPOSE
%  Extract the information of elements for calfem:
%  including the element dofs, global dofs and coordinates.
%
% INPUT: El : the element connectivity cell or matrix
%            example:elementNumber nodeNumber
%                    {[1            n1 n2;
%                      2            n3 n2;...];
%                     [5            n6 n3;
%                      6            n2 n7;...]}
%        ElTp : the element type cell or matrix
%              example:{'beam2e' 'bar2e'}
%        Coord : global coordinate matrix
%              in sequence of natural numbers
%              [x1 (y1 z1);
%               x2 (y2 z2);...];
% OUTPUT:  EDof: topology matrix , dim(t)= nie x ned+1
%                         nie= number of identical elements
%                         ned= number of element dof's 
%          GDof: the global dofs number matrix for nodes in 
%                the sequence of natural numbers
%                [x1 y1 z1 Mx1 My1 Mz1;
%                 x2 y2 z2 Mx2 My2 Mz2;...]
%          Ex,Ey,Ez : element coordinate matrices
%          Ex=[x1 x2 ...xnen;    one row for each element
%              ...     ...  ;
%              nel     ...  ]  
%              dim= nel x nen ;   nel:number of elemnts
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
    GDof = zeros(1,6*size(Coord,1));
    if strcmpi(class(El),'cell')
	   ElNb = size(El,2);
	   for i =1:ElNb
		 [EDof{i},GDof]=DofRaw(El{i},ElTp{i},GDof);
	   end
	   GDof =ReNumbGDof(GDof);
	   for i =1:ElNb
	      %Renumber element global dof.
	      EDof{i}(:,2:end) = GDof(EDof{i}(:,2:end));
		  [Ex{i},Ey{i},Ez{i}]=CoordElNd(El{i},Coord);
	   end
	else
		[EDof,GDof]=DofRaw(El,ElTp,GDof);
		GDof =ReNumbGDof(GDof);
		%Renumber element global dof.
		EDof(:,2:end) = GDof(EDof(:,2:end)); 
		[Ex,Ey,Ez]=CoordElNd(El,Coord);
	end
	GDof= (reshape(GDof,[6,size(Coord,1)]))';
%--------------------------end--------------------------------	  