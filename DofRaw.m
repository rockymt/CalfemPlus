  function [EDofR,GDof]=DofRaw(El,tp,GDof)
% [EDofR,GDof]=DofRaw(El,type,GDof)
%-------------------------------------------------------------
% PURPOSE
%  Create raw global dof for element of type "tp".
%
% INPUT: El : the element connectivity matrix
%             example:elementNumber nodeNumber
%             [1            n1 n2;
%              2            n3 n2;...]
%        type : the element type cell or matrix
%              example:'beam2e' 
%        GDof : the inital global dofs vector
%              [x1 y1 z1 Mx1 My1 Mz1 x2 y2 z2 Mx2 ...]
%
% OUTPUT:  EDofR : raw topology matrix
%          GDof : raw global dofs vector
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-04-23
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  
  switch tp
    case 'bar2e'
	  dofn1 = [El(:,2)*6-5 El(:,2)*6-4];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4];
	  EDofR = [El(:,1) dofn1 dofn2];
	  GDof([dofn1 dofn2])=1;
	case 'bar3e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6-3];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6-3];
	  EDofR = [El(:,1) dofn1 dofn2];
	  GDof([dofn1 dofn2])=1;
	case 'beam2e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6];
	  EDofR = [El(:,1) dofn1 dofn2];
	  GDof([dofn1 dofn2])=1;
	case 'beam3e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6-3 ...
          El(:,2)*6-2 El(:,2)*6-1 El(:,2)*6];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6-3 ...
          El(:,3)*6-2 El(:,3)*6-1 El(:,3)*6];
	  EDofR = [El(:,1) dofn1 dofn2];
	  GDof([dofn1 dofn2])=1;
	case {'plante' 'axisyme'}
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3];
	  GDof([dofn1 dofn2 dofn3])=1;
	case 'plani4e'
	  dofn1 = [El(:,2)*6-5 El(:,2)*6-4];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4];
	  GDof([dofn1 dofn2 dofn3 dofn4])=1;
	case 'plani8e'
	  dofn1 = [El(:,2)*6-5 El(:,2)*6-4];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4];
	  dofn5 = [El(:,6)*6-5 El(:,6)*6-4];
	  dofn6 = [El(:,7)*6-5 El(:,7)*6-4];
	  dofn7 = [El(:,8)*6-5 El(:,8)*6-4];
	  dofn8 = [El(:,9)*6-5 El(:,9)*6-4];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8])=1;
	case 'flw2te'
      dofn1 = [El(:,2)*6-5];
	  dofn2 = [El(:,3)*6-5];
	  dofn3 = [El(:,4)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3];
	  GDof([dofn1 dofn2 dofn3])=1;
	otherwise
	 'Your element type is undefined'
  end  

%--------------------------end--------------------------------
