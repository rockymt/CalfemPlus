  function [EDofR,GDof]=DofRaw(El,tp,GDof)
% [EDofR,GDof]=DofRaw(El,type,GDof)
%-------------------------------------------------------------
% PURPOSE
%  Create raw global dof for element of type "tp".
%
% INPUT: El : the element-node connectivity matrix
%             example:elementNumber nodeNumber
%             [1            n1 n2;
%              2            n3 n2;...]
%        tp : the string of element type name in cell or matrix
%              example:'beam2e' 
%        GDof : the inital global dofs vector
%              [x1 y1 z1 Mx1 My1 Mz1 x2 y2 z2 Mx2 ...]
%
% OUTPUT:  EDofR : raw topology matrix
%          GDof : raw global dofs vector
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2022-07-03
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
	case {'beam2e' 'beam2d'}
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
    case 'plani6e'
	  dofn1 = [El(:,2)*6-5 El(:,2)*6-4];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4];
	  dofn5 = [El(:,6)*6-5 El(:,6)*6-4];
	  dofn6 = [El(:,7)*6-5 El(:,7)*6-4];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6])=1;
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
	case 'tetra4e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6-3];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6-3];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4 El(:,4)*6-3];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4 El(:,5)*6-3];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4];
	  GDof([dofn1 dofn2 dofn3 dofn4])=1;
    case 'tetra10e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6-3];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6-3];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4 El(:,4)*6-3];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4 El(:,5)*6-3];
      dofn5 = [El(:,6)*6-5 El(:,6)*6-4 El(:,6)*6-3];
	  dofn6 = [El(:,7)*6-5 El(:,7)*6-4 El(:,7)*6-3];
	  dofn7 = [El(:,8)*6-5 El(:,8)*6-4 El(:,8)*6-3];
	  dofn8 = [El(:,9)*6-5 El(:,9)*6-4 El(:,9)*6-3];
      dofn9 = [El(:,10)*6-5 El(:,10)*6-4 El(:,10)*6-3];
	  dofn10 = [El(:,11)*6-5 El(:,11)*6-4 El(:,11)*6-3];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10])=1;
    case 'soli8e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6-3];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6-3];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4 El(:,4)*6-3];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4 El(:,5)*6-3];
      dofn5 = [El(:,6)*6-5 El(:,6)*6-4 El(:,6)*6-3];
	  dofn6 = [El(:,7)*6-5 El(:,7)*6-4 El(:,7)*6-3];
	  dofn7 = [El(:,8)*6-5 El(:,8)*6-4 El(:,8)*6-3];
	  dofn8 = [El(:,9)*6-5 El(:,9)*6-4 El(:,9)*6-3];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8])=1;
    case 'soli20e'
      dofn1 = [El(:,2)*6-5 El(:,2)*6-4 El(:,2)*6-3];
	  dofn2 = [El(:,3)*6-5 El(:,3)*6-4 El(:,3)*6-3];
	  dofn3 = [El(:,4)*6-5 El(:,4)*6-4 El(:,4)*6-3];
	  dofn4 = [El(:,5)*6-5 El(:,5)*6-4 El(:,5)*6-3];
      dofn5 = [El(:,6)*6-5 El(:,6)*6-4 El(:,6)*6-3];
	  dofn6 = [El(:,7)*6-5 El(:,7)*6-4 El(:,7)*6-3];
	  dofn7 = [El(:,8)*6-5 El(:,8)*6-4 El(:,8)*6-3];
	  dofn8 = [El(:,9)*6-5 El(:,9)*6-4 El(:,9)*6-3];
      dofn9 = [El(:,10)*6-5 El(:,10)*6-4 El(:,10)*6-3];
	  dofn10 = [El(:,11)*6-5 El(:,11)*6-4 El(:,11)*6-3];
	  dofn11 = [El(:,12)*6-5 El(:,12)*6-4 El(:,12)*6-3];
	  dofn12 = [El(:,13)*6-5 El(:,13)*6-4 El(:,13)*6-3];
	  dofn13 = [El(:,14)*6-5 El(:,14)*6-4 El(:,14)*6-3];
	  dofn14 = [El(:,15)*6-5 El(:,15)*6-4 El(:,15)*6-3];
      dofn15 = [El(:,16)*6-5 El(:,16)*6-4 El(:,16)*6-3];
	  dofn16 = [El(:,17)*6-5 El(:,17)*6-4 El(:,17)*6-3];
	  dofn17 = [El(:,18)*6-5 El(:,18)*6-4 El(:,18)*6-3];
	  dofn18 = [El(:,19)*6-5 El(:,19)*6-4 El(:,19)*6-3];
      dofn19 = [El(:,20)*6-5 El(:,20)*6-4 El(:,20)*6-3];
	  dofn20 = [El(:,21)*6-5 El(:,21)*6-4 El(:,21)*6-3];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10 ...
	  dofn11 dofn12 dofn13 dofn14...
	  dofn15 dofn16 dofn17 dofn18 dofn19 dofn20];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10 ...
	  dofn11 dofn12 dofn13 dofn14...
	  dofn15 dofn16 dofn17 dofn18 dofn19 dofn20])=1;
	case 'flw2te'
      dofn1 = [El(:,2)*6-5];
	  dofn2 = [El(:,3)*6-5];
	  dofn3 = [El(:,4)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3];
	  GDof([dofn1 dofn2 dofn3])=1;
    case 'flw2i4e'
	  dofn1 = [El(:,2)*6-5];
	  dofn2 = [El(:,3)*6-5];
	  dofn3 = [El(:,4)*6-5];
	  dofn4 = [El(:,5)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4];
	  GDof([dofn1 dofn2 dofn3 dofn4])=1;
	case 'flw2i6e'
	  dofn1 = [El(:,2)*6-5];
	  dofn2 = [El(:,3)*6-5];
	  dofn3 = [El(:,4)*6-5];
	  dofn4 = [El(:,5)*6-5];
	  dofn5 = [El(:,6)*6-5];
	  dofn6 = [El(:,7)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6])=1;
    case 'flw2i8e'
	  dofn1 = [El(:,2)*6-5];
	  dofn2 = [El(:,3)*6-5];
	  dofn3 = [El(:,4)*6-5];
	  dofn4 = [El(:,5)*6-5];
	  dofn5 = [El(:,6)*6-5];
	  dofn6 = [El(:,7)*6-5];
	  dofn7 = [El(:,8)*6-5];
	  dofn8 = [El(:,9)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8])=1;
	case 'flw3tet10e'
      dofn1 = [El(:,2)*6-5];
	  dofn2 = [El(:,3)*6-5];
	  dofn3 = [El(:,4)*6-5];
	  dofn4 = [El(:,5)*6-5];
	  dofn5 = [El(:,6)*6-5];
	  dofn6 = [El(:,7)*6-5];
	  dofn7 = [El(:,8)*6-5];
	  dofn8 = [El(:,9)*6-5];
      dofn9 = [El(:,10)*6-5];
	  dofn10 = [El(:,11)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10])=1;
	case 'flw3i20e'
      dofn1  = [El(:,2)*6-5];
	  dofn2  = [El(:,3)*6-5];
	  dofn3  = [El(:,4)*6-5];
	  dofn4  = [El(:,5)*6-5];
	  dofn5  = [El(:,6)*6-5];
	  dofn6  = [El(:,7)*6-5];
	  dofn7  = [El(:,8)*6-5];
	  dofn8  = [El(:,9)*6-5];
      dofn9  = [El(:,10)*6-5];
	  dofn10 = [El(:,11)*6-5];
	  dofn11 = [El(:,12)*6-5];
	  dofn12 = [El(:,13)*6-5];
	  dofn13 = [El(:,14)*6-5];
	  dofn14 = [El(:,15)*6-5];
	  dofn15 = [El(:,16)*6-5];
	  dofn16 = [El(:,17)*6-5];
	  dofn17 = [El(:,18)*6-5];
	  dofn18 = [El(:,19)*6-5];
      dofn19 = [El(:,20)*6-5];
	  dofn20 = [El(:,21)*6-5];
	  EDofR = [El(:,1) dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10 ...
	  dofn11 dofn12 dofn13 dofn14...
	  dofn15 dofn16 dofn17 dofn18 dofn19 dofn20];
	  GDof([dofn1 dofn2 dofn3 dofn4...
	  dofn5 dofn6 dofn7 dofn8 dofn9 dofn10 ...
	  dofn11 dofn12 dofn13 dofn14...
	  dofn15 dofn16 dofn17 dofn18 dofn19 dofn20])=1;
	otherwise
	 'Your element type is undefined'
  end  

%--------------------------end--------------------------------
