  function [d,Q]=solvep(K,f,bc)
% a=solvep(K,f)
% [a,Q]=solvep(K,f,bc)
%-------------------------------------------------------------
% PURPOSE
%  Solve static FE-equations considering boundary conditions 
%  with penalty method.
%  This function can handle multipoint constraints.
% INPUT: K : global stiffness matrix, dim(K)= nd x nd
%        f : global load vector, dim(f)= nd x 1
%
%        bc : boundary condition matrix
%             If boundary conditions are in form of "u0 = d_p0" 
%             and "u1 = d_p1+u2*d_p2+u3*d_p3+..." etc.
%             where
%             p0, p1, p2, p3...------ the dof index of d
%             u0, u1, u2, u3...------ the coefficients
%            bc should be [p0 u0;
%                          p1 u1;
%                         -p2 u2;
%                         -p3 u3;
%                           ...];
%            negative index means a multipoint constraint follows
% OUTPUT:  a : solution including boundary values
%          Q : reaction force vector
%              dim(a)=dim(Q)= nd x 1, nd : number of dof's
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-02-09
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
  if nargin==2 ; 
     d=K\f ;

  elseif nargin==3;
     C = abs(max(max(K)))*1E4;
     [nd,nd]=size(K);
     fdof=[1:nd]';
%
     d=zeros(size(fdof));
     Q=zeros(size(fdof));
%
     pdof=bc(:,1);
     dp=bc(:,2);
     id = find(pdof>0);
     idd = [id(2:end)-1; length(pdof)];
     pdof = abs(pdof);
     K1 = K; f1 = f;
     for i = 1:length(id)
        dpp = dp(id(i):idd(i));
        dp0 = dpp(1);  dpp(1) = 1;
        pdoff = pdof(id(i):idd(i));
        K1(pdoff,pdoff) = K1(pdoff,pdoff) + C*dpp*dpp';
        f1(pdoff) =  f1(pdoff) + C*dp0*dpp;
     end
%
     d=K1\f1;
%
  end  
     
  Q=K*d-f; 
%--------------------------end--------------------------------
