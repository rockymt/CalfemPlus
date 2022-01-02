function f=plant2T2f(f,Tbd,GDof,Coord,thick)
% [f]=plant2T2f(f,Tbd,GDof,Coord,thick)
%-------------------------------------------------------------
% PURPOSE
%  Create and assemble traction force.
%
% INPUT: f : the input global force vector
%        Tbd :  The row vector Tbd(1:7) contains linear 
%               distributed load information on one edge of an element 
%               Tbd(1): the load type: 
%                1--loads in x and y direction
%                2--pressure(positive if pressure on right side of Tbd(2)-->Tbd(3), 
%                            negtive if tension on right side of Tbd(2)-->Tbd(3)) 
%                   and shear stress(positive if shear stress along Tbd(2)-->Tbd(3),
%                                    negtive if shear stress along Tbd(3)-->Tbd(2)) 
%               Tbd(2:3): the node numbers of the edge
%               Tbd(4:5): the loads of node on node of column 2
%               Tbd(6:7): the loads of node on node of column 3
%               example: Tbd=[2 1  10  -400 0 -400 0];
%                        means uniform tension with value 400  
%                        acting on the edge of nodes 1 and 10.
%          GDof: the global dofs number matrix for nodes in 
%                the sequence of natural numbers
%                [x1 y1 z1 Mx1 My1 Mz1;
%                 x2 y2 z2 Mx2 My2 Mz2;...]
%          Coord : global coordinate matrix
%                  in sequence of natural numbers
%                 [x1 (y1 z1);
%                  x2 (y2 z2);...];
%          thick: thickness
% OUTPUT:  f : the output global force vector with traction force
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2019-03-11
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%------------------------------------------------------------- 

	id1=Tbd(2);
	id2=Tbd(3);

    idCo1=Coord(Coord(:,1)==id1,2:3);
    x_id1=idCo1(1); y_id1=idCo1(2);
    
	idCox=Coord(Coord(:,1)==id2,2:3);
    x_id2=idCox(1); y_id2=idCox(2);
    
	l=sqrt((x_id2-x_id1)^2+(y_id2-y_id1)^2);
    Tx1=Tbd(4); Ty1=Tbd(5);
	Tx2=Tbd(6); Ty2=Tbd(7);
	if Tbd(1)==2  %pressure on right side
        p1=Tx1;p2=Tx2;t1=Ty1;t2=Ty2;
    elseif Tbd(1)==3 %pressure on left side
        p1=-Tx1;p2=-Tx2;t1=Ty1;t2=Ty2;
    end
    if Tbd(1)~=1
        c=(y_id2-y_id1)/l; s=(-x_id2+x_id1)/l;         
        Tx1=-c*p1-s*t1; Tx2=-c*p2-s*t2;
        Ty1=-s*p1+c*t1; Ty2=-s*p2+c*t2;
    end
    
    f(GDof(id1,1:2))=f(GDof(id1,1:2))+thick*l/6.0*[2*Tx1+Tx2;2*Ty1+Ty2];
    f(GDof(id2,1:2))=f(GDof(id2,1:2))+thick*l/6.0*[Tx1+2*Tx2;Ty1+2*Ty2];
    