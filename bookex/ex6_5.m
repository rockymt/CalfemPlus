% example ex6_5
%----------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2017-04-23
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------

% REFERENCES
% TR Chandrupatla & AD Belegundu
% Introduction to finite elements in engineering 
%----------------------------------------------------------------
clear
%-----------Global coordinate matrix-----------------------------
Coord=[1 0   0 ;
       2 0   0 ;
	   3 0   0 ;
	   4 0   0 ;
       5 50  40;
       6 60  20;
       7 100 20;
       8 85  40;
	   9 70  60];
 Coord(:,2:3)=Coord(:,2:3)*1e-3;
%-----------Element connectivity matrix--------------------------	
 El=[1 6 7 8;
     2 5 8 9
     5 1 3 2
     7 2 5 4];
%-----------Element type----------------------------------------
 ElTp = 'plante';

%-------- Extract topology matrix Edof &-----------------------
%-------- element coordinate matrices -------------------------
 [EDof,GDof,Ex,Ey]=ExtractElInfo(El,ElTp,Coord(:,2:end));
 
%----- Element properties  -------------------------------------- 
 ptype=1; thick=10e-3;
 ep=[ptype thick];
 
%--------------------traction force------------------------ 
 nDof=max(max(GDof));
 f=zeros(nDof,1);
 TBd=[2 7 8 1e6 0 2e6 0;
     2 8 9 2e6 0 3e6 0];

 for ii=1:size(TBd,1)
     f=plant2T2f(f,TBd(ii,:),GDof,Coord,thick);
 end
 disp(f)
 
 %   TBd=[2 8 7 -2e6 0 -1e6 0;
%        2 9 8 -3e6 0 -2e6 0];

% TBd=[2 8 7 0 1e6 0 1e6;
%      2 9 8 0 1e6 0 1e6];