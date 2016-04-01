  function [NewCoord,NewEl,SnI,NewSn]=ReNumbNod(OldSn,OldCoord,OldEl)
% [NewCoord,NewEl,SnI,NewSn]=ReNumbNod(OldSn,OldCoord,OldEl)
%-------------------------------------------------------------
% PURPOSE
%      Renumbering the nodes with sequence of 1,2,3,4....
% INPUT:
%       OldSn: old Node number
%       OldCoord: old global coordinate matrix
%       OldEl : the old element connectivity cell or matrix
%            example:elementNumber nodeNumber
%                    {[1            n1 n2;
%                      2            n3 n2;...];
%                     [5            n6 n3;
%                      6            n2 n7;...]}
% OUTPUT:  NewCoord : new global coordinate matrix
%              in sequence of natural numbers
%              [x1 (y1 z1);
%               x2 (y2 z2);...];
%         NewEl : the renumbered element connectivity cell or matrix
%         SnI : the connectivity vector between 
%                the old and new node number
%         NewSn :   the new node number vector
%-------------------------------------------------------------
% LAST MODIFIED: Yan LIU  2016-03-29
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
    [NewNSn SnI] = sort(OldSn);
	NewCoord = OldCoord(SnI,:);
    NewSn(NewNSn) = [1:length(OldSn)]';
    NewEl = OldEl;
	
    if strcmpi(class(OldEl),'cell')
	   ElNb = size(OldEl,2);
	   for i =1:ElNb
	     NewEl{i}(:,2:end) = NewSn(OldEl{i}(:,2:end));
	   end
	else
	    NewEl(:,2:end) = NewSn(OldEl(:,2:end));
	end
%--------------------------end--------------------------------	   
	   