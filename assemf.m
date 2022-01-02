 function [f]=assemf(edof,f,fe)
% [f]=assemf(edof,f,fe)
%-------------------------------------------------------------
% PURPOSE
%  Assembel fe into the global force vector f 
%  according to the topology matrix edof.
%
% INPUT: edof:  topology matrix (n x m, n is 
%                         number of elements) 
%        f   :  the global force vector
%        fe  :  element force vector
%
% OUTPUT:  f :  the new global force vector
%-------------------------------------------------------------

% LAST MODIFIED: Yan LIU  2016-04-23
% Copyright (c)  School of Civil Engineering.
%                Ludong University
%-------------------------------------------------------------
    for i = 1:size(edof,1)
	     t= edof(i,2:end);
         f(t)=f(t)+fe;
    end
%--------------------------end--------------------------------
