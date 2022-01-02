function [El]=tetra10_Ansys_sort(El_Ansys)
El(:,[1 [4 3 1 2 9 6 7 10 8 5]+1])=El_Ansys(:,[1 (1:10)+1]);