function [TEl] = plant2BdToEl(TBd,TEl)
nbd=size(TBd,1);
TEl=zeros(nbd,1);
for i = 1:nbd
   [i1 j1]=find(TEl(:,2:end)==TBd(i,1));
   [i2 j2]=find(TEl(:,2:end)==TBd(i,2));
   TEl(i)=intersect(i1,i2);
end