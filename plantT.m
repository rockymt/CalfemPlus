function [Te] = plantT(ex,ey,ep,el,Tbd)

Te=zeros(6,1);
id1=find(el(2:end)==Tbd(1));
id2=find(el(2:end)==Tbd(2));
Tx1=Tbd(3); Ty1=Tbd(4);
Tx2=Tbd(5); Ty2=Tbd(6);
thick=ep(2);
x_id1=ex(id1); x_id2=ex(id2);
y_id1=ey(id1); y_id2=ey(id2);
l=sqrt((x_id2-x_id1)^2+(y_id2-y_id1)^2);
Te(2*id1-1:2*id1)=thick*l/6.0*[2*Tx1+Tx2;2*Ty1+Ty2];
Te(2*id2-1:2*id2)=thick*l/6.0*[Tx1+2*Tx2;Ty1+2*Ty2];


