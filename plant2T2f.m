function f=plant2T2f(f,TBd,El,EDof,Ex,Ey,ep)

 TEl = plant2Bd2El(TBd(:,2:end),El);
 for j=1:size(TEl)
    i=TEl(j);
    Te=plantT(Ex(i,:),Ey(i,:),ep,El(i,:),TBd(j,:));
	f=assemf(EDof(i,:),f,Te);
 end
 
function [TEl] = plant2Bd2El(TBd,El)
	nbd=size(TBd,1);
	TEl=zeros(nbd,1);
	for i = 1:nbd
	   [i1 j1]=find(El(:,2:end)==TBd(i,1));
	   [i2 j2]=find(El(:,2:end)==TBd(i,2));
	   TEl(i)=intersect(i1,i2);
	end

function [Te] = plantT(ex,ey,ep,el,Tbd)
	Te=zeros(6,1);
	id1=find(el(2:end)==Tbd(2));
	id2=find(el(2:end)==Tbd(3));
	Tx1=Tbd(4); Ty1=Tbd(5);
	Tx2=Tbd(6); Ty2=Tbd(7);
	thick=ep(2);
	x_id1=ex(id1); x_id2=ex(id2);
	y_id1=ey(id1); y_id2=ey(id2);
	l=sqrt((x_id2-x_id1)^2+(y_id2-y_id1)^2);
	if Tbd(1)==2
        p1=Tx1;p2=Tx2;t1=Ty1;t2=Ty2;
		if (id2-id1==1 || id2-id1==-2)
			c=(y_id2-y_id1)/l; s=(-x_id2+x_id1)/l;         
			Tx1=-c*p1-s*t1; Tx2=-c*p2-s*t2;
			Ty1=-s*p1+c*t1; Ty2=-s*p2+c*t2;
		else
		    c=(y_id1-y_id2)/l; s=(-x_id1+x_id2)/l;
			Tx1=-c*p1+s*t1; Tx2=-c*p2+s*t2;
			Ty1=-s*p1-c*t1; Ty2=-s*p2-c*t2;
		end
	end
	Te(2*id1-1:2*id1)=thick*l/6.0*[2*Tx1+Tx2;2*Ty1+Ty2];
	Te(2*id2-1:2*id2)=thick*l/6.0*[Tx1+2*Tx2;Ty1+2*Ty2];