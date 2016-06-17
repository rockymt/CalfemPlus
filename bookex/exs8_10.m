function exs8_10
ex=[1 5 3 0.6 2.5 4 1.5 0.8];
ey=[1 1 6 3.0 1 3.5 4.1 2.0];
xi=linspace(-1,1,4);Eta=xi;
[Xi,Eta]=meshgrid(xi,Eta);
[N]=shapefunc8n4e(Xi(:),Eta(:));
X=N*ex'; Y=N*ey';
plot(X,Y,'*');
end
function [N]=shapefunc8n4e(Xi,Eta)
  N(:,1)=-(1-Xi).*(1-Eta).*(1+Xi+Eta)/4; N(:,5)=(1-Xi.*Xi).*(1-Eta)/2;
  N(:,2)=-(1+Xi).*(1-Eta).*(1-Xi+Eta)/4; N(:,6)=(1+Xi).*(1-Eta.*Eta)/2;
  N(:,3)=-(1+Xi).*(1+Eta).*(1-Xi-Eta)/4; N(:,7)=(1-Xi.*Xi).*(1+Eta)/2;
  N(:,4)=-(1-Xi).*(1+Eta).*(1+Xi-Eta)/4; N(:,8)=(1-Xi).*(1-Eta.*Eta)/2;
end

