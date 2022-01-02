clear
xi=-1:0.1:1;eta=-1:0.1:1;
[Xi,Eta]=meshgrid(xi,eta);
N=(1-Xi).*(1-Eta)/4;
surf(Xi,Eta,N)
contour(Xi,Eta,N)