clc
clear
x=load('M.txt');
X=load('Dscfl3.txt');
Xf=X(:,1);
Yf=X(:,2);
xw=x(:,1);
yw=x(:,2);
[M,N]=size(x);
zw=zeros(M,1);
Ncx=1;
Nfx=1;
Cx=785.835;
Cy=642.75;
dx=0.0044708;
dy=0.0044708;
sx=1;
[R, T, f, k1] = Tsai(Xf, Yf, xw, yw, zw, Ncx, Nfx, dx, dy, Cx, Cy, sx)
