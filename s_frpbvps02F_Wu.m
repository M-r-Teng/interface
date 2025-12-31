n=1000;L=200;%/mm
global s0 lam2
s0=0.05;%/mm
%lam2=1.13706*10^(-4);
lam2=25*10^(-6);
lamda=sqrt(lam2);
xmesh = linspace(0,L,n);%离散
solinit = bvpinit(xmesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);

n=1-exp(-sol.y(1,n));
bun=sqrt(1-n^2);
Ys=log((n*cosh(lamda*sqrt(1-n^2)*xmesh)+1)/(1-n^2));
F0=n*bun*sinh(lamda*L*bun)/(1+n*cosh(lamda*L*bun))



subplot(2,2,1);
plot(sol.x, sol.y(1,:), '-*');
subplot(2,2,2);
plot(sol.x, sol.y(2,:), 'r-o');


subplot(2,2,3);
plot(L-xmesh, Ys, 'b-');


function dydx = bvpfcn(x,y) % 控制方程
global lam2
dydx = zeros(2,1);
dydx = [y(2)
       lam2*exp(-y(1))*(1-exp(-y(1)))];%=y'',本构
end
%--------------------------------
function res = bcfcn(ya,yb) % 边界条件
global s0
res = [ya(1)-s0
       yb(2)]; %global
end
%--------------------------------
function g = guess(x) %  试函数 y  y'
g = [1/(x+1)
     -1/(x+1)^2];
end
%--------------------------------