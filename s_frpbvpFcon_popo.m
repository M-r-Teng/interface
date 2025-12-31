
clc;clear all;
global lam2 alpha;
Ef=20.39*1000 ;tf=2.54;bf=80;
F=26*1000;%单位N，1Mpa=1N/mm2
alpha=-F/(Ef*bf*tf);
lam2=1/(Ef*tf);
dn=51;L=250;%/mm 此处的dn略微影响求解精度，和系数np区分开
xmesh = linspace(0,L,dn);%离散
solinit = bvpinit(xmesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);

function dydx = bvpfcn(x,y) % 控制方程
global lam2
np=3.01;tmax=4.24;smax=0.101;
dydx = zeros(2,1);
dydx = [y(2)
       lam2*tmax*y(1)/smax*np/(np-1+(y(1)/smax)^np)];%=y'',本构 若有其他要变的参数，也需global
end
%--------------------------------
function res = bcfcn(ya,yb) % 边界条件
global alpha
res = [ya(2)-alpha
       yb(2)]; 
end
%--------------------------------
function g = guess(x) %  试函数 y  y'
g = [1/(x+1)
     -1/(x+1)^2];
end
%--------------------------------
