
clc;clear all;
global s0 lam2
Ef=20.39*1000 ;tf=2.54;bf=80;
lam2=1/(Ef*tf);
s0=0.2;
n=51;L=250;%/mm 此处的n略微影响求解精度
xmesh = linspace(0,L,n);%离散
solinit = bvpinit(xmesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);
%再反求F
np=3.01;tmax=4.24;smax=0.101;
leny=size(sol.y,2); %可以加密
for i=1:leny %也可用向量运算
tao(i)=bf*tmax*sol.y(1,i)/smax*np/(np-1+(sol.y(1,i)/smax)^np);
end
F=trapz(sol.x,tao);%有误差,预计略偏小？

function dydx = bvpfcn(x,y) % 控制方程
global lam2
np=3.01;tmax=4.24;smax=0.101;
dydx = zeros(2,1);
dydx = [y(2)
       lam2*tmax*y(1)/smax*np/(np-1+(y(1)/smax)^np)];%=y'',本构 若有其他要变的参数，也需global
end
%--------------------------------
function res = bcfcn(ya,yb) % 边界条件
global s0
res = [ya(1)-s0
       yb(2)]; 
end
%--------------------------------
function g = guess(x) %  试函数 y  y'
g = [1/(x+1)
     -1/(x+1)^2];
end
%--------------------------------