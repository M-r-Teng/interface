clear;
%Hugo在Dai基础上加有侧向约束的本构
global s0 lam2 A B tres
A=1.3;B=5.07;tmax=3.58;
tres=1;
Ef=20.39*1000 ;tf=2.54 ;bf=40 ;
lam2=1/(Ef*tf);

s0=3; 
n=101;L=300;%/mm 此处的n略微影响求解精度
xmesh = linspace(0,L,n);%离散
solinit = bvpinit(xmesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);

leny=size(sol.y,2); %可以加密
for i=1:leny %也可用向量运算
s=sol.y(1,i);
tao(i)=A*A*B*(exp(-B*s)-exp(-2*B*s))+tres*(1-exp(-B*s));
end
F=bf*trapz(sol.x,tao)%有误差,预计略偏小

% plot(sol.x(1,:),sol.y(1,:));
% hold on;
% plot(sol.x(1,:),tao);

A=[sol.x
 sol.y(1,:)
tao]'

function dydx = bvpfcn(x,y) % 控制方程
global lam2 A B tres
dydx = zeros(2,1);
dydx = [y(2)
    lam2*(A*A*B*(exp(-B*y(1))-exp(-2*B*y(1)))+tres*(1-exp(-B*y(1))))];%=y'',本构 若有其他要变的参数，也需global
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