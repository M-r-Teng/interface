clear;

global sl lam2
sl=2; lam2=0.005^2;
n=51;L=2000;%/mm 此处的n略微影响求解精度
xmesh = linspace(0,L,n);%离散
solinit = bvpinit(xmesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);

leny=size(sol.y,2); %可以加密
for i=1:leny %也可用向量运算
tao(i)=sqrt(lam2)*exp(-sol.y(1,i))*(1-exp(-sol.y(1,i)));  %tao=Eftf*y'' 设Et=1
end
F=trapz(sol.x(1,:),tao)%有误差,预计偏小

 plot(sol.x(1,:),sol.y(1,:));
 hold on;
 plot(sol.x(1,:),tao);

 show=[sol.x(1,:);fliplr(sol.y(1,:));fliplr(tao)]';

function dydx = bvpfcn(x,y) % 控制方程
global lam2
dydx = zeros(2,1);
dydx = [y(2)
       lam2*exp(-y(1))*(1-exp(-y(1)))];%=y'',本构 若有其他要变的参数，也需global
end
%--------------------------------
function res = bcfcn(ya,yb) % 边界条件
global sl
res = [ya(1)-sl
       yb(2)]; 
end
%--------------------------------
function g = guess(x) %  试函数 y  y'
g = [1/(x+1)
     -1/(x+1)^2];
end
%--------------------------------