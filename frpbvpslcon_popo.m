clear;
global  lam2

%lam2=1.13706*10^(-4);
lam2=25*10^(-6);
lamda=sqrt(lam2);
dotnum=100; slmax=4;ddot=slmax/dotnum;%可变，单位mm。对此的分步大小只会影响F-L曲线的疏密，但不影响每个F的精度。
for i=1:dotnum %分num步
F(i)=frplcon(ddot*i); 
end
plot((1:dotnum)*ddot,F,'-o')
xlabel('sl');
ylabel('F');

function F= frplcon(sinput)
global s0 lam2
s0=sinput;
n=100;L=1200;%/mm 此处的n略微影响求解精度
xmesh = linspace(0,L,n);%离散
solinit = bvpinit(xmesh, @guess);
sol = bvp5c(@bvpfcn, @bcfcn, solinit);

leny=size(sol.y,2); %可以加密
for i=1:leny %也可用向量运算
tao(i)=lam2*exp(-sol.y(1,i))*(1-exp(-sol.y(1,i)));  %tao=Eftf*y'' 设Et=1
end
F=trapz(sol.x,tao);%有误差,预计略偏小
end


function dydx = bvpfcn(x,y) % 控制方程
global lam2
dydx = zeros(2,1);
dydx = [y(2)
       lam2*exp(-y(1))*(1-exp(-y(1)))];%=y'',本构 若有其他要变的参数，也需global
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