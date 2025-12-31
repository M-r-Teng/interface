clear;clc;
global lam2 L
lam2=5.877*10^(-5);
lamda=sqrt(lam2);
L=100;%长度mm 可调
srange=1;
n=100;deltas=srange/n;

for i=1:n %分n步
[sl(i),F(i)]=frp0con(deltas*i); %如n=1000,srange=5,从s0=0算到0.005*1000=5mm
end

plot(sl,F,'-o')
xlabel('sl');
ylabel('F');

function [sl,F]= frp0con(sinput)
global L lam2
dn=200;%划分份数,另一个n
s0=sinput;%/mm
xmesh = linspace(0,L,dn);
[t,y] = ode45(@vdp1,xmesh,[s0; 0]); %求解初值问题
leny=size(y,1); %可以加密
sl=y(leny,1);
n=1.978;%popo
for i=1:leny %也可用向量运算
tao(i)= 1.007*10^(5)*0.169*lam2*4.31597*y(i,1)/0.033*n/(n-1+(y(i,1)/0.033)^n); 
end
F=trapz(t,tao);%有误差,预计略偏小

end

function dydt = vdp1(t,y)
global lam2
n=1.978;
dydt = [y(2);lam2*4.31597*y(1)/0.033*n/(n-1+(y(1)/0.033)^n)*bool(y(1))]; %本构，y''=
end
function y=bool(x) %  分段sumax
sumax=0.348;
if x>sumax
    y=0;
else
    y=1;
end
end
