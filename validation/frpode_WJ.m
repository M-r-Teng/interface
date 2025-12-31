%ode45计算F-sl曲线 验证WuJiang
clear;
global lam2 Ef tf alpha beta bf
fco=46.05;Ef=248.3*10^3;tf=0.167;bf=50;
alpha=0.103;
kw=1.474;
beta=0.134*sqrt(Ef*tf)/(kw*fco^0.082)
lam2=1/(Ef*tf);%约为2E-5

dotnum=200; sinmax=1;%可变，单位mm。对此的分步大小只会影响F-L曲线的疏密，但不影响每个F的精度。
geo=1.15; %等比数列的比例
s0in(dotnum)=sinmax;
for i=dotnum-1:-1:1 %分num步
s0in(i)=s0in(i+1)/geo;
end
for i=1:dotnum %分num步
[sl(i),F(i)]=frps0con(s0in(i)) ;
end
plot(sl,F,'-o')
xlabel('sl');
ylabel('F');
sl=sl';F=F';

function [sl,F]= frps0con(sinput)
       global s0 lam2 Ef tf alpha beta bf
    s0=sinput;
    L=250;%/mm 
   
[xmesh,sx] = ode45(@(t,y) odefcn(t,y),[0 L],[s0; 0]);

    leny=size(sx,1); 
    sx=sx(:,1);
    for i=1:leny %也可用向量运算
        tao(i)=Ef*tf*alpha/beta^2*exp(-sx(i)/alpha)*(1-exp(-sx(i)/alpha));  
    end
    F=trapz(xmesh,tao)*bf;%有误差,预计略偏小 
    sl=sx(leny);
end

function dydt = odefcn(t,y)
global lam2 Ef tf alpha beta
  dydt = zeros(2,1);
  dydt(1) = y(2);
  dydt(2) = lam2*Ef*tf*alpha/beta^2*exp(-y(1)/alpha)*(1-exp(-y(1)/alpha));
end
