%ode45计算F-sl曲线 Wu

clear;
global  lam2

lam2=0.005^2;
dotnum=50; sinmax=5;%可变，单位mm。对此的分步大小只会影响F-L曲线的疏密，但不影响每个F的精度。
geo=1.1 %等比数列的比例
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
sl=sl' ;F=F';
function [sl,F]= frps0con(sinput)
       global s0 lam2
    s0=sinput;
    
    L=1000;%/mm 实际步长会变
   
[xmesh,sx] = ode45(@(t,y) odefcn(t,y),[0 L],[s0; 0]);

    leny=size(sx,1); 
    sx=sx(:,1);
    for i=1:leny %也可用向量运算
        tao(i)=lam2*exp(-sx(i))*(1-exp(-sx(i)));  %tao=Eftf*y'' 设Et=1
    end
    F=trapz(xmesh,tao)/0.005;%有误差,预计略偏小 
    sl=sx(leny);
end

function dydt = odefcn(t,y)
global lam2
  dydt = zeros(2,1);
  dydt(1) = y(2);
  dydt(2) = lam2*exp(-y(1))*(1-exp(-y(1)));
end
