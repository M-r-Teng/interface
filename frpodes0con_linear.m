clear;clc;
       global lam2
    lam2=0.005^2;
   
    for i=1:200 %分步
[sl(i),F(i)]=frp0con(0.005*i); %分步输入sinput
end
    plot(sl,F,'-o')
xlabel('sl');
ylabel('F');
    sl=sl';F=F';
    
  function [sl,F]= frp0con(sinput)
global lam2
  s0=sinput;
 L=50;%/mm 
 
  sol = ode45(@(t,y) odefcn(t,y),[0 L],[s0; 0]);
  xmesh=linspace(0,L,1000);
sx = deval(sol,xmesh);
    leny=size(sx,2);
    sx=sx(1,:);
    for i=1:leny %也可用向量运算
        tao(i)=bull(sx(i));   
    end
    F=trapz(xmesh,tao); %kN
    sl=sx(leny);

  end
 
function dydt = odefcn(t,y)
global lam2
  dydt = zeros(2,1);
  dydt(1) = y(2);
  dydt(2) = lam2*bull(y(1));
end

function y = bull(x)%本构放在这里
smax=0.5;tmax=4;
if x<smax
    y=tmax/smax*x;
else
    y=0;
end
end
