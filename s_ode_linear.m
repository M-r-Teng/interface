clear;clc;
       global s0 lam2
    s0=0.05;
    lam2=0.005^2;
    L=100;%/mm 实际步长会变
  
[xmesh,sx] = ode45(@(t,y) odefcn(t,y),[0 L],[s0; 0]);

    leny=size(sx,1);
    sx=sx(:,1);
    for i=1:leny %也可用向量运算
        tao(i)=bull(sx(i)); 
    end
    F=trapz(xmesh,tao)  
 plot(xmesh, sx);
 tao=tao';
 showsx=sx(leny)
 
function dydt = odefcn(t,y)
global lam2
  dydt = zeros(2,1);
  dydt(1) = y(2);
  dydt(2) = lam2*bull(y(1));
end

function y = bull(x)%本构放在这里
smax=0.5;tmax=4;
if x<=smax
    y=tmax/smax*x;
else
    y=0;
end
end
