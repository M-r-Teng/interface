clear;
       global s0 lam2
    s0=1;
    L=2000;%/mm 实际步长会变 
    n=50;
   xmesh = linspace(0,L,n);%离散
[xmesh,sx] = ode45(@(t,y) odefcn(t,y),[0 L],[s0; 0]);

    leny=size(sx,1);
    sx=sx(:,1);
    for i=1:leny %也可用向量运算
        tao(i)=sqrt(lam2)*exp(-sx(i))*(1-exp(-sx(i)));  %tao=Eftf*y'' 设Et=1
    end
    F=trapz(xmesh,tao)%有误差,预计略偏小 
    sl=sx(leny);

 plot(xmesh, sx);

 show=[xmesh,sx,tao'];
 
function dydt = odefcn(t,y)
global lam2
  dydt = zeros(2,1);
  dydt(1) = y(2);
  dydt(2) = lam2*exp(-y(1))*(1-exp(-y(1)));
end