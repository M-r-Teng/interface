%ode45计算F-sl曲线 验证WuJiang
clear;clc;
global lam2 Ef tf alpha beta bf
fco=46.05;Ef=238.1*10^3;tf=0.167*3;bf=50;
alpha=0.103;
kw=1.474;
beta=0.134*sqrt(Ef*tf)/(kw*fco^0.082)
lam2=1/(Ef*tf);%约为2E-5

       global s0 lam2 Ef tf alpha beta bf
   s0=1.03E-06;%输入s0
    L=250;%/mm 
   
[xmesh,sx] = ode45(@(t,y) odefcn(t,y),[0 L],[s0; 0]);

    leny=size(sx,1); 
    sx=sx(:,1);
    for i=1:leny %也可用向量运算
        tao(i)=Ef*tf*alpha/beta^2*exp(-sx(i)/alpha)*(1-exp(-sx(i)/alpha));  
    end
    F=trapz(xmesh,tao)*bf%有误差,预计略偏小 
    sl=sx(leny); 
    
        for i=2:leny 
        strain(i)=(sx(i)-sx(i-1))/(xmesh(i)-xmesh(i-1));  
    end
strain(1)=0;

 for i=1:leny %左右颠倒
    Lx(i)=L-xmesh(i);
 end
 plot(Lx,strain,'-o')
xlabel('x');
ylabel('strain');
strain=fliplr(strain)';Lx=fliplr(Lx)';
strain(1)*1000%输出看端部应变

function dydt = odefcn(t,y)
global lam2 Ef tf alpha beta
  dydt = zeros(2,1);
  dydt(1) = y(2);
  dydt(2) = lam2*Ef*tf*alpha/beta^2*exp(-y(1)/alpha)*(1-exp(-y(1)/alpha));
end
