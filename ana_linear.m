clear;clc;
tmax=4 ;smax=0.5 ;
Ef=40*1000 ;tf=1 ;bf=1;
Lb=50;
alpha=sqrt(tmax/(smax*Ef*tf));
n=100;
xmesh = linspace(0,smax,n);

for i=1:n
    F(i)=alpha*Ef*tf*bf*sinh(alpha*Lb)/cosh(-alpha*Lb)*xmesh(i);%N
end
plot(xmesh,F,'-o')
xlabel('sl');
ylabel('F');
xmesh=xmesh';F=F';