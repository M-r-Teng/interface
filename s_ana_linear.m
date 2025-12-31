clear;
tmax=4 ;smax=0.5 ;
Ef=40*1000 ;tf=1 ;bf=1;
F=54.7377;%1Mpa=1N/mm2 F控制，单位N？
Lb=100;
alpha=sqrt(tmax/(smax*Ef*tf));%/mm-1
n=101;
xmesh = linspace(0,Lb,n);

for i=1:n
    cons=F/(alpha*Ef*tf*bf)/sinh(alpha*Lb);%常数提出来统一算，快很多
    s(n+1-i)=cons*cosh(alpha*xmesh(i)-alpha*Lb);
end
plot(xmesh,s,'-o')
xlabel('x');
ylabel('s');

show = [xmesh',s'];