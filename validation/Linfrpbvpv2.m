
for i=1:10
F(i)=Linfrpsol(0.3*i);
close all
end
smesh = linspace(0,30,100);
%plot(smesh,F, '-');
%  F=Linfrpsol(19)
function F = Linfrpsol(sduanin)
%参数定义
global Ep tp tcmax B dmax sdmax scm mumax smumax kmu1 kd1 lamda2 linet1 linet2 linet3 sduan;
Ep=33.2;tp=12.7;
tcmax=6.9349;B=12.1398;
dmax=0.2;sdmax=3.94;scm=0.0571;
mumax=2;smumax=0.4;
kmu1=1.58;kd1=1.49;
lamda2=1/(Ep*tp);
sduan=sduanin;

%数值解
xmesh = linspace(0,10,2000);
solinit = bvpinit(xmesh, @guess);
options = bvpset('RelTol',1e-2,'Stats','on') %可改e-6
%sol = bvp4c(@bvpfcn, @bcfcn, solinit);
sol = bvp4c(@bvpfcn, @bcfcn, solinit,options);

figure
subplot(2,2,1);
plot(sol.x, sol.y(1,:), '-o');
xlabel('L/mm','FontSize',12');
ylabel('s/mm','FontSize',12');
subplot(2,2,2);
plot(sol.x, sol.y(2,:), '-');
xlabel('L/mm','FontSize',12');
ylabel('ds/dL','FontSize',12');
%原始t(s)图
subplot(2,2,3);
x = [linspace(0,scm,100) linspace(scm,0.18,100) linspace(0.18,5,500) linspace(5,100,500)];%此处只用到前边一段，需要精度高。后段相对平缓，前密后疏.100mm后本构的小力暂时忽略。
for i=1:1200
if x(i)>scm %区分第二阶段
    kinp=(x(i)+1.22)/(0.0256*x(i)+0.0032);
    mus(i)=kmu1/smumax*x(i)/(kmu1-1+(x(i)/smumax)^kmu1)*mumax;
    sigmaip(i)=kinp*dmax*((x(i)-scm)/(sdmax-scm))*kd1/(kd1-1+((x(i)-scm)/(sdmax-scm))^kd1);
    tf(i)=mus(i)*sigmaip(i);
else
    tf(i)=0;
end
end
%y=4*tcmax*(exp(-B*x)-exp(-B*2*x));
y=4*tcmax*(exp(-B*x)-exp(-B*2*x))+tf;%理想t
%x0=x(1:100);%用原解析式
x1=x(101:200);
y1=y(101:200);
func=@(oria,x1)oria(1)*x1.^11+oria(2)*x1.^10+oria(3)*x1.^9+oria(4)*x1.^8+oria(5)*x1.^7+oria(6)*x1.^6+oria(7)*x1.^5+oria(8)*x1.^4+oria(9)*x1.^3+oria(10)*x1.^2+oria(11)*x1;%定义待拟合函数，a为系数,x为变量，使得零次项=0
oria=[0 0 0 0 0 0 0 0 0 0 0];%定义系数初始点
fita= lsqcurvefit(func,oria,x1,y1);%得到拟合系数
linet1=[fita 0];
linet2=polyfit(x(201:700),y(201:700),15);
linet3=polyfit(x(701:1200),y(701:1200),11);
ypol0=y(1:100);
ypol1=polyval(linet1,x(101:200));
ypol2=polyval(linet2,x(201:700));
ypol3=polyval(linet3,x(701:1200));
ypol=[ypol0 ypol1 ypol2 ypol3];
plot(x,y,'-');
hold on;
plot(x,ypol,'*'); %图c：理想和拟合后τ-s本构
xlabel('s/mm','FontSize',12');
ylabel('τ/MPa','FontSize',12');
%close all;
lenthsol=size(sol.x);%重取长度，也因为网格数可能会自己变，不同于最初预设
for i=1:lenthsol(2) %分段求对应L的剪力t
   if sol.y(1,i) < scm
       tsol(i) = 4*tcmax*(exp(-B*sol.y(1,i))-exp(-B*2*sol.y(1,i)));
   else if sol.y(1,i) < 0.18
           tsol(i) = polyval(linet1,sol.y(1,i));
       else if sol.y(1,i) < 5
       tsol(i) = polyval(linet2,sol.y(1,i));
           else
               tsol(i) = polyval(linet3,sol.y(1,i));
           end
       end
   end
end
subplot(2,2,4);
plot(sol.x,tsol, '-');
xlabel('L/mm','FontSize',12');
ylabel('t/Mpa','FontSize',12');
F=trapz(sol.x,tsol);%差分求合力

%--------------------------------参数
function dydx = bvpfcn(x,y) % 控制方程
dydx = zeros(2,1);
dydx = [y(2)
       %lamda2*((4*tcmax*(exp(-B*y(1))-exp(-B*2*y(1))))+(y(1)>scm)*(kmu1/smumax*y(1)/(kmu1-1+(y(1)/smumax)^kmu1)*mumax)*(kinp*dmax*((y(1)-scm)/(sdmax-scm))*kd1/(kd1-1+((y(1)-scm)/(sdmax-scm))^kd1)))];
        lamda2*((y(1)>scm)*(y(1)<0.18)*polyval(linet1,y(1))+(y(1)>=0.18)*(y(1)<5)*polyval(linet2,y(1))+(y(1)>5)*polyval(linet3,y(1))+(y(1)<scm)*(y(1)>0)*(4*tcmax*(exp(-B*y(1))-exp(-B*2*y(1)))))];
       %y1是y,y2是y’,dy2就是y''=-f(y) %按条件分段，y1>0需要判断，使得本构的负半段都是0，有利于解稳定
end
%--------------------------------
function res = bcfcn(ya,yb) % 边界条件
res = [ya(1)-sduan
       yb(2)];
end
%--------------------------------
function g = guess(x) % y，y'打靶值
g = [ 1/(x+1)
    -1/(x+1)^2];
end

end