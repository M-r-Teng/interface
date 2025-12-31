clear;clc;

Ef=20.39*1000 ;tf=2.54 ;bf=80 ;
tmax=4.24;smax=0.101;np=3.01;

Lb=250;%长度mm 可调
n=51;%相当于论文中的n+1
xmesh = linspace(0,Lb,n);
h=Lb/(n-1);%变量命名基本和Hugo2013论文一致
gamma=tmax/(Ef*tf*smax)*np/(np-1);

for j=1:100
F=0.28*j*1000;%1Mpa=1N/mm2
alpha=-F/(Ef*bf*tf);
lamda=2+gamma*h^2;
Gsk=tmax/(Ef*tf*smax);%Gs前的常数

JA=zeros(n,n);%或diag
sk=zeros(n,1);
count=0;
b=zeros(n,1);
for i=2:n-1
   JA(i,i-1)=1; JA(i,i)=-lamda; JA(i,i+1)=1; 
end
JA(1,1)=-lamda/2;JA(1,2)=1;
JA(n,n-1)=1;JA(n,n)=-lamda/2;

%迭代求解
%for count=1:880
delta=1;%预设相对误差限值
while delta>10^(-6)
for i=1:n
b(i)=Gsk*sk(i)*np/(np-1+(sk(i)/smax)^np)-gamma*sk(i);
end
b(1)=b(1)/2+alpha/h;
b(n)=b(n)/2;
count=count+1;
sk2=JA\b*h^2;%x = A\b 的计算方式与 x = inv(A)*b 不同，论文没说具体程序上怎么做
delta=abs(norm(sk2)-norm(sk))/norm(sk);%模之差
sk=sk2;
end
sl(j)=sk(1) %单个加载端输出
end
