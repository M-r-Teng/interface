clear;clc;
n=100;%份数
srange=7;
deltas=srange/n;
global lam2

lam2=25*10^(-6);
lamda=sqrt(lam2);

%画F-s0
% for j=2:2
%     L=500*j;
%     %变化s0，得到F-s0图
%     for i=1:n
%         ns0(i)=1-exp(-i*deltas);
%         buns0(i)=sqrt(1-ns0(i)^2);
%         Fs0(i)=ns0(i)*buns0(i)*sinh(lamda*L*buns0(i))/(1+ns0(i)*cosh(lamda*L*buns0(i)));
%          sL(i)=log((ns0(i)*cosh(lamda*L*buns0(i))+1)/(1-ns0(i)^2));
%     end
% 
%     plot((1:n)*deltas, Fs0, 'b-');
%      plot((1:n)*deltas, sL, 'r-');
%     hold on;
% end

%画F-sL
% for j=2 %变j可以得系列F-s0图,j=1:*
%     L=0.5*j;
%     %变化s0，得到F-s0图
%     for i=1:n
%         ns0(i)=1-exp(-i*deltas);
%         buns0(i)=sqrt(1-ns0(i)^2);
%         sL(i)=log((ns0(i)*cosh(lamda*L*buns0(i))+1)/(1-ns0(i)^2));
%         Fs0(i)=ns0(i)*buns0(i)*sinh(lamda*L*buns0(i))/(1+ns0(i)*cosh(lamda*L*buns0(i)));
%     end
% 
%     plot(sL, Fs0, 'r-');
%     hold on;
% end
% sL=sL';Fs0=Fs0';
% %画sx分布
s0=0.000579616686757801;
ns0=1-exp(-s0);
buns0=sqrt(1-ns0^2);
deltax=10;
for i=1:200
     sl(i)=log((ns0*cosh(lamda*deltax*i*buns0)+1)/(1-ns0^2));
end
 plot(deltax:deltax:2000, sl, 'r*-');

 show=[deltax:deltax:2000;sl]';