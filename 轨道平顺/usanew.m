function [q1,q2,q3]=usanew()     %函数名
fs=1000;      %采样频率
l=0.5:1/fs:300;     %波长范围
Av=0.2095;
Aa=0.0762;
Wc=0.8245;
Ws=0.8209;
k=0.25;
Lc=2*pi/Wc;
Ls=2*pi/Ws;        %参数设置
%美国五级线路
S1=100*k*Av./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %高低不平顺
S2=100*k*Aa./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %方向不平顺
S3=100*2*k*Av./(pi*l.^2*Lc.^2.*(1./l.^2+1/Lc.^2).*(1./l.^2+1/Ls.^2));   %水平及轨距不平顺
%美国六级线路
Av2=0.0339;
Ws2=0.4380;
Ls2=2*pi/Ws2;   %参数设置
S4=100*k*Av2./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %高低不平顺
S5=100*k*Av2./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %方向不平顺
S6=100*2*k*Av2./(pi*l.^2*Lc.^2.*(1./l.^2+1/Lc.^2).*(1./l.^2+1/Ls2.^2));   %水平及轨距不平顺
figure;
loglog(l,l.^2.*S4,'r');   %绘制波长与功率谱密度的双对数坐标图
title('美国六级高低不平顺','fontsize',25);   %题名，并设置字号为25
ylabel('功率谱密度/[mm2/(rad/m)]','fontsize',20);   %y轴名称
grid on   %加网格
xlabel('波长/m','fontsize',20);   %x轴名称
figure
loglog(l,l.^2.*S5,'r');
title('美国六级方向不平顺','fontsize',25);
ylabel('功率谱密度/[mm2/(rad/m)]','fontsize',20);
grid on
xlabel('波长/m','fontsize',20);
figure
loglog(l,l.^2.*S6,'r');
title('美国六级水平及轨距不平顺','fontsize',25);
ylabel('功率谱密度/[mm2/(rad/m)]','fontsize',20);
grid on
xlabel('波长/m','fontsize',20);
