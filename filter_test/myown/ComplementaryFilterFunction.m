% =========================================================================
%
%                  互补滤波器的功能测试
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 11月18日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.发现并没有相加为1的结果出现，可能是在结果上有问题，也可能是文档写错了
%        2.因为还要考虑前面的滤波器，比如Fs或者Bz等
%        3. 
%--------------------------------------------------------------------------


clear all;
close all;
v = 100;
T = 0.25/v;
fs = 1./T;
f = 0.0001:0.0001:4*v/2;

%% H3s 低通滤波器
Omega1 = 0.76;
Omega2 = 6.10;
Omega3 = 0.38;

tp1 = Omega1*Omega2^2*Omega3^3;
tp2 = Omega2*Omega3*(Omega1*Omega2 + Omega2*Omega3 + Omega1 * Omega3);
tp3 = Omega2^2*Omega1;
b = [tp2,tp1];
a = tp3*[1,Omega3,Omega3^2];

s = 2*pi*f*1j;
temp1 = s.^2 + Omega3.*s+Omega3^2;
H3s = (tp1 + tp2.*s) ./(tp3*( temp1 ));

figure1 = figure('Color',[1 1 1]);
semilogx(f , 20*log10(abs(H3s)),'LineWidth',1);



%%
z_1 = exp(-1j*2*pi.*f/fs);
H3z = (Omega3*T)^2 + (Omega3*T)*(1+Omega3/Omega1+Omega3/Omega2)*(1-z_1);
H3z = H3z./( 1 + Omega3*T +(Omega3*T)^2  - ( (2+Omega3*T)*z_1 - z_1.^2 ) );
hold on;
semilogx(f , 20*log10(abs(H3z)),'LineWidth',1);
% legend 模拟滤波器 数字滤波器

%%
% figure;plot(f , 180*(angle(H3z))/pi,'LineWidth',1);

%% 高通滤波器
% Omega1 = 10^5/2^17;
% Omega2 = 10^5/2^14;
% Omega3 = 10^5/2^18;

zeta = Omega1+Omega2+Omega3;
nang = Omega2^2 + Omega2*Omega3 + Omega2*Omega1 + Omega1*Omega3 +Omega3^2;
kai = (Omega1 + Omega3)*Omega2^2 + Omega2*Omega1*Omega3 + (Omega1 + Omega2)*Omega3^2;

temp2 = s.^2 + Omega2.*s+Omega2^2;
P2s = s.^5 + zeta*s.^4 + nang*s.^3 + kai*s.^2;
P2s = P2s./( Omega1.*temp1.*temp2 )./s;

%% plot
hold on;
semilogx(f , 20*log10(abs(P2s)),'LineWidth',1);
semilogx(f , 20*log10(abs(-P2s + H3z*2)),'LineWidth',1);
legend 低通滤波器 高通滤波器 并联滤波器
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('f /Hz');
ylabel('Mag /dB');

%%
figure;plot(f , 180*(angle(P2s))/pi,'LineWidth',1);

