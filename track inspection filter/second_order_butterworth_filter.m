close all
clear
clc


%% 是不是就是Fs和Gz两个滤波器？
%% 
% 二阶butterworth抗混滤波器
% 时域幅频特性
a = 10^5 / 2^14;
f = 0.0001 : 0.0001 :10;        % 模拟频率Hz
omiga = 2 * pi .* f;              % 模拟角频率 rad/sec
Fs1   = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);
figure;
semilogx(f,20*log10(Fs1));
xlabel('Ψ(Hz)');
ylabel('dB ');
title('二阶抗混滤波器时域幅频特性曲线');

% 二阶butterworth抗混滤波器
% 空间域幅频特性
k = 10^5 / 2^14;
a = k/2;
b = (k^2-a^2)^(1/2);

%空间频率
pesi = 0.0001 : 0.0001 : 10;     %空间波长 1m ~ 1000m
v   = 16/3.6;
Fs1 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 36/3.6;
Fs2 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
v   = 128/3.6;
Fs3 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));


figure
semilogx(pesi,20*log10(Fs1),pesi,20*log10(Fs2),pesi,20*log10(Fs3));
ylabel('dB ');
xlabel('Ψ(1/m) ');
title('二阶抗混滤波器空间域幅频特性曲线');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            



%% 二阶抗混滤波器对应的补偿滤波器
%% 应该是重新设计了滤波器
% 补偿滤波器
w = 0.0001 : 0.0001 :10;           %空间波长 1m ~ 1000m
v = 16/3.6;
t = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw1 = a1 - a2*exp(-1j*2*pi*0.25.*w) + a3*exp(-2j*2*pi*0.25.*w);

v  = 36/3.6;
t  = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw2 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

v  = 128/3.6;
t  = 0.25/v;
a1 = ( 1 + a*t + a^2*t*t )/( t^2*(a^2+b^2) );
a2 = 2*(1-a^2*t^2) / ( t^2*(a^2+b^2) ) ;
a3 = ( 1 - a*t + a^2*t*t ) / ( t^2*(a^2+b^2) ) ;
Djw3 = a1-a2*exp(-1j*2*pi*0.25.*w)+a3*exp(-2j*2*pi*0.25.*w);

figure
semilogx(w,20*log10(Djw1),w,20*log10(Djw2),w,20*log10(Djw3));
xlabel('Ψ（1/m） ');
ylabel('dB ');
title('二阶补偿滤波器空间域幅频特性曲线');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');            

figure;
semilogx(w,20*log10(Djw1.*Fs1),w,20*log10(Djw2.*Fs2),'b',w,20*log10(Djw3.*Fs3));
xlabel('Ψ（1/m） ');
ylabel('dB ');
title('抗混-补偿滤波器级联幅频特性曲线');
h1=legend('1:v=16km/h','2:v=36km/h','3:v=128km/h');
set(h1,'Box','off');           %不显示方框




%% 这个滤波器表现了抗混叠+互补滤波器的效果，组合而成的低频滤波器可以有效的抗混叠
%%将留下波长为1m以上的，而这部分正是我们想要的
%%这个对比不一定对，所以暂时不去考虑就是了

figure;semilogx(1./w,20*log10(Djw1.*Fs1),1./w,20*log10(Djw2.*Fs2),'b',1./w,20*log10(Djw3.*Fs3));
xlabel('\lamda \m');
%%所以还是低频滤波器




%%
% hold on;
% a = 10^5 / 2^14;
% f = 0.0001 : 0.0001 :10;        % 模拟频率Hz
% omiga = 2 * pi .* f;              % 模拟角频率 rad/sec
% Fs1   = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);
% semilogx(1./f,20*log10(Fs1));
% legend;
% xlabel('WaveLength m')
% ylabel('Mag dB')


% %%
% figure;
% %空间频率
% pesi = 0.0001 : 0.0001 : 10;     %空间波长 1m ~ 1000m
% v   = 3.6/3.6;
% Fs1 = (a^2+b^2) ./(-(2*pi*pesi*v).^2 + 2*1j*a*2*pi.*pesi*v + (a^2+b^2));
% semilogx(1./w,20*log10(Fs1));
% hold on;
% a = 10^5 / 2^14;
% f = 0.0001 : 0.0001 :10;        % 模拟频率Hz
% omiga = 2 * pi .* f;              % 模拟角频率 rad/sec
% Fs1   = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);
% semilogx(1./f,20*log10(Fs1));
% legend;
% xlabel('WaveLength m')
% ylabel('Mag dB')
