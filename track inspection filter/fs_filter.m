close all;
clear all;
%一阶模拟抗混滤波器
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;  % 模拟频率Hz
omiga =  2*pi.*f;         % 模拟角频率 rad/sec
H = omiga1 ./ (1j.*(omiga) + omiga1);
figure
semilogx( f , 20*log10(H) );
xlabel('模拟频率 (Hz)');
ylabel('dB');            % 截止频率为 0.76/2/pi = 0.121Hz
title('一阶抗混滤波器');


%%
%不同速度下，一阶模拟抗混滤波器幅频特性
%速度为16km/h 36km/h 128km/h
lamda = 1:0.001:10000;        %波长
pesi  = 1./lamda;             %空间频率
v  = 16; 
v = v/3.6;
H1 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 36;
v = v/3.6;
H2 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 128;
v = v/3.6;
H3 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);

figure
semilogx(lamda,20*log10(H1),lamda,20*log10(H2),lamda,20*log10(H3));
semilogx(pesi,20*log10(H1),pesi,20*log10(H2),pesi,20*log10(H3));
title ('抗混滤波器空间域频响特性')
xlabel('空间频率(1/m)' );
ylabel('幅值(dB)');
legend('16km/h','36km/h','128km/h');

%% 互补滤波器




