close all
clear 

%一阶模拟抗混滤波器
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;  % 模拟频率Hz
omiga =  2*pi.*f;         % 模拟角频率 rad/sec
H = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
semilogx( f , 20*log10(H) );
xlabel('模拟频率 (Hz)');
ylabel('dB');            % 截止频率为 0.76/2/pi = 0.121Hz
title('一阶抗混滤波器');


%% 讨论速度的影响
%不同速度下，一阶模拟抗混滤波器幅频特性
%速度为16km/h 36km/h 128km/h
lamda = 1:0.001:10000;        %波长
pesi  = 1./lamda;             %空间频率
v  = 16; v = v/3.6;
H1 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 36;v = v/3.6;
H2 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);
v  = 128;v = v/3.6;
H3 = omiga1 ./ (1j*(2*pi*v.*pesi) + omiga1);

figure1 = figure('Color',[1 1 1]);
semilogx(lamda,20*log10(H1),lamda,20*log10(H2),lamda,20*log10(H3));
% semilogx(pesi,20*log10(H1),pesi,20*log10(H2),pesi,20*log10(H3));
title ('抗混滤波器空间域频响特性')
xlabel('空间频率(1/m)' );
ylabel('幅值(dB)');
legend('16km/h','36km/h','128km/h');
 
%% 这里的方法值得学习一下
%% s = 1j*(2*pi*v.*pesi)
%% z^-1 = exp(-1j*2*pi.*pesi*v*Ti)
%% z = exp(1j*2*pi.*pesi*v*Ti);
%% 
%% 上面是有速度v的表达式
%% 如果没有v，则如下：
%% 参考论文：数字滤波技术在轨道检测中的应用
%% 
% % 一阶数字滤波器
% % 将模拟滤波器转化为数字滤波器 反向差分法
% T = 0.002;
% v = 16;
% H4 = omiga1*T  ./ (1+omiga1*T-exp(-1j.*2*pi*v.*pesi*T));
% v  = 36;
% H5 = omiga1*T  ./ (1+omiga1*T-exp(-1j.*2*pi*v.*pesi*T));
% v  = 128;
% H6 = omiga1*T  ./ (1+omiga1*T-exp(-1j.*2*pi*v.*pesi*T));
% figure
% semilogx(lamda,20*log10(H4),lamda,20*log10(H5),lamda,20*log10(H6));
% title ('数字抗混滤波器');
% xlabel('波长(m)' );
% ylabel('幅值(dB)');
% legend('16km/h','36km/h','128km/h');


%%
%一阶数字补偿滤波器
% 超越方程的系数不利于简化
v    = 16;
v    = v/3.6;
Ti   = 0.25/v;
par1 = omiga1 * Ti;
B1 = ( (1+par1/2)-(1-par1/2).*exp(-1j*2*pi.*pesi*v*Ti) ) ./ par1;

v    = 36;
v    = v/3.6;
Ti   = 0.25/v;
par1 = omiga1 * Ti;
B2 = ( (1+par1/2)-(1-par1/2).*exp(-1j*2*pi.*pesi*v*Ti) ) ./ par1;

v    = 128;
v    = v/3.6;
Ti   = 0.25/v;
par1 = omiga1 * Ti;
B3 = ( (1+par1/2)-(1-par1/2).*exp(-1j*2*pi.*pesi*v*Ti) ) ./ par1;

% figure;
hold on
% semilogx(pesi,20*log10(B1),pesi,20*log10(B2),pesi,20*log10(B3));
semilogx(lamda,20*log10(B1),lamda,20*log10(B2),lamda,20*log10(B3));
title ('补偿滤波器空间域频响特性');
xlabel('空间频率(1/m)' );
ylabel('幅值(dB)');
legend('16km/h','36km/h','128km/h');


%%
%模拟抗混滤波器+数字补偿滤波器
figure1 = figure('Color',[1 1 1]);
% semilogx(pesi,20*log10(B1.*H1),pesi,20*log10(B2.*H2),pesi,20*log10(B3.*H3));
semilogx(lamda,20*log10(B1.*H1),lamda,20*log10(B2.*H2),lamda,20*log10(B3.*H3));
title ('模拟抗混、数字补偿滤波器空间域频响特性');
% xlabel('空间频率(1/m)' );
xlabel('波长(m)' );
ylabel('幅值(dB)');
legend('16km/h','36km/h','128km/h');

%% 注释
%%同样也是一个低通滤波器，把1m以上的波长都留下来？


% %数字抗混滤波器+数字补偿滤波器
% figure
% semilogx(lamda,20*log10(B1.*H4),lamda,20*log10(B2.*H5),lamda,20*log10(B3.*H6));
% title ('数字抗混滤波器+数字补偿滤波器');
% xlabel('波长(m)' );
% ylabel('幅值(dB)');
% legend('16km/h','36km/h','128km/h');




