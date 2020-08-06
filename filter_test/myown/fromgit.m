clear all;
close all;

fs = 44100;  %# sampling frequency
fc = 1000;  %# corner frequency of the lowpass 低通频率

% # coefficients of analog lowpass filter 模拟低通滤波器系数
Qinf = 0.8;
sinf = 2*pi*fc ;
C = 1e-6;
L = 1/(sinf^2*C);
R = sinf*L/Qinf;%% 这个应该是涉及到无源函数的性质
B = [0, 0, 1];
A = [L*C, R*C, 1];
% # cofficients of digital filter
T = 1/fs;
b = [T^2, 2*T^2, T^2];
a = [(4*L*C+2*T*R*C+T^2), (-8*L*C+2*T^2), (4*L*C-2*T*R*C+T^2)];

% # compute frequency responses
[ Hd,Om ] = freqz(b, a, 1024,fs);%#离散的
[H,tmp ] = freqs(B, A);%# 连续的模拟信号
figure;
semilogx(tmp/pi/2 , 20*log10(H),'LineWidth',2);
hold on; semilogx(Om , 20*log10(Hd),'LineWidth',2);
set(gca,'Fontname','Times New Roman','fontsize',16);
xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')