%一阶二阶滤波器时域上的福频特性
% 人为设定时间频率为  t= 0.002
clear all;
close all;
% 一阶模拟
Q1 = 2^17;
Q1 = 100000/Q1;
b = [0 Q1];
a =[1 Q1];
[h0 f0] = freqs(b,a);
figure();suptitle('一阶模拟抗混叠时域福频特性曲线');subplot(2,1,1);
semilogx(f0/(2*pi),20*log10(abs(h0)));
xlabel('频率（Hz）');ylabel('幅值(dB)');
subplot(2,1,2);
semilogx(f0/(2*pi),angle(h0)*180/pi);
xlabel('频率（Hz）');ylabel('相位');
% 一阶时域上
Q1 = 2^17;
Q1 = 100000/Q1;
t = 0.002;
b = [Q1*t 0];
a = [1+Q1*t -1];
[h,f] = freqz (b,a,15000,500);
figure();suptitle('一阶数字抗混叠时域福频特性曲线');subplot(2,1,1);
semilogx(f,20*log10(abs(h)));
xlabel('频率（Hz）');ylabel('幅值(dB)');
subplot(2,1,2);
semilogx(f,angle(h)*180/pi);
xlabel('频率（Hz）');ylabel('相位');
% 补偿
bc = [1+(Q1*t/2)  (-1)+(Q1*t/2)];
ac = [Q1*t 0];
bc = [ 1+(Q1*t/2)  (-1)+(Q1*t/2)];
ac = [Q1*t 0];
[h1,f1] = freqz (bc,ac,15000,500);
figure();suptitle('一阶数字补偿时域福频特性曲线');subplot(2,1,1);
semilogx(f1,20*log10(abs(h1)));
xlabel('频率（Hz）');ylabel('幅值(dB)');
subplot(2,1,2);
semilogx(f1,angle(h1)*180/pi);
xlabel('频率（Hz）');ylabel('相位');

% ----二阶模拟-----
Q = 2^14;
Q = 100000/Q;
b = [0 0 Q^2];
a = [1 Q Q^2];
[h0 f0] = freqs(b,a);
figure();suptitle('二阶模拟抗混叠时域福频特性曲线');subplot(2,1,1);
semilogx(f0/(2*pi),20*log10(abs(h0)));
xlabel('频率（Hz）');ylabel('幅值(dB)');
subplot(2,1,2);
semilogx(f0/(2*pi),angle(h0)*180/pi);
xlabel('频率（Hz）');ylabel('相位');
% ------二阶时域上的-----
Q = 2^14;
Q = 100000/Q;
t = 0.002;
b0 = (Q*t)^2;
a0 = (Q*t)^2+Q*t+1;
a1 = -(Q*t+2);
a2 = 1;
%b = [0 0 b0];
%a = [a2 a1 a0];
b = [b0 0 0];
a = [a0 a1 a2];
[h,f] = freqz (b,a,1500,500);
figure();suptitle('二阶数字抗混叠时域福频特性曲线');subplot(2,1,1);
semilogx(f,20*log10(abs(h)),'r');
xlabel('频率（Hz）');ylabel('幅值(dB)');
subplot(2,1,2);
semilogx(f,angle(h)*180/pi,'r');
xlabel('频率（Hz）');ylabel('相位');
% -----补偿-----
a = Q/2; aa = a^2;bb = 3*(Q^2)/4;bbb = Q*sqrt(3/4);
b2 = [(1+a*t+aa*t*t)*bbb -2*(1-aa*t*t)*bbb (1-a*t+aa*t*t)*bbb];
a2 = [(aa+bb)*t*2 0 0];

b2 = [1+a*t+aa*t*t -2*(1-aa*t*t) 1-a*t+aa*t*t];
a2 = [(aa+bb)*t*t 0 0];
[h2 f2] = freqz(b2,a2,1000,500);
figure();suptitle('二阶数字补偿时域福频特性曲线');subplot(2,1,1);
semilogx(f2,20*log10(abs(h2)),'r');
xlabel('频率（Hz）');ylabel('幅值(dB)');
subplot(2,1,2);
semilogx(f2,angle(h2)*180/pi,'r');
xlabel('频率（Hz）');ylabel('相位');

