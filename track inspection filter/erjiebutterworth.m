% 二阶butterworth抗混滤波器
%幅度谱
a =10^5 / 2^14;
f  = 0.001 : 0.001 :10;  % 模拟频率Hz
omiga = f.* 2*pi;        % 模拟角频率 rad/sec
Fs1 = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);

a = 10^5 / 2^13;
f  = 0.001 : 0.001 :10;  % 模拟频率Hz
omiga = f.* 2*pi;        % 模拟角频率 rad/sec
Fs2 = a^2 ./((1j.*omiga).^2 + 1j.*omiga*a + a^2);

figure
plot( omiga,Fs1,'r',omiga,Fs2,'g');


% figure
% plot( omiga,20*log10(Fs1),'r',omiga,20*log10(Fs2),'g');
% xlabel('Ψ ');
% ylabel('dB ');
% grid on


