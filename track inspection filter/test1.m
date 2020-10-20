clear all;
close all;
B = [1,2];
A = [3,4];
[ H , tmp ] = freqs(B, A);%# 连续的模拟信号
figure;
semilogx(tmp/pi/2 , 20*log10(H),'--','LineWidth',1);
set(gca,'Fontname','Times New Roman','fontsize',16);
xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
f = 0.1:0.001:100;
Omega = 2*pi*f;
Cs = (1j*Omega*B(1) + B(2) )./( 1j*Omega*A(1) + A(2) );
hold on
semilogx( f , 20*log10(Cs),'LineWidth',1 );
xlabel('模拟频率 (Hz)');
ylabel('dB');   