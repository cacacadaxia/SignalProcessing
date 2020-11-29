


clear all;
close all;

wc = 1;%%低通滤波器，截止频率为wc
b = wc;%%低通
a = [1,wc];
[h_low,f] = freqs(b,a);
b = [1,0];%%高通
[h_high,f] = freqs(b,a);
figure1 = figure('Color',[1 1 1]);
semilogx(f,20*log10(abs(h_low)));hold on;
semilogx(f,20*log10(abs(h_high)));
semilogx(f,20*log10(abs(h_low+h_high)));


set(gca,'Fontname','Times New Roman','fontsize',16);
grid on;