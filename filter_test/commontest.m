clear all;

% FMIctrl中的滤波器幅频频特性
% ---------- 10 Hz(对于什么？)  -------
fs = 500;
N = 80000;
b10 = [40000 0 0];
a10 = [4010000 -7600000 3610000];
[h10 f10]= freqz(b10,a10,N,'whole',fs);

% % 
mag = 20*log10(abs(h10));    % get magnitude of spectrum in dB
phase = angle(h10)/pi*180;     % get phase in deg.
figure,
subplot(2,1,1),semilogx(f10,mag)
xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
grid on
subplot(2,1,2),semilogx(f10,phase,'r')
xlabel('Frequency (Hz)'),ylabel('Phase (deg.)')
grid on
suptitle('10Hz');

% ----------20 Hz-----------
coef1 = 40000;coef2= 1800000;
coef3=810000 ;coef4=10340000 ;
b20 = [coef1 0 0];
a20 = [coef4 -coef2 coef3];
figure();
[h20 f20]= freqz(b20,a20,N,'whole',fs);
subplot(2,1,1);semilogx(f20,20*log10(abs(h20)));xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
subplot(2,1,2);semilogx(f20,angle(h20)*180/pi);xlabel('Frequency (Hz)'),ylabel('Phase (deg.)')
suptitle('20Hz');grid on;
%-----------1 Hz-------------


%{
Q1 = 2^17;
Q1 = 100000/Q1;
b = [0 Q1];
a = [1 Q1];
[h,f] = freqs(b,a);
%plot((f/(2*pi)),20*log10(abs(h)));
figure();
semilogx((f/(2*pi)),20*log10(abs(h)));
figure();
semilogx((f/(2*pi)),angle(h)*180/pi);
%}
%axis([0 log(0.21) 0 1.1]);
%}
%{
s=tf('s');
Q = 2^14;
Q = 100000/Q;
G= (Q^2)/(s^2 + Q*s+ Q^2); 
bode(G)
%}
%{
Q = 2^14;
Q = 100000/Q;
t = 0.002;
b0 = (Q*t)^2;
a0 = (Q*t)^2+Q*t+1;
a1 = -(Q*t +2);
a2 = 1;
b = [0 0 b0];
a = [a2 a1 a0];
[h,f] = freqz (b,a,1500,500);

figure();
subplot(2,1,1);
%plot(f,20*log10(abs(h)));
%semilogx((f/(2*pi)),20*log10(abs(h)));
semilogx(f,20*log10(abs(h)));
xlabel('频率（Hz）');ylabel('幅值(dB)');

if Q == (100000)/(1900)
    title('10Hz截频滤波器的幅相频特性');
elseif Q == (100000)/(900)
    title('20Hz截频滤波器的幅相频特性');
else 
    title('1Hz截频滤波器的幅相频特性');
end

subplot(2,1,2);%figure (2);
semilogx(f,angle(h)*180/pi);
xlabel('频率（Hz）');ylabel('相位');
%}




Q = 2^14;
Q = 100000/Q;
t = 0.002;
Qt = Q*t;
b = (Q*t)^2;
a0 = b;
a1 = Q;
a2 = 1;


%作图

plot (w,Hf);
plot ((w/2/pi) ,20*log10(abs(Hs)),'color',[.6 .6 .6],'linewidth',3);%hold on;
grid;axis([0 500 -60 5]);hold on;
plot(fz,20 *log10(abs(Hz)),'k');
legend('模拟滤波器','数字滤波器');
xlabel('频率/Hz');ylabel('幅值/db');
title('不进行预畸的数字滤波器与模拟滤波器响应曲线比较');
set(gcf,'color','w');box on;


figure (1);
plot (fz,abs(Hz));hold on;
figure(2);
plot (w,abs(Hs));hold on;

xlabel('频率/Hz');ylabel('幅值/db');
title('巴特沃斯滤波器的幅值响应');
set(gcf,'color','w');



