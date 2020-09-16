clear
close all
global Fs
Fs = 360;
load '118m.mat'%mit数据库第118条数据
signal = val(1,100000:111600)/200;
%% 采用FIR I型设计20Hz以下的低通滤波器
fp=14; fs=18;
detap = 0.01;detas = 0.01;
[M,beta] = selectFirFilterN(fp,fs,detap,detas);
N = M+1;
w = kaiser(N,beta);
hd = FIRItypeIdealpulse(fp,fs,N,'low');
h = hd.*w';
    % 设计的滤波器
omega = linspace(0,pi,512);
mag = freqz(h,[1],omega);
figure
plot(omega/(2*pi)*Fs,20*log10(abs(mag)));
title('FIR低通(14hz以下)滤波器频率相应');
xlabel('频率');
ylabel('增益(dB)');
%% 采用FIR I型设计8Hz以上的高通滤波器
fp2 = 8; fs2=4;
detap2 = 0.01; detas2 = 0.01;
[M2,beta2] = selectFirFilterN(fp2,fs2,detap2,detas2);
N2 = M2+1;
w2 = kaiser(N2,beta2);
hd2 = FIRItypeIdealpulse(fp2,fs2,N2,'high');
h2 = hd2.*w2';
    % 设计的滤波器
omega = linspace(0,pi,512);
mag = freqz(h2,[1],omega);
figure
plot(omega/(2*pi)*Fs,20*log10(abs(mag)));
title('FIR高通(8Hz以上)滤波器频率相应');
xlabel('频率');
ylabel('增益(dB)');
%% 信号滤波
sigFiltered = filter(h,[1],signal);
sigFiltered2 = filter(h2,[1],sigFiltered);
figure
subplot(3,1,1);
plot(signal,'r');
subplot(3,1,2);
hold on
plot(sigFiltered(M/2:end),'b');
subplot(3,1,3);
hold on
plot(sigFiltered2(M/2+M2/2:end),'g');
ylim([-2,2]);
title('信号滤波');
%% 傅里叶变换画出滤波后的频谱
data = FilteredSignal;
M = length(data);
N = M*2-1;
X = fft(data,N);
f = [0:M-1]*Fs/N;
figure
Xabs = abs(fftshift(X));
plot(f(1:end/2),Xabs(M:end-M/2));
title('滤波后的信号频谱');