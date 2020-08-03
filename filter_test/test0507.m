clear all;
close all;
t = 0:0.01:3;               % ��ʵ����ʱ��
f1 = 0;                      % Ƶ��
f2 = 200;           
f3 = 0;                    % �趨�������ź�
f4 = 0;

F = @(t)(sin(2*pi*f1*t) + sin(2*pi*f2*t)+ exp(j*2*pi*f3*t) + exp(j*2*pi*f4*t));    % �źź���
y = F(t);                   % �����ź�
% figure;subplot(3,1,1);plot(t , y);         % �ź���ʵͼ
fs = 1000;                  % ������
dtc = 1/fs;                 % �������ʱ��
tc = 0:dtc:4;               % ����ʱ������
yc = F(tc);          % �����ź�����
%% ����Ҷ�任�Լ���ͼ
figure;
N = length(yc);
x = (-N/2+1:N/2)/N*fs;
semilogy(x , abs(fftshift(fft(yc))));


%%
Fs = fs;  % Sampling Frequency

Fstop = 50;              % Stopband Frequency
Fpass = 100;             % Passband Frequency
Dstop = 0.0001;          % Stopband Attenuation
Dpass = 0.057501127785;  % Passband Ripple
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop, Fpass]/(Fs/2), [0 1], [Dstop, Dpass]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});

Hd = dfilt.dffir(b);
		
yf = conv( b , yc);		% �˲�����ź�
%%
figure;plot(abs(3*yf));hold on;plot(abs(yc))

figure;
N = length(yf);
x = (-N/2+1:N/2)/N*fs;
semilogy(x , abs(fftshift(fft(yf))));

%%
N = 512;

num = [c,1];
den = [1,c];
[h1 , ftp] = freqz(num,den,N,Fs);

mag = 20*log10(abs(h1));    % get magnitude of spectrum in dB
phase = angle(h1)/pi*180;     % get phase in deg.

figure,
subplot(2,1,1),plot(ftp,mag)
xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
grid on

subplot(2,1,2),plot(ftp,phase,'r')
xlabel('Frequency (Hz)'),ylabel('Phase (deg.)')
grid on