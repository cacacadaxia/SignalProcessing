	

clear, close all

%% initialize parameters
% 载波频率
samplerate = 1000;   % in Hz 采样率
N = 512;             % number of points, must be even, better be power of 2

%% define a and b coeffients of H (transfer function)
a = [1 -0.2 0.8];       % denominator terms
b = [0.2 0.5 0];      % numerator terms

%% 或者下面：
N = 512;
[h1 , ftp] = freqz(b,a,N,samplerate);

mag = 20*log10(abs(h1));    % get magnitude of spectrum in dB
phase = angle(h1)/pi*180;     % get phase in deg.
figure,
subplot(2,1,1),semilogx(ftp,mag)
xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
grid on
subplot(2,1,2),semilogx(ftp,phase,'r')
xlabel('Frequency (Hz)'),ylabel('Phase (deg.)')
grid on

% N = 512;
% 	a = 1;
% 	H = fft(b,N)/fft(a,N);   % H矩阵
% 	mag = 20*log10(abs(H));    % get magnitude of spectrum in dB 幅值
% 	phase = angle(H)*2*pi;     % get phase in deg.相位
% 	faxis = samplerate/2*linspace(0,1,N/2);  % the axis of frequency
% 	%% plot the spectrum of H
% 	figure,
% 	subplot(2,1,1),plot(faxis,mag(1:N/2))
% 	xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
% 	grid on
% 	subplot(2,1,2),plot(faxis,phase(1:N/2),'r')
% 	xlabel('Frequency (Hz)'),ylabel('Phase (deg.)')
% 	grid on