


% ̽��bodeͼ�Ļ�ͼ��ʽ
clear, close all

%% initialize parameters
% �ز�Ƶ��
samplerate = 1000;   % in Hz ������
N = 512;             % number of points, must be even, better be power of 2

%% define a and b coeffients of H (transfer function)
a = [1 -0.2 0.8];   % denominator terms
b = [0.2 0.5 0];      % numerator terms

%% �������棺
N = 512;
[h1 , ftp] = freqs(b,a,N);

mag = 20*log10(abs(h1));    % get magnitude of spectrum in dB
phase = angle(h1)/pi*180;     % get phase in deg.
figure,
subplot(2,1,1),semilogx(ftp,mag)
xlabel('Frequency (Hz)'),ylabel('Magnitude (dB)')
grid on
phase(1:330) = phase(1:330)-360;
subplot(2,1,2),semilogx(ftp,phase,'r')
xlabel('Frequency (Hz)'),ylabel('Phase (deg.)')
grid on

%%
figure;bode(b,a)
