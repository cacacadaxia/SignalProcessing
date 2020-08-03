clear;
close all;
% ----------------------------------------------------------------------
%  step 1 �����ź�
% ------------------------------------------------------------
sr = 100;                       % �����ʣ�1s�ڷ��ŵ�ĸ���(������ʱ��)
t = 0:1/sr:1/sr*(40-1);                  % ��Ժ�����ʱ��
f_sig = [100:0.5:250 ];           % ʵ��
sig = 0;
sigmode = 'real';
switch sigmode
    case 'real'
        for i = 1:length(f_sig)
            sig = sig + sin(2*pi*f_sig(i)*t)*(1+1j);          % ģ��OFDM�ź�
        end
    case 'imag'
        for i = 1:length(f_sig)
            sig = sig + exp(1j*2*pi*f_sig(i)*t);
        end
end
sig = exp(1j*2*pi*30*t)+exp(1j*2*pi*(-25)*t);
f_B = max(f_sig)-min(f_sig);                        % Ƶ�ʣ�ͬʱҲ�Ǵ���B
% figure;plot(t,sig);hold on
% Ϊʲô�����ʺ��źŴ���û��ϵ����Ϊ����sin�ź�
% ���������ͨ�ĵ����źţ�QPSK����BPSK������ô�źŵĴ���B����sr/2�����Բ����Ļ���ֻ��Ҫ2*B=sr�Ϳ��Բ�����
% Ϊʲô����Ϊʵ�ź��븴�ź�֮�������
% ------------------------------------------------
% step 2 ���Բ�ֵ+����
% ------------------------------------------------
fs = 200;
t_cy = 1:1/fs:10;               % ������ʱ��

sig_fj = interp1(t,sig,t_cy,'spline');      % �����ĺ���
% -----------------------------------------------
% hanning windows
sig = conv(hanning(200),sig);

% fft ʵ���źŵ�Ƶ��
NN = length(sig);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig)))),'r-o');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')
% fft �������Ƶ��
% NN = length(sig_fj);figure;plot((-NN/2+1:NN/2)/NN*fs,20*log10(abs(fftshift(fft(sig_fj)))));title('������');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')












