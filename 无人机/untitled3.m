clear;
close all;
% ----------------------------------------------------------------------
%  step 1 生成信号
% ------------------------------------------------------------
sr = 100;                       % 符号率，1s内符号点的个数(决定了时钟)
t = 0:1/sr:1/sr*(40-1);                  % 针对函数的时钟
f_sig = [100:0.5:250 ];           % 实在
sig = 0;
sigmode = 'real';
switch sigmode
    case 'real'
        for i = 1:length(f_sig)
            sig = sig + sin(2*pi*f_sig(i)*t)*(1+1j);          % 模拟OFDM信号
        end
    case 'imag'
        for i = 1:length(f_sig)
            sig = sig + exp(1j*2*pi*f_sig(i)*t);
        end
end
sig = exp(1j*2*pi*30*t)+exp(1j*2*pi*(-25)*t);
f_B = max(f_sig)-min(f_sig);                        % 频率，同时也是带宽B
% figure;plot(t,sig);hold on
% 为什么符号率和信号带宽没关系？因为这是sin信号
% 如果对于普通的调制信号（QPSK或者BPSK），那么信号的带宽B就是sr/2，所以采样的话，只需要2*B=sr就可以采样。
% 为什么？因为实信号与复信号之间的区别
% ------------------------------------------------
% step 2 线性插值+采样
% ------------------------------------------------
fs = 200;
t_cy = 1:1/fs:10;               % 采样的时钟

sig_fj = interp1(t,sig,t_cy,'spline');      % 采样的函数
% -----------------------------------------------
% hanning windows
sig = conv(hanning(200),sig);

% fft 实际信号的频谱
NN = length(sig);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig)))),'r-o');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('幅度谱 dB')
% fft 采样后的频谱
% NN = length(sig_fj);figure;plot((-NN/2+1:NN/2)/NN*fs,20*log10(abs(fftshift(fft(sig_fj)))));title('采样后');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('幅度谱 dB')












