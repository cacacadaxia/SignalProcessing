clear;close all;
% ----------------------------------------------------------------------
%  step 1 生成信号
% ------------------------------------------------------------
sr = 1000;                       % 符号率，1s内符号点的个数(决定了时钟)
t = 1:1/sr:10;                  % 针对函数的时钟
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

f_B = max(f_sig)-min(f_sig);                        % 频率，同时也是带宽B
% figure;plot(t,sig);hold on
% 为什么符号率和信号带宽没关系？因为这是sin信号
% 如果对于普通的调制信号（QPSK或者BPSK），那么信号的带宽B就是sr/2，所以采样的话，只需要2*B=sr就可以采样。
% 为什么？因为实信号与复信号之间的区别
% ------------------------------------------------
% step 2 线性插值+采样
% ------------------------------------------------
t_BS = 20;                       % 模拟采样偏差，时钟偏差(前后的采样时钟不准确)(这个具体怎么理解？)
fs   = t_BS*f_B;                % 采样频率 >2B
% fs   = 500;                     % 采样频率
t_cy = 1:1/fs:10;               % 采样的时钟
sig_fj = interp1(t,sig,t_cy,'spline');      % 采样的函数
% -----------------------------------------------
% fft 实际信号的频谱
NN = length(sig);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig)))));axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('幅度谱 dB')
% fft 采样后的频谱
NN = length(sig_fj);figure;plot((-NN/2+1:NN/2)/NN*fs,20*log10(abs(fftshift(fft(sig_fj)))));title('采样后');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('幅度谱 dB')
% ------------------------------------------------
% 频谱搬移

delta_f = 100;                   % 需要搬移的频率，也可以认为是频偏(正表示向右边搬移)
sig_by = sig.*exp(1j*2*pi*delta_f/sr.*(1:length(sig)));     % 频谱搬移
NN = length(sig_by);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig_by)))));
% ----------------------------------------------------------------------
% FPGA对于样点的处理
% ------------------------------------------------------------
% 成型滤波 
alfa  = 0.35;               % 滚降
delay = 6;                  % 延时
Ns    = 4;                     % 滤波成型倍数
fs2   = Ns*sr;                % 另一种采样
HopSendOutRcos = rcosflt(sig,sr,fs2,'fir/sqrt',alfa,delay); % 升余弦成型
HopSendOutRcos = HopSendOutRcos.';
NN = length(HopSendOutRcos);figure;plot((-NN/2+1:NN/2)/NN*fs2,20*log10(abs(fftshift(fft(HopSendOutRcos)))));title('滤波成型后');axis([-inf,inf,-inf,inf])
% 上采样，产生周期性延展，相当于把信号变成离散的，离散信号的频谱是周期性延展的
sigUpsam = upsample(sig,Ns);        
NN = length(sigUpsam);figure;plot((-NN/2+1:NN/2)/NN*fs2,20*log10(abs(fftshift(fft(sigUpsam)))));title('滤波成型后');axis([-inf,inf,-inf,inf])
% 
% % 滤波后与原信号重合
% % -------------------------------------------------------------
% 接收机部分的ADC芯片对模拟信号在采样上的偏差
offset = 1*sr/4;
samplefre = fs2 + offset;
time1 =fs2/samplefre.*(1:length(HopSendOutRcos));
signalsample = interp1((1:length(HopSendOutRcos)),HopSendOutRcos,time1,'spline');  
NN = length(signalsample);figure;plot((-NN/2+1:NN/2)/NN*fs2,20*log10(abs(fftshift(fft(signalsample)))));title('加上ADC采样偏差后')
