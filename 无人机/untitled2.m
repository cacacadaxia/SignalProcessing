clear;close all;
% ----------------------------------------------------------------------
%  step 1 �����ź�
% ------------------------------------------------------------
sr = 1000;                       % �����ʣ�1s�ڷ��ŵ�ĸ���(������ʱ��)
t = 1:1/sr:10;                  % ��Ժ�����ʱ��
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

f_B = max(f_sig)-min(f_sig);                        % Ƶ�ʣ�ͬʱҲ�Ǵ���B
% figure;plot(t,sig);hold on
% Ϊʲô�����ʺ��źŴ���û��ϵ����Ϊ����sin�ź�
% ���������ͨ�ĵ����źţ�QPSK����BPSK������ô�źŵĴ���B����sr/2�����Բ����Ļ���ֻ��Ҫ2*B=sr�Ϳ��Բ�����
% Ϊʲô����Ϊʵ�ź��븴�ź�֮�������
% ------------------------------------------------
% step 2 ���Բ�ֵ+����
% ------------------------------------------------
t_BS = 20;                       % ģ�����ƫ�ʱ��ƫ��(ǰ��Ĳ���ʱ�Ӳ�׼ȷ)(���������ô��⣿)
fs   = t_BS*f_B;                % ����Ƶ�� >2B
% fs   = 500;                     % ����Ƶ��
t_cy = 1:1/fs:10;               % ������ʱ��
sig_fj = interp1(t,sig,t_cy,'spline');      % �����ĺ���
% -----------------------------------------------
% fft ʵ���źŵ�Ƶ��
NN = length(sig);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig)))));axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')
% fft �������Ƶ��
NN = length(sig_fj);figure;plot((-NN/2+1:NN/2)/NN*fs,20*log10(abs(fftshift(fft(sig_fj)))));title('������');axis([-inf,inf,-inf,inf]);xlabel('fs Hz');ylabel('������ dB')
% ------------------------------------------------
% Ƶ�װ���

delta_f = 100;                   % ��Ҫ���Ƶ�Ƶ�ʣ�Ҳ������Ϊ��Ƶƫ(����ʾ���ұ߰���)
sig_by = sig.*exp(1j*2*pi*delta_f/sr.*(1:length(sig)));     % Ƶ�װ���
NN = length(sig_by);figure;plot((-NN/2+1:NN/2)/NN*sr,20*log10(abs(fftshift(fft(sig_by)))));
% ----------------------------------------------------------------------
% FPGA��������Ĵ���
% ------------------------------------------------------------
% �����˲� 
alfa  = 0.35;               % ����
delay = 6;                  % ��ʱ
Ns    = 4;                     % �˲����ͱ���
fs2   = Ns*sr;                % ��һ�ֲ���
HopSendOutRcos = rcosflt(sig,sr,fs2,'fir/sqrt',alfa,delay); % �����ҳ���
HopSendOutRcos = HopSendOutRcos.';
NN = length(HopSendOutRcos);figure;plot((-NN/2+1:NN/2)/NN*fs2,20*log10(abs(fftshift(fft(HopSendOutRcos)))));title('�˲����ͺ�');axis([-inf,inf,-inf,inf])
% �ϲ�����������������չ���൱�ڰ��źű����ɢ�ģ���ɢ�źŵ�Ƶ������������չ��
sigUpsam = upsample(sig,Ns);        
NN = length(sigUpsam);figure;plot((-NN/2+1:NN/2)/NN*fs2,20*log10(abs(fftshift(fft(sigUpsam)))));title('�˲����ͺ�');axis([-inf,inf,-inf,inf])
% 
% % �˲�����ԭ�ź��غ�
% % -------------------------------------------------------------
% ���ջ����ֵ�ADCоƬ��ģ���ź��ڲ����ϵ�ƫ��
offset = 1*sr/4;
samplefre = fs2 + offset;
time1 =fs2/samplefre.*(1:length(HopSendOutRcos));
signalsample = interp1((1:length(HopSendOutRcos)),HopSendOutRcos,time1,'spline');  
NN = length(signalsample);figure;plot((-NN/2+1:NN/2)/NN*fs2,20*log10(abs(fftshift(fft(signalsample)))));title('����ADC����ƫ���')
