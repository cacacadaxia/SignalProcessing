


% =========================================================================
%
%                  ����ǰ���˲���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��18��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���Է����źŻ���û�з����ı䣬��ô������ʱ��ȥ����ǰ�˵��˲�����
%        2.��ô���Եõ����ۣ�ǰ���˲������Բ��ã���û�й�ϵ
%       3. 
%--------------------------------------------------------------------------


close all;
clear all;

%% �����ź�

v = 62.5;
Ti = 0.25/v;
fs2 = 1/Ti;
fs = 500;

f1 = 1;
f2 = 237;
f3 = 3;
t1 = 0:1/fs:500-1/fs;
sig = sin(2*pi*f1*t1) + sin(2*pi*f2*t1) +sin(2*pi*f3*t1)  ;

%% �����˲���
signal = sig;
for i = 1:length(signal)
    tbs = 1/fs*1e5;
    sig_B(i) = B(signal(i),tbs);
end
sig2 = sig_B(1:2:end);

for i = 1:length(sig2)
    tbs = 1/fs2*1e5;
    sig_C(i) = C(sig2(i) ,tbs)*fs2/0.76;
end

%% �Ա��ź�
t2 = 0:1/fs2:500-1/fs2;
plot_mag(sig,fs);
plot_mag(sig_B,fs);
plot_mag(sig_C,fs2);

% figure;plot(t1,sig);
% hold on;plot(t2,sig_C);
%%
err = sig(1:2:end) - sig_C;
% figure;plot(err);

%% ����
function plot_mag(sig,fs)
N = length(sig);
x = (-N/2:N/2-1)/N*fs;
figure1 = figure('Color',[1 1 1]);
plot(x,20*log10(abs(fftshift(fft(sig)))));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('f /Hz');
ylabel('Mag /dB');
end

function out = C(x_k,tbs)
omega1 = 0.76;
persistent x;
if isempty(x)
    x = zeros(2,1);
end
%%
x(2) = x_k;
y = x(2)-x(1) + tbs/2^18*(x(2)+x(1));
y = y/omega1;
%%
x(1) = x(2);
out = y;
end

function out = B(x_k,tbs)
persistent y;
if isempty(y)
    y = zeros(2,1);
end
y(2) = ( y(1)*2^17 + tbs*x_k )/(2^17 + tbs );
% tbs*x_k
%%
y(1) = y(2);
out = y(2);
end

