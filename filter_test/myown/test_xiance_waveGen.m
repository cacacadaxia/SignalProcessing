% =========================================================================
%
%                  �����Ҳⷨ���ԵĲ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��24��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���Ҳⷨ�������ַ��Ĵ��ݺ����Ա�
%        2.
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;
del_x = 0.25;
ftime = @(t,v)( 1*(2*1e-3*sin(2*pi/8.*v*t) + 1*1e-3*sin(2*pi/100*v.*t) + 0.7*1e-3*sin(2*pi/12*v.*t) + 1e-3*sin(2*pi/30*v*t) ) );
fspace = @(x)(1*(2*1e-3*sin(2*pi/8 *x) + 1e-3*sin(2*pi/10*x)+0.7*1e-3*sin(2*pi/12*x) + 1e-3*sin(2*pi/30*x) + ...
    1e-3*sin(2*pi/70*x)));
x = 0:del_x:1600;
N = length(x);
longwave = fspace(x);

%% �Ҳ�
L = 10;
xcRes = [];
for i = 1:length(x)-L*4
    xcRes(i) = xiance(longwave,i,L);
end
figure;plot(-xcRes);hold on;
% plot(longwave(20:end));
fspace2 = @(x)(1*(2*1e-3*sin(2*pi/8 *x).*(1-cos(pi*L/8)) + 1e-3*sin(2*pi/10*x).*(1-cos(pi*L/10))+0.7*1e-3*sin(2*pi/12*x).*(1-cos(pi*L/12)) + 1e-3*sin(2*pi/30*x).*(1-cos(pi*L/30)) + ...
    1e-3*sin(2*pi/70*x).*(1-cos(pi*L/70))  ));
xiancedata = fspace2(x);
plot(xiancedata(20:end));
legend �Ҳⷨ ���� �Ҳ�ģ��;
%%��֤�Ҳⷨ���߼��ǲ���ͨ˳�ģ�����ǶԵ�



%% �۲�Ƶ��
% plot_mag(longwave,'ԭʼ����');
% plot_mag(xcRes,'�Ҳⷽ���������','hold');
% legend ԭʼ����Ƶ�� �Ҳ���Ƶ��;

%% 
load('xiancefilter.mat');
xcRes1 = conv(abs(hn),xcRes);
xcRes1(1:400) = [];
% figure;plot(xcRes1(1:end));hold on;plot(longwave(20:end));
% plot_mag(longwave,'ԭʼ����');
% plot_mag(xcRes1,'�Ҳⷽ���������','hold');
% legend ԭʼ����Ƶ�� �Ҳ���Ƶ��

%%
function out = xiance(wave,i,L)
out = (wave(i+L*4)+wave(i))/2 - wave(i+2*L);
end

%%
function plot_mag(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end
fs = 4;     %% 0.25mΪһ���������
signal_data = signal_data./1e3;
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end));
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
title(tit);
grid on;
end