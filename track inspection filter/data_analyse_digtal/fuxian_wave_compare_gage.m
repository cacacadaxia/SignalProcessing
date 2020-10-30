
% =========================================================================
%
%                  ���ֹ�������㷨����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ȸ�����Ϊ�򵥵Ĺ��򲿷�
%        2.
%        3. 
%--------------------------------------------------------------------------

% load_txt;
close all;
clear all;
filepath = 'data/0916_1337_x/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
%% ���ĶԱ�
% ****************�����趨*********************
delay = 418;

% ***************step1 ģ�ʹ*************************
rou_l = fmctrl_data(:,9);
rou_r = fmctrl_data(:,11);

rou_l = rou_l/(129.01);
rou_r = rou_r/(129.01);

gage = wave_out(:,5);
gage = gage/129.01;
sensor_gage = rou_l + rou_r;
sensor_gage = -sensor_gage/2;

offset = -12.8363;
sensor_gage = offset + sensor_gage;
sensor_gage = [zeros(delay,1) ; sensor_gage];
% ------------�����˲�---------
filter_1 = [1/3,1/3,1/3];
sensor_gage = conv(sensor_gage,filter_1);
figure;plot(sensor_gage);
hold on;plot(gage);
legend '������' '���'

% **************Ƶ�׷���****************
signal_data = gage;
figure;
fs = 4;     %% 0.25mΪһ���������
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title('������');

signal_data = sensor_gage;
figure;
fs = 4;     %% 0.25mΪһ���������
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title('��������')
% ***************step2 �˲���*************************

% ***************step3 ���ζԱ�*************************


