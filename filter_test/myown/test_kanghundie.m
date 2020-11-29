

% =========================================================================
%
%                  ���Կ�����˲����Ŀ������Ч��
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��20��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.Ŀǰ��û���ܹ��ҵ�˵������ķ���
%        2.
%       3. 
%--------------------------------------------------------------------------



clear all;
close all;
f1 = 10;
f2 = 220;   %%�����
fs = 500;
t = 0:1/fs:1000;
sig1 = sin(2*pi*f1*t) + sin(2*pi*f2*t);
figure;plot(t,sig1);

N = length(t);
x = (-N/2:N/2-1)/N*fs;
figure;plot(x,20*log10(abs(fftshift(fft(sig1)))));
grid on;

%% 
b = load('filter_test.mat');%%180Hz
L = (length(b.Num)-1)/2;
sig2 = conv(b.Num,sig1);
sig2(1:L) = [];
sig2(end-L+1:end) = [];

N = length(t);
x = (-N/2:N/2-1)/N*fs;
figure;plot(x,20*log10(abs(fftshift(fft(sig2)))));
grid on;

%%���˵��

%%
sig3 = sig2(1:2:end);

N = length(sig3);
x = (-N/2:N/2-1)/N*fs/2;
figure;plot(x,20*log10(abs(fftshift(fft(sig3)))));
grid on;

%%
% fs = 250;
% t = 0:1/fs:1000;
% sig2 = sin(2*pi*f1*t) + sin(2*pi*f2*t);
% sig3 = sin(2*pi*f1*t);
% 
% figure;plot(sig2-sig3);

%% function

