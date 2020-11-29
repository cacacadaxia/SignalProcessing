
% =========================================================================
%
%                  ���ֹ�������㷨����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��23��(�޸���11.28)
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.��ƫ�˲�����������˵��˲�����һ���ģ��ڿռ����ϲ��ܲ������ٶȵı仯
%        2.��ôʱ�������أ���һ��
% 
% 
%--------------------------------------------------------------------------

%% ��֤��ƫ�˲����ڿռ������ǲ��������ٶȱ仯
%%����϶��ڲ�ͬ�ٶ��¶���һ���ģ���Ϊ��Ӧ�Ĳ�����Ҳ��һ����������󲨳���Ӧ�ľ���һ����
close all;
clear all;
wd = 0.01;
b = [1-wd,-1+wd];
a = [1,-1+wd];
figure1 = figure('Color',[1 1 1]);
v = 100/3.6; t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
semilogx(lamda,20*log10(abs(h)));hold on;
v = 200/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
semilogx(lamda,20*log10(abs(h)));

v = 300/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
hold on;
semilogx(lamda,20*log10(abs(h)));
xlabel('\lambda m');ylabel('Mag dB');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
legend 100km/h 200km/h 300km/h
% 
Ls = 0.25;
psi = f./v;
Rz = (1-wd)*(1- exp(-1j*2*pi*psi*Ls))./( 1-(1-wd)*exp(-1j*2*pi*psi*Ls) );
hold on;semilogx(lamda,20*log10(Rz),'r--','LineWidth',0.5);
%%��������˲��������յ�����ʱ��ĸ��Ŷ��Խ�ֹƵ�ʲ����仯���ڿռ���������
%%��ô���е��˲�������Ҫ�ڿռ�����ȥ�����������ʣ�������˲���������ݱȽ���Ҫ�ĵط�
 
%% �������ݶ��źŽ�����֤�������˲���֮�������0�߸���
lambda1 = 10;
k = 0:1e5;
sig = sin(2*pi*0.25/lambda1.*k);
sig = sig + rand(1,length(sig));
sig = sig + k.*1e-4;
plot_mag(sig,'�����ź�');
sig_filter = filter(b,a,sig);%%�����˲���
plot_mag(sig_filter,'�����ź�','hold');
figure;plot(k./4,sig);
hold on;plot(k./4,sig_filter);
%% ��ʱ�����ϸı������𣿵�Ȼ���ڿռ��򲻱䣬��ô��ʱ�����һ����ı�
    % close all;
    % clear all;
    % wd = 0.01;
    % b = [1-wd,-1+wd];
    % a = [1,-1+wd];
    % v = 100/3.6;
    % t1 = 0.25/v;
    % [h,f] = freqz(b,a,10000,1/t1);
    % lamda = v./f;
    % figure1 = figure('Color',[1 1 1]);
    % semilogx(f,20*log10(abs(h)));
    % 
    % v = 200/3.6;
    % t1 = 0.25/v;
    % [h,f] = freqz(b,a,10000,1/t1);
    % lamda = v./f;
    % hold on;
    % semilogx(f,20*log10(abs(h)));
    % xlabel('f Hz');ylabel('Mag dB');
    % set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
%% �����˲���(ͨ��fdatool��Ƶ��˲���)
% b = load('filter1.mat');
% [h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],10000,4);
% figure1 = figure('Color',[1 1 1]);semilogx(1./f,20*log10(h));grid on;xlabel('���� /m');
% ylabel('��ֵ /dB')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% figure;semilogx(1./f,angle(h));grid on;
%% ����
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
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title(tit);
grid on;
end



