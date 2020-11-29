
% =========================================================================
%
%                  复现轨道检测的算法部分
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月23日(修改于11.28)
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.解偏滤波器与其他后端的滤波器是一样的，在空间域上不受采样的速度的变化
%        2.那么时间域上呢？不一致
% 
% 
%--------------------------------------------------------------------------

%% 验证解偏滤波器在空间域上是不是随着速度变化
%%这里肯定在不同速度下都是一样的，因为对应的采样率也不一样，所以最后波长对应的就是一样的
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
%%这个数字滤波器并不收到采样时间的干扰而对截止频率产生变化，在空间域上来看
%%那么所有的滤波器都需要在空间域上去检验它的性质，这才是滤波器相关内容比较重要的地方
 
%% 生成数据对信号进行验证，加上滤波器之后基本在0线附近
lambda1 = 10;
k = 0:1e5;
sig = sin(2*pi*0.25/lambda1.*k);
sig = sig + rand(1,length(sig));
sig = sig + k.*1e-4;
plot_mag(sig,'正弦信号');
sig_filter = filter(b,a,sig);%%加上滤波器
plot_mag(sig_filter,'正弦信号','hold');
figure;plot(k./4,sig);
hold on;plot(k./4,sig_filter);
%% 在时间域上改变特性吗？当然，在空间域不变，那么在时间域就一定会改变
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
%% 长波滤波器(通过fdatool设计的滤波器)
% b = load('filter1.mat');
% [h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],10000,4);
% figure1 = figure('Color',[1 1 1]);semilogx(1./f,20*log10(h));grid on;xlabel('波长 /m');
% ylabel('幅值 /dB')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% figure;semilogx(1./f,angle(h));grid on;
%% 函数
function plot_mag(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end
fs = 4;     %% 0.25m为一个采样间隔
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



