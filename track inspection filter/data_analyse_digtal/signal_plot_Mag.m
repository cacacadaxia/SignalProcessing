% =========================================================================
%
%                  查看plotMag的方法对不对
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月10日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.对比不同的长波滤波器的副频改变的区别
%        2.
%       3. 
%--------------------------------------------------------------------------



close all;
clear all;
N = 1e5;
signal = randn(N,1)*200;
fs = 4;

P = mean(abs(signal).^2);
E = sum(abs(signal).^2);
Mag = 10*log10(E);

%%
sig_2 = longwave_filter(signal,281,71,281,491);
plot_mag(signal,'滤波前');
plot_mag(sig_2,'滤波后');


b = load('filter1.mat');
sig_3 = conv(signal,b.Num);
plot_mag(sig_3,'fdatool滤波后');

b = load('filter_400m.mat');
sig_4 = conv(signal,b.Num);
plot_mag(sig_4,'fdatool滤波后');


b = load('filter_1000m.mat');
[h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],10000,4);
figure1 = figure('Color',[1 1 1]);semilogx(1./f,20*log10(h));grid on;xlabel('波长 /m');
ylabel('幅值 /dB')
set(gca,'Fontname','Times New Roman','fontsize',16);


%%
% %% 25m长波滤波加积分
% alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
% sscal = 0.000825;
% sbsci = 0.019802;
% fscal = 0.1;
% sbsc = 101.000;
% % 数组设定
% Num = 768;
% amcol_array = zeros(Num,1);
% amcol_arraytmp = zeros(Num,1);
% in1 = 533; in2 = 432;in4=382;in6=331;in7=230;
% in = 539;
% 
% for i = 1:length(signal)           %%简单积分，肯定是不对的
%     %%
% %     amcol = aln(:,3);
%    amcol_array(in) = signal(i);
%     alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
%     elupp = alu;
%     elup = elup + elupp;
%     elu = elu + elup;
%     emco = - amcol_array(in4);
%     als = als + amcol_array(in2) - amcol_array(in6);
%     alss = alss + als;
%     alss = alss + sbsc*emco;
%     alsss = alsss + alss;
%     xtemp = (alsss*sbsci - sscal*elu)*fscal;
%     yL(i,1) = xtemp;
%     
%     %% 更新数组
%     in1 = mod(in1,Num)+1;
%     in2 = mod(in2,Num)+1;
%     in4 = mod(in4,Num)+1;
%     in6 = mod(in6,Num)+1;
%     in7 = mod(in7,Num)+1;
%     in = mod(in,Num)+1;
%     %%
%     save(i,1) = alu;
%     save(i,2) = elu;
%     save(i,3) = als;
%     save(i,4) = alsss;
%      
% end
% plot_mag(yL,'滤波后');



%% function
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
set(gca,'Fontname','Times New Roman','fontsize',14);
title(tit);
grid on;
end

function plot_mag2(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure;
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