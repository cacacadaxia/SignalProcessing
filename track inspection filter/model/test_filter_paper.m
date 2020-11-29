
% =========================================================================
%
%                  测试论文中的算法
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.
%        2.
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;


%% 截止频率
wc = 0.76;%%低通滤波器，截止频率为wc
b = wc;%%低通
a = [1,wc];
[h_low,f] = freqs(b,a);
b = [1,0];%%高通
[h_high,f] = freqs(b,a);
figure1 = figure('Color',[1 1 1]);
semilogx(f,20*log10(abs(h_low)));hold on;
semilogx(f,20*log10(abs(h_high)));
semilogx(f,20*log10(abs(h_low+h_high)));

set(gca,'Fontname','Times New Roman','fontsize',14);
grid on;
xlabel('\psi Hz');ylabel('Mag dB');

%% 1/s
[h,f] = freqs(1,[1,0]);
figure1 = figure('Color',[1 1 1]);
semilogx(f,20*log10(abs(h)));hold on;
set(gca,'Fontname','Times New Roman','fontsize',14);
grid on;
xlabel('\psi Hz');ylabel('Mag dB');

