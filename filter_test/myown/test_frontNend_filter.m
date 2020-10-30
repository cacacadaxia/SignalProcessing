% =========================================================================
%
%                  测试
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月22日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.
%        2.
%       3. 
%--------------------------------------------------------------------------
clear all
close all

Omega2 = 10^5/2^14;
b = [Omega2^2];
a = [1,Omega2,Omega2^2];
[H,Om] = freqs(b,a);

figure1 = figure('Color',[1 1 1]);
semilogx(Om/2/pi , 20*log10(abs(H)),'LineWidth',2);
hold on
%%
%%修改T的值并不会改变滤波器的性质
T = 1/500;
b = [Omega2^2*T^2];
a = [(1+Omega2*T+(Omega2*T)^2) ,  -(2+Omega2*T) , 1];
[h,f] = freqz(b,a,800000,1/T);
% figure;
semilogx( f , 20*log10(abs(h)),'--r','LineWidth',1);
xlabel('f /Hz')
ylabel('Mag /dB');
set(gca,'Fontname','Times New Roman','fontsize',16);

%%为什么截止频率是1Hz？

%% 空间域？
%%这里的理解不一定对？
v = 100;v = v/3.6;lamda = v./f;
figure;
semilogx( lamda , 20*log10(abs(h)),'--r','LineWidth',1);
hold on;
v = 300;v = v/3.6;lamda = v./f;
semilogx( lamda , 20*log10(abs(h)),'--r','LineWidth',1);
xlabel('波长 /m')
ylabel('Mag /dB');
set(gca,'Fontname','Times New Roman','fontsize',16);




