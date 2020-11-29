
% =========================================================================
%
%                  矩形窗
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.窗函数的画图，参考文档「高速线路轨道长波不平顺检测技术研究.docx」
%        2.这个很重要，解释了各种符号之间的关系
%        3. 几种画图的方法，非常重要
% 
%   1102：
%   1.长波滤波器+积分
%    11.28:
%   1.加上美国的长波滤波器，未完成
%--------------------------------------------------------------------------

close all;
clear all;

%%
figure1 = figure('Color',[1 1 1]);
M = 30;
N = 2*M+1;
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %空间域频率
Omega = 2*pi*pesi*0.25;
Wn = sin(N.*Omega/2)./sin(Omega./2)./N;
semilogx(lamda,20*log10(abs(Wn)));

%% freqz
b = ones(1,N)/N;
[h,f] = freqz(b,[zeros(1,N-1),1],10000,4);
hold on;semilogx(1./f,20*log10(h));
grid on;

%% z
K = 30;
M = 2*K+1;
% lamda = 1:0.1:1000;
% pesi  = 1./lamda;           %空间域频率
pesi = 0.001:0.0001:1;
lamda = 1./pesi;
temp = exp(1j*2*pi.*pesi*0.25);
%%换一种表示，就是两个矩形窗
DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
           ./( ( 1- temp.^(-1) ) ) ...
            /M ;
% hold on;
semilogx( lamda , 20*log10(DjpesiDz) ,'LineWidth',1);
xlabel('\lambda /m');ylabel('Mag /dB')
set(gca,'Fontname','Times New Roman','fontsize',16);

%% 长波滤波器（美国）(有点问题)
k = 70;
m = floor(0.829*k);
n = floor(0.646*k);
z_1 = temp.^(-1);
Tz = 1-1/k/m/n.*(1-z_1.^k).*(1-z_1.^m).*(1-z_1.^n)./(1-z_1).^3;
figure1 = figure('Color',[1 1 1]);
semilogx( lamda , 20*log10(Tz) ,'LineWidth',1);
xlabel('\lambda /m');ylabel('Mag /dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
grid on;




%% 二阶差分组合
% H_integral = 1./( 1 - 2*temp.^(-1) + temp.^(-2) );
% Hz = DjpesiDz.*H_integral;
% hold on;grid on;
% % semilogx(lamda , 20*log10(Hz),'LineWidth',1);
% legend ;
%% fdatool
% b = load('filter1.mat');
% [h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],length(lamda)*2,4);
% figure1 = figure('Color',[1 1 1]);semilogx(1./f(1:length(f)/2),20*log10(h(1:length(f)/2)));grid on;xlabel('波长 /m');
% % figure1 = figure('Color',[1 1 1]);semilogx(1./f(length(f)/2+1:end),20*log10(h(length(f)/2+1:end)));grid on;xlabel('波长 /m');
% ylabel('幅值 /dB')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% Hz2 = h(1:length(f)/2).*H_integral.';
% hold on;
% semilogx(lamda , 20*log10(abs(Hz2)));
% %%这里得到结论：信号二阶差分的截止波长与信号本身截止波长是一回事。对二阶差分进行滤波就相当于对信号进行滤波。
% %%所以在分析信号经过长波滤波器的频谱时，发现有这样的结果。




