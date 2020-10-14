
% =========================================================================
%
%                  复现轨道检测的算法部分
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.参考：解便滤波器 / 论文
%        2.
%         %% 解偏滤波器
%         %% 解偏滤波器不受到速度的影响
%         %% 为什么？因为这个滤波器并不是随速度变化截止频率的，并不是和前端滤波器一致的
%         %% 这样画图都不收到影响
% 2.1: 因为前端采集信号的特殊性，也是按照时间序列采集然后再选点
%        3. 
%--------------------------------------------------------------------------

close all;
clear all;
wd = 0.001;
b = [1-wd,-1+wd];
a = [1,-1+wd];
% b = [1,2,3];
% a = [4,5];

v = 100/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
figure;
semilogx(lamda,20*log10(abs(h)));


v = 200/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
hold on;
semilogx(lamda,20*log10(abs(h)));
legend 1 2

v = 300/3.6;
t1 = 0.25/v;
[h,f] = freqz(b,a,10000,1/t1);
lamda = v./f;
hold on;
semilogx(lamda,20*log10(abs(h)));
% legend 1 2 3

%% 采样率不是一个问题吗？
lamda = 1:0.001:10000;
psi = 1./lamda;
wd = 0.001;
Ls = 0.25;
Rz = (1-wd)*(1- exp(-j*2*pi*psi*0.25))./( 1-(1-wd)*exp(-j*2*pi*psi*0.25) );
hold on;semilogx(lamda,20*log10(Rz),'LineWidth',1);

%% 模仿freqs画图，不是，另一种表示方法
%%这是啥意思？没看明白，在哪里找到这个公式的
tmp = 2*(1-wd)^2*(1-cos(2*pi.*psi*0.25))./( 1+ (1-wd)^2 -2*(1-wd)*cos(2*pi.*psi*0.25) );
Rz2 = 10*log10(tmp);
%%bode
figure;semilogx(1./psi , 20*log10(Rz),'r--','LineWidth',1);
hold on;semilogx(1./psi , Rz2,'LineWidth',1);

%% 模拟了一个freqz的过程
psi = 0:0.001:900;
fs = 1000;
wd = 0.001;
Ls = 0.25;
Rz = (1-wd)*(1- exp(-j*2*pi*psi/fs))./( 1-(1-wd)*exp(-j*2*pi*psi/fs) );
figure;semilogx(psi,20*log10(Rz),'LineWidth',1);
[h,f] = freqz(b,a,1e5,fs);
hold on;semilogx(f,20*log10(h),'--r','LineWidth',2)


