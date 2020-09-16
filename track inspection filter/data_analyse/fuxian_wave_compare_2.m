
% =========================================================================
%
%                  复现轨道检测的算法部分
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.首先复现最为简单的轨向部分
%        2.加上滤波器的部分，并进行对比。主要是轨向
%        3. 
%--------------------------------------------------------------------------

% load_txt;
close all;
clear all;

load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
%% 轨距的对比
% ****************参数设定*********************
delay = 418;%%一直是这个数吗？
G = 9.8;
ht = 0;

tbs = fmctrl_data(:,end);
tbs_s = tbs/1e5;
% ***************step1 模型搭建*************************
rou_l = fmctrl_data(:,9);
rou_r = fmctrl_data(:,11);

rou_l = rou_l/(129.01);
rou_r = rou_r/(129.01);

for i = 3:length(rou_l)
    rou_l_dot2(i) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end


% ******************step2 加速度计的滤波 **********************************
ay = fmctrl_data(:,5);
% ay = ay/(32768/(2*9.83));

% ---------------------- 经过滤波器 --------------------
for i = 1:length(ay)
    ay_Fz(i) = F(ay(i),tbs(i));         %%filter out
    ay_Rz(i) = R(ay_Fz(i),tbs(i));      %%filter out
    ay_Gz(i) = G(ay_Rz(i),tbs(i));      %%filter out
end


figure;plot(ay_Gz);hold on;
plot(ay);legend 1 2


%% 频谱观察
plot_mag(ay,'滤波前')
plot_mag(ay_Gz,'滤波后')


%% 积分
for i = 3:length(sita_b)
    sita_b_dot2 = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz.*tbs.^2 - G .* sita_b .* tbs.^2 + ht*sita_b_dot2;
yL_dot2 = camo - rou_l_dot2;
yR_dot2 = camo + rou_r_dot2;


x_dot = 0;
x = 0;
for i = 3:length(yL_dot2)           %%简单积分，肯定是不对的
    x_dot = x_dot + yL_dot2(i);
    x = x + x_dot;
    yL(i) = x;
end

%% 轨向的波形频谱

aln_l = wave_out(:,3);
aln_r = wave_out(:,4);
plot_mag(aln_l,'左轨向波形');
% figure;plot(x,aln_l,x,aln_r);

function out = F(x,tbs)
%% 滤波器设定
%% 对于x[3]，其与同理
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2
%% 这个滤波器存在延时，几十不等
persistent y;
if isempty(y)
    y = zeros(3,1);
end
y(3) = ( y(2)*(2*2^28+2^14*tbs) - y(1)*2^28+tbs^2*x  )/(2^28 + 2^14*tbs + tbs^2);
%% 更新temp量
y(1) = y(2);
y(2) = y(3);
out = y(3);
end

function out = R(x_k,tbs)
wd = 0.001;

persistent y x;
if isempty(y)
    y = zeros(2,1);
    x = zeros(2,1);
end
x(2) = x_k;
x_dot = x(2)-x(1);
y(2) = (1 - wd) * (x_dot + y(1));

%% 更新temp量
x(1) = x(2);
y(1) = y(2);
out = y(2);

end

function out = G(x_k,tbs_k)
persistent x y1 tbs;
if isempty(x)
%     y = zeros(2,1);
    x = zeros(3,1);
    y1 = zeros(2,1);
    tbs = zeros(2,1);
end
%% 更新k时刻
x(3) = x_k;
tbs(2) = tbs_k;

%% 离散滤波
x_dot = x(3)-x(2);
x_dot_2 = x(3) - 2*x(2) + x(3);
y1(2) = (x(3) - x(2)) *tbs(2)/2^15 + x_dot;

y = (y1(2)+y1(1))* (tbs(2)-tbs(1))/2^16 + x_dot_2;

%% 更新temp量
y1(1) = y1(2);
x(1) = x(2);
x(2) = x(3);
tbs(1) = tbs(2);
out = y;

end

function plot_mag(signal_data , tit)
figure;
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

end










