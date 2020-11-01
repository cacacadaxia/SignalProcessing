
% =========================================================================
%
%                  计算差分误差
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月20日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.探讨一下积分的问题
%        2.
%        3. 
%--------------------------------------------------------------------------


%% 生成数据
clear all;
close all;
global L;
L = 3;
del_x = 0.25;               %%单位为m
waveLen = 300;              %%单位为m
x = 0:del_x:waveLen*10;     %%空间域序列
waveMag = 40;               %%mm 长波的幅值
longwave = waveMag*1e-3*sin(2*pi/waveLen*x);
% figure;plot(x,longwave);xlabel('采样间隔 /0.25m');ylabel('幅值 /mm');

%% 为什么去找长波的切线呢？
%%因为轨面曲线的切线就象征着构架计算出来的点头角速度
wavediff = waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen* x );

%% 时间序列
v = 100/3.6;                                             %%m/s
dt = del_x/v;                                       %%采样时间间隔
t = 0:dt:(length(x)-1)*dt;
T = dt*waveLen/del_x;                               %%时间采样周期
z_true = waveMag*1e-3*sin(2*pi/T*t);                %%时间序列
z_v = waveMag*1e-3*2*pi/T*cos(2*pi/T*t);            %%速度序列
z_acc = -waveMag*1e-3*2*pi/T*2*pi/T*sin(2*pi/T*t);  %%加速度计的值


%% 加速度计
z_dot =  wavediff(1)/4;
z = 0;
for i = 1:length(z_acc)
    z_dot = z_dot + z_acc(i)*dt*dt;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end

% figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
l1 = z_save-longwave;



%% 梯形积分
%%这里的处理可能有点问题
z_dot =  wavediff(1)/4;
z = 0;
for i = 1:length(z_acc)-1
    z_dot = z_dot + (z_acc(i) + z_acc(i+1))/2*dt*dt;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end

% figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
l2 = z_save - longwave;

%%
figure1 = figure('Color',[1 1 1]);plot(l1*1e3);hold on;plot(l2*1e3);
legend 一阶积分 梯形积分
set(gca,'Fontname','Times New Roman','fontsize',14);
%%更换积分的方法确实会对结果产生影响，暂时不去分析


%% function
%%这个函数比较耗时
function out = cauDis(x1,x2)
    out = norm(x1 - x2);
end
function [out,out2 ] = cau(x_pos, f)
global L;

Dis = L;
Fx = [x_pos+Dis , f(x_pos+Dis)];
Bx = [x_pos , f(x_pos)];
out2 = cauDis(Fx,Bx);
out = (f(x_pos + Dis) - f(x_pos)) / Dis;
end
