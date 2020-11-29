

% =========================================================================
%
%                  长波检测方法讨论
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 11月25日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.调整一下不平顺，让其变得复杂一些
%        2.陀螺主要收到白噪声大小的影响，但是加速度计收到速度影响比较大(JC-21数字化惯性导航系统技术指标JC-21数字化惯性导航系统技术指标.docx)
%        3. 这里并没有讨论"两组测距组件"对于结果的影响
% 
% 11.25
%       1.验证方案的可行性（陀螺替代的方法+求俯仰角的方法）
% 
%--------------------------------------------------------------------------

clear all;
close all;
% lamda1 = 8;
del_x = 0.25;
ftime = @(t,v)(1*(2*1e-3*sin(2*pi/8.*v*t)+1*1e-3*sin(2*pi/10*v.*t)+0.7*1e-3*sin(2*pi/12*v.*t) + 1e-3*sin(2*pi/30*v*t)));
fspace = @(x)(1*(2*1e-3*sin(2*pi/8 *x) + 1e-3*sin(2*pi/10*x)+0.7*1e-3*sin(2*pi/12*x) + 1e-3*sin(2*pi/30*x)));
ftime = @(t,v)(1*(2*1e-3*sin(2*pi/8.*v*t) ));
fspace = @(x)(1*(2*1e-3*sin(2*pi/8 *x) ));
x = 0:del_x:1600;
N = length(x);
longwave = fspace(x);
% figure1 = figure('Color',[1 1 1]);plot(x,longwave);set(gca,'Fontname','Times New Roman','fontsize',14);

%%
syms x_syms t_syms v_syms;
f_tp = diff(fspace(x_syms));
fspace_dot = @(x_syms)eval(f_tp);
wavediff = fspace_dot(x);
f_tp = diff(ftime(t_syms,v_syms) , t_syms);
f_z_v = @(t_syms,v_syms)eval(f_tp);
f_tp2 = diff(f_tp , t_syms);
f_z_acc = @(t_syms,v_syms)eval(f_tp2);

v = 100/3.6;
dt = del_x/v;                                       %%采样时间间隔
t = 0:dt:(length(x)-1)*dt;
z_v = f_z_v(t,v);

%% 加噪声
sigma_acc = 1e-4 * 9.8  / sqrt(dt);
sigma_acc = 0;
z_acc = f_z_acc(t,v);
% z_acc = z_acc + randn(1,length(t)) * sigma_acc;
% figure;plot(z_acc);hold on;plot( f_z_acc(t,v))
%%对acc加噪声完全不行


%% 加速度进行积分
z_dot =  wavediff(1)/4;%%有点飘
z_dot = 6.921463910193095e-04 - 0.0000001;
for i = 1:length(z_acc)
    z_save(i) = double_integral(z_acc(i)*dt*dt,z_dot);
end

% figure;plot(x(3:end),z_save(3:end)*1e3);hold on;plot(x , longwave*1e3);legend 1 2
figure;plot(z_save(1:end)*1e3);hold on;plot(longwave(1:end)*1e3);legend 1 2
l1 = z_save - longwave;
% figure;plot(l1*1e3);

%% 陀螺仪进行积分
% wavediff_k_1 %%w*del_t
% pitch = atan(wavediff);
pitch = (longwave(2:end)-longwave(1:end-1))/0.25;
wy = ( pitch(2:end) - pitch(1:end-1) )/dt;  %%陀螺仪测出来的数据
wy = [0,wy ];                               %%这个影响很多

sigma_gyro = 1/180*pi/3600 / sqrt(dt);
sigma_gyro = 0;
wy = wy + randn(1,length(wy)) *sigma_gyro;
z_dot = pitch(1)/4;
z_dot = pitch(1)/4;
z = 0;
for i = 1:length(wy)
    z_dot = z_dot + wy(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;%%无论条不调整延时都存在误差
end
l2 = z_save(1:end-1) - longwave(2:end);
figure;plot(z_save*1e3);hold on;plot( longwave(2:end)*1e3);legend 1 2

figure1 = figure('Color',[1 1 1]);plot(l2*1e3,'r','LineWidth',0.5);hold on;plot(l1*1e3,'k');legend 陀螺计算二阶差分 加速度计计算二阶差分;
grid on;
xlabel('采样点 /0.25m');ylabel('误差 /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);
%%从理论上来说这两者应该是相等的，没有太多的差别的
%%但是这里产生的差别就是积分带来的误差？
%%这个可以作为用来标定的参考方法

%%这里是不是有问题，怎么误差又变得这么小？
% 可能是因为模型相关没想清楚

%% 角速度计算
% wavediff;
%%半个点的延时
% pit = (longwave(2:end)-longwave(1:end-1))/0.25;
% figure;plot(pit(1:end));hold on;
% plot(wavediff(1:end))
% legend 1 2;

%% 新的方案
L = 5;
pitch2 = (fspace(x + L) - fspace(x))/L;
wy2 = ( pitch2(2:end) - pitch2(1:end-1) )/dt;%%陀螺仪测出来的数据
wy2 = [0 , wy2];
figure;plot(pitch2);hold on;plot(pitch(2*5:end)/2)

list = (pitch2(2:end)-pitch2(1:end-1))*0.25;
dottmp = pitch2(1)/4;
for i = 1:length(list)
    res(i) = double_integral2(list(i),dottmp);
end
Res2 = (fspace(x) + fspace(x+L))/2;
Res2 = fspace(x) - fspace(x+L);
figure;plot(res+0.002);
hold on;
plot(Res2)
legend 1 2

%% 观察频谱
plot_mag(longwave,'原始波形');
plot_mag(res,'新方案','hold');
plot_mag(Res2,'z1+z2/2');




%% 处理过延时的问题去积分
%%积分函数没有考虑延时的问题
function out = double_integral(x_k , x_dot_0)
persistent ydot ydot_1 y_1;
if isempty(ydot)
    ydot = 0;
    ydot_1 = x_dot_0;%%初值很重要
    y_1 = 0;
end
ydot = ydot_1 + x_k;
y_k = y_1 + ydot;
out = y_1;
%%
y_1 = y_k;
ydot_1 = ydot;
end

function out = double_integral2(x_k , x_dot_0)
persistent ydot ydot_1 y_1;
if isempty(ydot)
    ydot = 0;
    ydot_1 = x_dot_0;%%初值很重要
    y_1 = 0;
end
ydot = ydot_1 + x_k;
y_k = y_1 + ydot;
out = y_1;
%%
y_1 = y_k;
ydot_1 = ydot;
end

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
% N = 1024
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data,N)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title(tit);
grid on;
end
