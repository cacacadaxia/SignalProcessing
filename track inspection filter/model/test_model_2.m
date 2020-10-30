
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
%  功能： 1.
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

%% 一阶差分求和的固定误差如下：
% z = 0;
% for i = 1:length(z_v)
%     z = z + z_v(i)*dt;
%     z_save(i) = z;
% end
% % figure;plot(x,z_save);hold on;plot(x , longwave);legend 1 2
% l1 = z_save-longwave;
%%这里说明了用一阶差分去求和会有固定的误差，大概在1/50左右


%% 陀螺仪
% wavediff_k_1 %%w*del_t
pitch = atan(wavediff);
wy = ( pitch(2:end) - pitch(1:end-1) )/dt;  %%陀螺仪测出来的数据
wy = [0,wy ];                               %%这个影响很多

z_dot = wavediff(1)/4;
z = 0;
for i = 1:length(wy)
    z_dot = z_dot + wy(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end
l2 = z_save - longwave;
% figure1 = figure('Color',[1 1 1]);plot(x,l2*1e3,'r');hold on;plot(x,l1*1e3,'k--');legend 陀螺计算二阶差分 加速度计计算二阶差分;
xlabel('采样点 /0.25m');ylabel('误差 /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);
%%从理论上来说这两者应该是相等的，没有太多的差别的
%%但是这里产生的差别就是积分带来的误差？
%%这个可以作为用来标定的参考方法

%% 用两个测距的组件
f = @(x)(waveMag*1e-3*sin(2*pi/waveLen*x));

for i = 1:length(x)
    [pitch2(i) , out2(i)] = cau(x(i),f);
end
% figure1 = figure('Color',[1 1 1]);plot(pitch2);hold on;plot(pitch(7:end));legend 1 2;
%%这说明了什么？这说明了"存在3m基长的轨道的切角与实际的切角误差"非常小，可以忽略不计
%%只不过两者之间有延时。为什么，因为实际上也是切线，是(z1+z2)/2的切线




%% 
% (z1+z2)/2
f_tp = @(x)(( waveMag*1e-3*sin(2*pi/waveLen*x) + waveMag*1e-3*sin(2*pi/waveLen*(x+L)) )/2);
wavediff_tp = (waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen*x) + waveMag*1e-3*2*pi/waveLen*cos(2*pi/waveLen*(x+3)))/2;
% figure;plot((atan(wavediff_tp) - pitch2)/pi*180);
%%这里说明了实际上测的是(z1+z2)/2的结果

%% 测试新的方案
det1 = sin(pitch2)*L;

num = 0;
% det1 = [zeros(1,num),det1(1:end-num)];
% det1 = [det1(1+num:end),zeros(1,num)];
% wavediff_k_1 %%w*del_t

wy2 = ( pitch2(2:end) - pitch2(1:end-1) )/dt;%%陀螺仪测出来的数据
wy2 = [0 , wy2];

% z_dot = wavediff_tp(1)/4;%%只是这个还不行
z_dot = 2.093017317643780e-04;
z = 0;
for i = 1:length(wy2)
    z_dot = z_dot + wy2(i) * dt * del_x;
    z = z + z_dot;
    z_dot_save(i) = z_dot;
    z_save(i) = z;
end
det2 = z_save * 2;
y1 = (det2 - det1)/2;
l3 = y1 - longwave;
% figure1 = figure('Color',[1 1 1]);plot(x,(y1 - longwave)*1e3+1.05);
xlabel('采样点 /0.25m');ylabel('误差 /mm');
set(gca,'Fontname','Times New Roman','fontsize',16);
%%这里有幅值的差别，明显是存在问题的，所以需要重新测量


%% 对比所有的方案
% y1 = (y1(2:end)+y1(1:end-1))/2;
% y1 = [0,y1];
figure1 = figure('Color',[1 1 1]);plot( x , l3*1e3 + 1.25);
hold on;
plot(x,l2*1e3,'r');hold on;plot(x,l1*1e3,'k--');legend 两组测距组件  陀螺计算二阶差分 加速度计计算二阶差分;
xlabel('采样点 /0.25m');ylabel('误差 /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);

%%经过调整之后，发现其误差与一阶差分求和的结果是一样的
%%这样对比虽然产生了一定的误差结果，但是总感觉有点问题，怎么才能和长波联系上呢？




%%
figure1 = figure('Color',[1 1 1]);plot(x(1:end-6),(z_save(1:end-6) - longwave(7:end))*1e3);
figure;plot(z_save);hold on;plot(longwave);legend 1 2


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
