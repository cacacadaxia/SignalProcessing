
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
%        2.加上滤波器的部分，并进行对比
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
delay = 418;



% ***************step1 模型搭建*************************
rou_l = fmctrl_data(:,9);
rou_r = fmctrl_data(:,11);

rou_l = rou_l/(129.01);
rou_r = rou_r/(129.01);

gage = wave_out(:,5);
gage = gage/129.01;
sensor_gage = rou_l + rou_r;
sensor_gage = -sensor_gage/2;

offset = -12.8363;
sensor_gage = offset + sensor_gage;
sensor_gage = [zeros(delay,1) ; sensor_gage];
% ------------三点滤波---------
filter_1 = [1/3,1/3,1/3];
sensor_gage = conv(sensor_gage,filter_1);
figure;plot(sensor_gage);
hold on;plot(gage);
legend '传感器' '轨距'


