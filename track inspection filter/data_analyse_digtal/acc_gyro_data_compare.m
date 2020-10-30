
% =========================================================================
%
%                  惯组数据处理
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月24日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.
%        2.
%        3.
%
%--------------------------------------------------------------------------
close all
clear all
%% 读取数据
N = 10000;
filepath = 'data/0916_1337_x/';
tmp3 = textread([filepath,'fmctrl_data_1337.txt']);
if length(tmp3)>N
    tmp3 = tmp3(1:N,:);
end
filepath = 'data/0916_1304_s/';
tmp4 = textread([filepath,'fmctrl_data_1337.txt']);
if length(tmp4)>N
    tmp4 = tmp4(1:N,:);
end
%% 惯组单位换算
gyroll_1 = tmp3(:,1)/3276.8/180*pi/1.31;
gyroll_2 = tmp4(:,1)/3276.8/180*pi/1.31;
TBS1 = tmp3(:,16)/1e5;
TBS2 = tmp4(:,16)/1e5;
figure;plot(gyroll_1);hold on;plot(gyroll_2);
figure;plot(gyroll_1 - gyroll_2);

figure;plot(TBS1);hold on;plot(TBS2);%%为什么速度能保持一致呢？








