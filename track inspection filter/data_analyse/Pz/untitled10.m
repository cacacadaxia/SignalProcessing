



% =========================================================================
%
%                  验证滤波器
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.验证sita_b的准确性
%        2.
%        3. 
%--------------------------------------------------------------------------



clear all;
N = 10000;
%% 读取数据
tmp = textread('Hz_filter_inc.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end

tmp2 = textread('Ps3_filter_wx.txt');
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end

tmp3 = textread('tmp_zhongjian_1337.txt');
if length(tmp3)>N
    tmp3 = tmp3(1:N,:);
end
gpxbr = tmp3(:,1);
%% 分析数据
figure;plot(tmp(:,3) + tmp2(:,3));
hold on ;
plot(gpxbr);
legend 1 2

figure;plot(tmp3(:,2));hold on;
plot(tmp(:,3))
legend high low



