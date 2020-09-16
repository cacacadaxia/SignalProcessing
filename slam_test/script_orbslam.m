

 % =========================================================================
%
%                  slam对比
%
% =========================================================================
%
%　(C)2003-2010 CARC
%   版本：V1.0
%   日期：2020年 8月13日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.对比orb-slam2的性能
%        2.
%        3.
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------


clear;
close all;

data2 = textread('groundtruth.txt');
data1 = textread('KeyFrameTrajectory.txt');
data3 = textread('CameraTrajectory.txt');


% 匹配
t_true = data2(:,1);
data_true_ag=[];
indexlist = [];
% data1(:,1) = data1(:,1) + 1305031102;
for i=1:length(data1)
    t1 = data1(i,1);
    [~,index] = min(abs(t_true-t1));
    data_true_ag = [data_true_ag;
        t_true(index),data2(index,6),data2(index,7),data2(index,8)];
    indexlist(i) = index;
end
% figure;plot(data_true_ag(:,2)-data_true_ag(1,2));hold on;plot(-data1(:,6));legend true ag
% figure;plot(data_true_ag(:,3)-data_true_ag(1,3));hold on;plot(data1(:,7));legend true ag
% figure;plot(data_true_ag(:,4)-data_true_ag(1,4));hold on;plot(data1(:,8)-1);legend true ag

%%
% 匹配
t_true = data2(:,1);
data_true_ag=[];
indexlist = [];
% data1(:,1) = data1(:,1) + 1305031102;
for i=1:length(data3)
    t1 = data3(i,1);
    [~,index] = min(abs(t_true-t1));
    data_true_ag = [data_true_ag;
        t_true(index),data2(index,6),data2(index,7),data2(index,8)];
    indexlist(i) = index;
end
figure;
subplot(3,1,1);plot(data_true_ag(:,2)-data_true_ag(1,2));hold on;plot(-data3(:,6));legend true ag
subplot(3,1,2);plot(data_true_ag(:,3)-data_true_ag(1,3));hold on;plot(data3(:,7));legend true ag
subplot(3,1,3);plot(data_true_ag(:,4)-data_true_ag(1,4));hold on;plot(data3(:,8)-1);legend true ag


