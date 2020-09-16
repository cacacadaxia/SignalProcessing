


 % =========================================================================
%
%                  my project compare
%
% =========================================================================
%
%��(C)2003-2010 CARC
%   �汾��V1.0
%   ���ڣ�2020�� 8��20��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ݸ�ʽ��text��timestamp, x, y, z
%        2.
%        3.
%        4.
%        5.
%        6.
%--------------------------------------------------------------------------

clear;
close all;

data2 = textread('groundtruth.txt');
data1 = textread('text.txt');


% ƥ��ʱ��
t_true = data2(:,1);
data_true_ag=[];
indexlist = [];
data1(:,1) = data1(:,1);
for i=1:length(data1)
    t1 = data1(i,1);
    [~,index] = min(abs(t_true-t1));
    data_true_ag = [data_true_ag;
        t_true(index),data2(index,6),data2(index,7),data2(index,8)];
    indexlist(i) = index;
end


figure;
subplot(3,1,1);plot(data_true_ag(:,2)-data_true_ag(1,2));hold on;plot(-data1(:,2));legend ground ag
subplot(3,1,2);plot(data_true_ag(:,3)-data_true_ag(1,3));hold on;plot(data1(:,3));
subplot(3,1,3);plot(data_true_ag(:,4)-data_true_ag(1,4));hold on;plot(data1(:,4));

