



% =========================================================================
%
%                  ��֤�˲���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.��֤sita_b��׼ȷ��
%        2.
%        3. 
%--------------------------------------------------------------------------



clear all;
N = 10000;
%% ��ȡ����
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
%% ��������
figure;plot(tmp(:,3) + tmp2(:,3));
hold on ;
plot(gpxbr);
legend 1 2

figure;plot(tmp3(:,2));hold on;
plot(tmp(:,3))
legend high low



