
% =========================================================================
%
%                  ����IMU������ģ��
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------


clear all;
close all;



%% bias 
% v = 10/3.6;
% dt = 0.25/v; %%���ⶼ��
% BS = 0.01;    %%o/��h-->6��Ӧ��(3600/100)
% sigma = BS*sqrt(dt)/sqrt(100);%%400s%%����������ߵİ�����
% %%������ô�ӣ������Ҫ����һ��
% for j = 1:10000
%     ARW = 0;
%     for i = 0:dt:100
%         ARW = ARW + randn*sigma;
%     end
%     ARW_save(j) = ARW;
% end
% std(ARW_save)*sqrt(dt) - sigma*sqrt(100)


%%
v = 1/3.6;
dt = 0.25/v; %%���ⶼ��
BS = 0.01;    %%o/��h-->6��Ӧ��(3600/100)
sigma = BS*sqrt(dt);%%400s%%����������ߵİ�����
%%������ô�ӣ������Ҫ����һ��
t_end = 1000;
for j = 1:10000
    ARW = 0;
    for i = 0:dt:t_end
        ARW = ARW + randn*sigma;
    end
    ARW_save(j) = ARW;
end
std(ARW_save)
BS*sqrt(t_end)