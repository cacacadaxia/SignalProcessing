

clear all;
close all;



%% bias 
v = 10/3.6;
dt = 0.25/v; %%���ⶼ��
BS = 0.01;    %%o/��h-->6��Ӧ��(3600/100)
sigma = BS*sqrt(dt)/sqrt(100);%%400s
%%������ô�ӣ������Ҫ����һ��
for j = 1:10000
    ARW = 0;
    for i = 0:dt:100
        ARW = ARW + randn*sigma;
    end
    ARW_save(j) = ARW;
end
std(ARW_save)



%% 
% v = 100/3.6;
% dt = 0.25/v; %%���ⶼ��
% sigma = 0.2;
% %%������ô�ӣ������Ҫ����һ��
% for j = 1:10000
%     ARW = randn*sigma;
%     ARW_save(j) = ARW;
% end
% std(ARW_save)