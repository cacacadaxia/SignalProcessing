

clear all;
close all;



%% bias 
v = 10/3.6;
dt = 0.25/v; %%任意都行
BS = 0.01;    %%o/√h-->6对应√(3600/100)
sigma = BS*sqrt(dt)/sqrt(100);%%400s
%%加噪怎么加？这个需要考虑一下
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
% dt = 0.25/v; %%任意都行
% sigma = 0.2;
% %%加噪怎么加？这个需要考虑一下
% for j = 1:10000
%     ARW = randn*sigma;
%     ARW_save(j) = ARW;
% end
% std(ARW_save)