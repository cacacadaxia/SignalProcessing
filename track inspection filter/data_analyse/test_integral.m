

% =========================================================================
%
%                  测试
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.差分的计算与结果对比
%        2.
%        3. 
%  进一步想法：
%           1.可以通过仿真的方法，对比分析误差量
%           2.怎么样去仿真呢？这是一个问题，是不是可以进行模态分析？
% 
% 10.14
%   1. 讨论二阶差分与积分直接的关系；积分实际上也是一种低通滤波的形式
% 
%--------------------------------------------------------------------------

%%
% d = textread('gplpe.txt');
% for i = 1:length(d)
%     if i == 1
%        x_dot = 0;
%        x_dot_2 = 0;
%        p_x_dot = 0;
%     else
%         x_dot = d(i) - d(i-1);
%         x_dot_2 = x_dot - p_x_dot;
%         p_x_dot = x_dot;
%         
%         x_dot_2_s(i) = x_dot_2;
%     end
% end
% 
% x_dot = 0;
% p_x_dot = 0;
% x = 0;
% p_x = 0;
% for i = 4:length(d)
%     
%     x_dot_2 = x_dot_2_s(i);
%     x_dot = x_dot + x_dot_2;
%     x = x + x_dot;
% 
%     
%     x_s(i) = x;
%     
% end
% 
% figure;
% plot([zeros(1,0), x_s]);
% hold on;
% plot(d)

%%
% clear all;
% close all;
% del_t = 0.01;
% t = 1:del_t:100;
% T = 10; 
% a = - 4*pi^2/T^2*sin( 2*pi/T*t );
% v = 2*pi/T;
% s = 0;
% for i = 1:length(a)
%    v(i+1) = v(i) + a(i)*del_t;
%    s(i+1) = s(i) + v(i)*del_t;
%   
% end
% figure;plot(s);

%%直接积分也会出问题的，那么怎么办？
%%因为省略了高阶的项，那也不应该啊，咋回事

%% 为什么用二阶差分的方法
clear all;
close all;
del_t = 0.01;
t = 1:del_t:100;
T = 10; 
a = - 4*pi^2/T^2*sin( 2*pi/T*t );
v = 2*pi/T;
s = 0;
s_dot = 0;
for i = 1:length(a)
   tmp = a(i) * del_t * del_t;
   s_dot(i+1) = s_dot(i) + tmp;
   s(i+1) = s(i) + s_dot(i);
end
figure;plot(s);

%% function 高通滤波器









