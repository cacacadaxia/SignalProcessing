


% =========================================================================
%
%                  测试非线性方法在位置估计中的应用
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.利用非线性优化的方法去观察超长波测量是不是可行的
%        2.
%       3. 
%--------------------------------------------------------------------------
clear all;
close all;
dt = 0.1;
t = 0:dt:100;
f = 1;
wave = sin(2*pi*f*t);
diff = 2*pi*f*cos(2*pi*f*t);

diff = diff + randn(1,length(diff))*5;
for i = 2:length(diff)
    z_out(i) = once_integral(diff(i-1)*dt , 0);
end
figure1 = figure('Color',[1 1 1]);
plot(z_out);hold on;plot(wave);
options = optimset('TolX',1e-6,'Algorithm','Levenberg-Marquardt',...
                   'Display','iter');
[a,resnorm]=lsqnonlin(@cau_E,z_out,[],[],options,diff);
plot(a,'k','LineWidth',1);
set(gca,'Fontname','Times New Roman','fontsize',14);


% 只是用这个优化的方法是不行的，只是幅值发生了变化
% 必须要对噪声进行建模才行，用卡尔曼滤波，但是又找不到观测量，所以也是不行的
% 那么怎么办呢？


%%
max_ = z_out(end);
for i = 1:length(z_out)
    z_out_2(i) = z_out(i) - max_/length(z_out)*i;
end
plot(z_out_2);
legend 加噪后 原始波形 优化波形 做平差;

%%
b = load('filter_10Hz.mat');
Len = (length(b.Num)-1)/2;

out = conv(b.Num,a);
out = out(Len:end-Len-1);
figure;plot(out,'LineWidth',1);hold on;plot(wave)
out = conv(b.Num,z_out);
out = out(Len:end-Len-1);
plot(out);
legend 1 2 3

%% 函数
function E = cau_E(x, diff )
E(1) = x(1) - 0;
E(2) = x(end) - 0;
dt = 0.1;
for i = 1:length(x)-1
    E(i+2) = x(i+1) - x(i) - dt * diff(i);
end
end
%% 两次积分
function y_k = once_integral( x_k , y_0 )
persistent  y_1;
if isempty(y_1)
    y_1 = y_0;
end

y_k = y_1 + x_k;
%%
y_1 = y_k;
end
    