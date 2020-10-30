
% =========================================================================
%
%                  模拟系统的高低部分
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月27日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.轨向部分
%        2.
%        3.
%        4. 
% 
% 
%--------------------------------------------------------------------------

%%
% load_txt;
close all;
% clear all;
filepath = 'data/1026/';
start_pos = 1;
N = 1e4;
% load_txt;
x = 0:0.25:0.25*(N-1);
x = x/1000;

%% 分析数据
gpvmp = output_wave(:,1);%%km
gpvft = output_wave(:,2);%%
gplpe3 = output_wave(:,43);
gprpe3 = output_wave(:,44);
gptbs = fmctrl_data(:,20);

gpan = fmctrl_data(:,13);
gpgl = fmctrl_data(:,14);
gpgr = fmctrl_data(:,15);

sita_b = aln(:,6);


%%
TBS = gptbs;
ay = gpan;
rou_l = gpgl;
rou_r = gpgr;

%%
% ****************参数设定*********************
delay = 418; %%一直是这个数吗？不一定，每个测量项目对应的参数都不一样

%% 确定这两个参数的含义
%%都是标定出来的量
G_par = 3.8259e-14*141500.0;
ht = 3.90398e-05*141500.0*0.268;
% ***************step1 模型搭建*************************
rou_l_dot2(1,1) = rou_l(1);
rou_r_dot2(1,1) = rou_r(1);
rou_l_dot2(2,1) = rou_l(2) - 2*rou_l(1);
rou_r_dot2(2,1) = rou_r(2) - 2*rou_r(1);
for i = 3:length(rou_l)
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end
% ---------------------- 经过滤波器 --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , TBS(i));
end
%% 求amcol与amcor
sita_b_dot2(2,1) = sita_b(2) - sita_b(1) - sita_b(1);
sita_b_dot2(1,1) = sita_b(1);
for i = 3:length( sita_b )
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = - ( ay_Gz - G_par .* sita_b .* TBS.^2 + ht * sita_b_dot2 );
camo = quzheng( camo );%%这里说明了取整带来了问题
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;

%% 对比
amcol_ref = aln(:,3);
figure;plot( amcol - amcol_ref );
yL = shortwave_filter( amcol ); %%25m以下短波
figure1 = figure('Color',[1 1 1]);plot( yL(9:end),'r','LineWidth',0.5);hold on;plot(aln(:,4),'-');legend matlab gj;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('采样点 /0.25m');ylabel('int /32768/10 inch')



%% 对比结果
% yL_2 = longwave_filter(amcol,281,71,281,491);%%都是N的值

%%
% plot_mag(yL,'25m');
% plot_mag(yL_2,'单边幅值谱,截止波长30m和70m对比','hold');
% legend 30m 70m;





%% 短波滤波器
function out = shortwave_filter(in)
amcol = in;
%% 25m长波滤波加积分
alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
sscal = 0.000825;
sbsci = 0.019802;
fscal = 0.1;
sbsc = 101.000;
% 数组设定
Num = 768;
amcol_array = zeros(Num,1);
amcol_arraytmp = zeros(Num,1);
in1 = 533; in2 = 432;in4=382;in6=331;in7=230;
in = 539;

for i = 1:length(amcol) 
    amcol_array(in) = amcol(i);
    alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
    elupp = alu;
    elup = elup + elupp;
    elu = elu + elup;
    emco = - amcol_array(in4);
    als = als + amcol_array(in2) - amcol_array(in6);
    alss = alss + als;
    alss = alss + sbsc*emco;
    alsss = alsss + alss;
    xtemp = (alsss*sbsci - sscal*elu)*fscal;
    yL(i,1) = xtemp;
    
    %% 更新数组
    in1 = mod(in1,Num)+1;
    in2 = mod(in2,Num)+1;
    in4 = mod(in4,Num)+1;
    in6 = mod(in6,Num)+1;
    in7 = mod(in7,Num)+1;
    in = mod(in,Num)+1;
    %%
    save(i,1) = alu;
    save(i,2) = elu;
    save(i,3) = als;
    save(i,4) = alsss;
end
%%
out = yL;
end
%% 前端滤波器
function out = C(x_k,tbs)
omega1 = 0.76;

persistent x;
if isempty(x)
    x = zeros(2,1);
end
%%
x(2) = x_k;
y = x(2)-x(1) + tbs/2^18*(x(2)+x(1));
y = y/omega1;
%%
x(1) = x(2);
out = y;
end

function out = F(x,tbs)
%% Fs的标准写法
%% 滤波器设定
%% 对于x[3]，其与同理
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2
%% 这个滤波器存在延时，几十不等
persistent y;
if isempty(y)
    y = zeros(3,1);
end
y(3) = ( y(2)*(2*2^28+2^14*tbs) - y(1)*2^28+tbs^2*x  )/(2^28 + 2^14*tbs + tbs^2);
%% 更新temp量
y(1) = y(2);
y(2) = y(3);
out = y(3);
end

function out = R(x_k,tbs)
wd = 0.001;

persistent y x;
if isempty(y)
    y = zeros(2,1);
    x = zeros(2,1);
end
x(2) = x_k;
x_dot = x(2)-x(1);
y(2) = (1 - wd) * (x_dot + y(1));

%% 更新temp量
x(1) = x(2);
y(1) = y(2);
out = y(2);

end
function out = G(x_k,tbs_k)
persistent x y1 tbs;
if isempty(x)
%     y = zeros(2,1);
    x = zeros(3,1);
    y1 = zeros(2,1);
    tbs = zeros(2,1);
end
%% 更新k时刻
x(3) = x_k;
tbs(2) = tbs_k;

%% 离散滤波
x_dot = x(3) - x(2);
x_dot_2 = x(3) - 2*x(2) + x(1);
y1(2) = (x(3) + x(2)) *tbs(2)/2^15 + x_dot;
y = (y1(2)+y1(1))* (tbs(2) + tbs(1))/2^16 + x_dot_2;


%% 更新temp量
y1(1) = y1(2);
x(1) = x(2);
x(2) = x(3);
tbs(1) = tbs(2);
out = y;

end

function plot_mag(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure;
end
fs = 4;     %% 0.25m为一个采样间隔
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title(tit);
grid on;
end
function out = quzheng(in)
out = zeros(2,1);
for i = 1:length(in)
    if in(i)>=0
        out(i) = floor(in(i));
    elseif in(i)<0
        out(i) = floor(in(i)) + 1;
    end
end
end



