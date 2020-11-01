
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
%  功能： 1.直接进行积分处理，那无疑可能会存在问题
%        2.高低这样怎么弄？因为与文档上不一致
%        3.桥梁的数据处理
%        4. 
% 
% 
%--------------------------------------------------------------------------

%%
% load_txt;
close all;
clear all;
filepath = 'data/1026/';
start_pos = 2e4;
N = 4e4;
load_txt;
x = 0:0.25:0.25*(N-1);
x = x/1000;

gpvmp = output_wave(:,1);%%km
gpvft = output_wave(:,2);%%
gplpe3 = output_wave(:,43);
gprpe3 = output_wave(:,44);
gppl = fmctrl_data(:,11);
gppr = fmctrl_data(:,12);
TBS = fmctrl_data(:,20);
gplae = output_wave(:,5);

%% 一整个桥的
start_mileage = 398.2;  end_mileage = 399.5;
index_start = find_index(gpvmp,gpvft,start_mileage);
index_end = find_index(gpvmp,gpvft,end_mileage);
x = gpvmp(index_start:index_end) + gpvft(index_start:index_end)/4000;

% figure1 = figure('Color',[1 1 1]);plot(x,gplpe3(index_start:index_end));xlabel('里程 /1km');ylabel('长波 /3276.8/1 inch')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% figure1 = figure('Color',[1 1 1]);plot(x,gplae(index_start:index_end));xlabel('里程 /1km');ylabel('长波 /3276.8/1 inch')
% set(gca,'Fontname','Times New Roman','fontsize',16);


%% 其中的一段
bridge0 = 398.43475;
bridge1 = 398.609;
bridge2 = 398.757;
bridge3 = 398.948;
bridge4 = 399.08225;
bridge5 = 399.2105;
index_start = find_index(gpvmp,gpvft,bridge1);
index_end = find_index(gpvmp,gpvft,bridge3);

% index_start = 1;
% index_end = 40000;

%% 短弦中支矩经过短波滤波器然后对比
x = gpvmp(index_start:index_end) + gpvft(index_start:index_end)/4000;
zL30mref = output_wave(index_start:index_end,3);
zL120mref = output_wave(index_start:index_end,43);
pmcol = gppro_data(index_start:index_end,1);
zL30m = shortwave_filter(pmcol);%%经过短波滤波器
    % % figure1 = figure('Color',[1 1 1]);plot(zL30mref(347:end));hold on;plot(zL30m(1:end));
    % % xlabel('里程 /0.25m');ylabel('长波 /3276.8/1 inch')
    % % set(gca,'Fontname','Times New Roman','fontsize',16);
    % % legend gj matlab 
    % % grid on

%% 70m
% 897	245	897	1561	
% zL200m = longwave_filter(pmcol,897,245,897,1561);%%200m
zL200m = longwave_filter(pmcol,481,121,481,841);%%120m
% zL200m = longwave_filter(pmcol,2245,615,2245,3909);%%350m
% figure1 = figure('Color',[1 1 1]);plot(zL120mref(347:end));hold on;plot(zL200m(1:end));%%30m
figure1 = figure('Color',[1 1 1]);plot(zL120mref(83:end));hold on;plot(zL200m(1:end)*1.21);%%120m
% figure1 = figure('Color',[1 1 1]);plot(zL120mref(1:end)/129.01);hold on;plot(zL200m(287:end)/129.01);%%200m
xlabel('里程 /0.25m');ylabel('长波 /mm');set(gca,'Fontname','Times New Roman','fontsize',16);legend gj matlab;grid on;
% figure;plot(zL120mref(83:end));hold on;plot(zL200m(1:end))
% plot_mag(zL120mref,'120m');
% plot_mag(zL200m,'zL200m','hold');


%% 直接对短线二阶中支矩进行积分
% zdot = 54.98484;
% zdot = 1.06931e2;
zdot = 2.191974e2;
z = 0;
for i = 1:length(pmcol)
    zdot = zdot + pmcol(i);
    z = z + zdot;
    zdot_save(i) = zdot;
    z_save(i) = z;
end
figure1 = figure('Color',[1 1 1]);plot(x,z_save/129);xlabel('里程 /km');ylabel('长波 /mm')
set(gca,'Fontname','Times New Roman','fontsize',16);grid on



%% 加速度计经过Gz然后再进行滤波

for i = 1:length(gppl)
    gpplGz(i) = G(gppl(i),TBS(i));
end
x = gpvmp + gpvft/4000;
% result = longwave_filter(gpplGz,539,147,539,937);%%120m
result = longwave_filter(gpplGz,2245,615,2245,3909);%%120m
result = result(1:end-1);
figure1 = figure('Color',[1 1 1]);plot(x , result);
set(gca,'Fontname','Times New Roman','fontsize',16);grid on;
xlabel('里程 /km');ylabel('加速度计积分')
hold on;
b0 = find_index(gpvmp,gpvft,bridge0);
b1 = find_index(gpvmp,gpvft,bridge1);
b2 = find_index(gpvmp,gpvft,bridge2);
b3 = find_index(gpvmp,gpvft,bridge3);
b4 = find_index(gpvmp,gpvft,bridge4);
b5 = find_index(gpvmp,gpvft,bridge5);
plot(gpvmp(b0) + gpvft(b0)/4000,result(b0),'ro','LineWidth',2);
plot(gpvmp(b1) + gpvft(b1)/4000,result(b1),'ro','LineWidth',2);
plot(gpvmp(b2) + gpvft(b2)/4000,result(b2),'ro','LineWidth',2);
plot(gpvmp(b3) + gpvft(b3)/4000,result(b3),'ro','LineWidth',2);
plot(gpvmp(b4) + gpvft(b4)/4000,result(b4),'ro','LineWidth',2);
plot(gpvmp(b5) + gpvft(b5)/4000,result(b5),'ro','LineWidth',2);


%% 画了一整个桥
% 
% figure1 = figure('Color',[1 1 1]);plot(x , gplpe3);
% set(gca,'Fontname','Times New Roman','fontsize',16);grid on;
% xlabel('里程 /km');ylabel('长波 /mm')
% hold on;
% b0 = find_index(gpvmp,gpvft,bridge0);
% b1 = find_index(gpvmp,gpvft,bridge1);
% b2 = find_index(gpvmp,gpvft,bridge2);
% b3 = find_index(gpvmp,gpvft,bridge3);
% b4 = find_index(gpvmp,gpvft,bridge4);
% b5 = find_index(gpvmp,gpvft,bridge5);
% plot(gpvmp(b0) + gpvft(b0)/4000,gplpe3(b0),'ro','LineWidth',2);
% plot(gpvmp(b1) + gpvft(b1)/4000,gplpe3(b1),'ro','LineWidth',2);
% plot(gpvmp(b2) + gpvft(b2)/4000,gplpe3(b2),'ro','LineWidth',2);
% plot(gpvmp(b3) + gpvft(b3)/4000,gplpe3(b3),'ro','LineWidth',2);
% plot(gpvmp(b4) + gpvft(b4)/4000,gplpe3(b4),'ro','LineWidth',2);
% plot(gpvmp(b5) + gpvft(b5)/4000,gplpe3(b5),'ro','LineWidth',2);
%%

%% 加速度计的积分，这里有点问题
gppl_ = gppl(index_start:index_end)/13853;
tbs = TBS(index_start:index_end)/1e5;
z = 0;
zdot = 3.008759058979370e-04;
zdot = 2.352158858659502e-04;
for i = 1:length(tbs)
    del = gppl_(i)*tbs(i)*tbs(i);
    zdot = zdot + del;
    z = z + zdot;
    zdot_save(i) = zdot;
    z_save(i) = z;
end
% figure1 = figure('Color',[1 1 1]);plot(x,z_save*1e3);xlabel('里程 /km');ylabel('长波 /mm');
% set(gca,'Fontname','Times New Roman','fontsize',16);grid on;

%%


%% 函数
%% 输入里程，找到index
function index = find_index(gpvmp,gpvft,pos)
par1 = floor(pos);
par2 = pos - floor(pos);
in3 = floor(par2*4000);
in1 = find(gpvmp == par1);
in2 = gpvft(in1(1));
index = in1(1) + in3 - in2;
end
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
    figure1 = figure('Color',[1 1 1]);
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



