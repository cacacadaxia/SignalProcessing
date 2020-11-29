% =========================================================================
%
%                  用陀螺仪取代加速度计测量长波
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月26日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.高低借鉴轨向的做法，对过程进行简化
%        2.高低基本长波的数值没啥问题
% 10月
%        3. 新增平差试验
%--------------------------------------------------------------------------

%%
close all;
clear all;
addpath('func');
start_pos = 1;
N = 2e4;
filepath = 'data/0916_1337_x/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
% gpxbr,hfcra,lfcrp
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
tmp2 = tmp2(start_pos:start_pos+N-1,:);
% pmcol\gplpe\gplpe1(70m)
tmp3 = textread([filepath , 'LongWaveResultForGppro_L.txt']);
tmp3 = tmp3(start_pos:start_pos+N-1,:);
% 
tmp4 = textread([filepath , 'fmctrl_data_1337.txt']);
tmp4 = tmp4(start_pos:start_pos+N-1,:);

%% 读取高低相关的数据
gypitch = fmctrl_data(:,2);
gpvlo = fmctrl_data(:, 8);
gpvro = fmctrl_data(:,10);
gppl  = fmctrl_data(:,7);   %%z轴向加速度
TBS   = fmctrl_data(:,end);
sita_b = tmp2(:,1);

%% 信号处理，获得差分值，前两个先不管
% ------------1. 参数设定--------------------------
% 增益值的作用是什么？这个让人很困惑？
% 接下来要怎么办？
canMagpar = -159900;
ararm_par = 2.9280e-5 * abs(canMagpar);
pgscale = 1.9646e-19 * abs(canMagpar);
% ------------2. 计算差分------------------------
for i = 3:length(gpvlo)
    tbs_ = TBS(i);
    gpvlo(i) = gpvlo(i)*2;
    gpvro(i) = gpvro(i)*2;
    gpvlo_dot2(i,1) = gpvlo(i) - 2*gpvlo(i-1) + gpvlo(i-2);
    gpvro_dot2(i,1) = gpvro(i) - 2*gpvro(i-1) + gpvro(i-2);
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
    
end
for i = 1:length(gppl)
    az_Gz(i,1) = G( gppl(i)  , TBS(i) );
end
az_Gz = quzheng( az_Gz );
amarm = sita_b_dot2.*ararm_par;
pgrav = sita_b.^2.*TBS.^2*pgscale;
% plmct = -(az_Gz + pgrav + amarm);
% prmct = -(az_Gz + pgrav - amarm);
disp(mean(abs(amarm)))
disp(mean(abs(az_Gz)))
plmct = -(az_Gz - pgrav + amarm);%%修改之后也不变化
prmct = -(az_Gz - pgrav - amarm);

% plmct = quzheng(plmct);
% prmct = quzheng(prmct);
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

pmcol = pmcol(300:end);
pmcor = pmcor(300:end);

%% 陀螺仪替代
pitchParameter = 2.0970;
sampleDistance = 0.25;
compf = -1;
pitch = gypitch;
for i = 1:length(pitch)
    pitch(i) = C(pitch(i) , TBS(i));%%Tbs包含在了C中
end
temp = sampleDistance * pitch * compf * pitchParameter;%%这个竟然是对的
% temp = 0;
plmct = - ( temp + amarm );
prmct = - ( temp - amarm );

% 其余的都是一样的
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

%% 短波滤波
zL_ref30m = tmp3(:,2);
zL_ref70m = tmp3(:,3);
zL_30m = shortwave_filter(pmcol);

%% 长波滤波并积分
zL_70m = longwave_filter(pmcol,281,71,281,491);
zR_70m = longwave_filter(pmcor,281,71,281,491);
zL_250m = longwave_filter(pmcol,1123,307,1123,1953);
zR_250m = longwave_filter(pmcor,1123,307,1123,1953);

% figure1 = figure('Color',[1 1 1]);plot(zL_250m);hold on;plot(zR_250m);set(gca,'Fontname','Times New Roman','fontsize',16);

%% fdatool 设计高通滤波器
b = load('filter1.mat');%%70m的长波
% b = load('filter_250m.mat');%%250m长波
% b = load('filter_300m.mat');%%300m长波
% b = load('filter_400m.mat');%%400m长波
% b = load('filter_1000m.mat');%%1000m长波


tmpN = (length(b.Num)-1)/2;
pmcol_70m_tmp = conv(b.Num,pmcol);
pmcol_70m_tmp(1:tmpN) = [];pmcol_70m_tmp(end-tmpN+1:end) = [];%%去除暂态点之后
z_dot = -9.6;   %%70m长波参数
% z_dot = -19.7;%%250m
% z_dot = -4.7;%%300m
% z_dot = -52.7;%%400m
% z_dot = -261;%%250m 加速度计
% z_dot = -132;
zL_70m_fdatool = 0;
for i = 1:length(pmcol_70m_tmp)
    z_dot = z_dot + pmcol_70m_tmp(i);
    zL_70m_fdatool = zL_70m_fdatool + z_dot;
    z_dot_save(i) = z_dot;
    zL_70m_fdatool_save(i) = zL_70m_fdatool;
end

% figure1 = figure('Color',[1 1 1]);plot( -zL_70m_fdatool_save/5);hold on;plot(zL_70m(247:end));
% legend fdatool设计的滤波器 窗函数滤波器(70m);set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('里程 /0.25m');ylabel('高低 / (32768/10 inch)');
figure1 = figure('Color',[1 1 1]);plot( -zL_70m_fdatool_save/5/129.01);hold on;plot(zL_70m(247:end)/129.01);
legend fdatool设计的滤波器 窗函数滤波器(70m);set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('里程 /0.25m');ylabel('高低 / mm');
%%这个滤波器是存在漂移的，用窗函数不仅计算量小，而且不产生信号的畸变和漂移，很神奇
%%这里的幅值产生变化就很奇怪，为什么会对幅值产生影响呢？
%%这里也是对这个二阶差分进行滤波
%%这里需要注意的地儿就是：zL_70m(247:end)才是参考的点的位置，这个很重要

%% 陀螺方案与fdatool对比
% figure;plot(zL_250m(1000:end));hold on;plot(-zL_70m_fdatool_save/5);legend 1 2
% figure;plot(zL_ref30m/129.01);

%%
% figure1 = figure('Color',[1 1 1]);plot(zL_30m);hold on;plot(zL_70m(89:end));legend 陀螺仪取代方法30m 陀螺仪取代方法70m;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('里程 /0.25m');ylabel('高低 / (32768/10 inch)');grid on;

%% 结果对比
pmcol_ref = tmp3(:,1);
zL_ref30m = tmp3(:,2);
zL_ref70m = tmp3(:,3);

% figure;plot(pmcol_ref,'k','LineWidth',0.5);hold on;plot(pmcol,'r--','LineWidth',0.5)
% figure;plot(pmcol_ref - pmcol);

%% 为什么会有一个系数呢？感觉很奇怪，这个系数是做什么的？
%%因为这里左右数据输出输出反了
% ratio = 1; figure1 = figure('Color',[1 1 1]);plot(ratio.*zL_30m);hold on;plot(zL_ref30m(1:end));legend 陀螺仪取代方法 原方法;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('里程 /0.25m');ylabel('高低 / (32768/10 inch)');title('30m')
% ratio = 1; figure1 = figure('Color',[1 1 1]);plot(ratio.*zL_70m);hold on;plot(zL_ref70m(177:end));legend 陀螺仪取代方法 原方法;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('里程 /0.25m');ylabel('高低 / (32768/10 inch)');title('70m')

%% 查看频谱，怎么感觉有点问题
plot_mag(yL_ref30m,'25m')
plot_mag(yL_ref70m,'70m','hold')
legend 1 2

%% 300m的桥里面有1200个点(假设这里是一个300m的桥，但是实际上这样考虑并不是很合理，因为
% 并不能保证一开始的就是桥墩的位置，只能假设)
% % 300m-->1200
% start = 4000;
% z_dot = 0;
% z = 0;
% len = 300*4;
% pmcolBridge = pmcol(start+1:start+len);
% pmcolBridge = pmcolBridge - mean(pmcolBridge);
% for i = 1:len
%     z_dot = z_dot + pmcolBridge(i);
%     z = z + z_dot;
%     z_save(i) = z;
%     z_dot_save2(i) = z_dot;
% end
% 
% figure1 = figure('Color',[1 1 1]);grid on;
% subplot(2,1,1);plot(z_dot_save2);title('一阶差分');
% xlabel('采样点 /0.25m');ylabel('int')
% set(gca,'Fontname','Times New Roman','fontsize',14);
% % subplot(2,1,2);
% figure;plot(z_save);title('高低不平顺');xlabel('采样点 /0.25m');ylabel('int');set(gca,'Fontname','Times New Roman','fontsize',14);
% % figure1 = figure('Color',[1 1 1]);plot(zL_ref70m(start:start+len-1));hold on;plot(zL_ref30m(start:start+len-1));set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% % xlabel('采样点 /0.25m');ylabel('数值 /32768/10 inch');legend 70m 30m
% 
% % max_ = z_save(end);
% % del = max_/len;
% % for i = 1:len
% %     z_save(i) = z_save(i) - del*i;
% % end
% % figure1 = figure('Color',[1 1 1]);
% % plot(z_save);hold on;plot(yL_ref70m(start:start+len-1));
% % xlabel('采样点 /0.25m');ylabel('数值 /32768/10 inch');set(gca,'Fontname','Times New Roman','fontsize',16);





% plot_mag()
%% 画图的函数
function plot_for_paper()
% figure1 = figure('Color',[1 1 1]);plot( -zL_70m_fdatool_save/5/129.01);hold on;plot(zL_70m(247:end)/129.01);
% legend fdatool设计的滤波器 窗函数滤波器(70m);set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('里程 /0.25m');ylabel('高低 / mm');
end
%% 函数
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

% out = y/tbs(1)/tbs(2);
out = y;
tbs(1) = tbs(2);


end

function out = find_index(in)
%% 用在滤波器的队列中
Num = 2048;
if mod(in,Num) == 0
    out = Num;
else
    out = mod(in,Num);
end
end

function out = point3filter(in)
%% 简单的三点滤波
persistent x;
if isempty(x)
    x = zeros(3,1);
end
x(3) = in;
out = sum(x)/3;
%%
x(1) = x(2);
x(2) = x(3);

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


function out = B1(x,tbs)
persistent y;
if isempty(y)
    y = zeros(2,1);
end
Omega1 = 10^5/2^17;
y(2) = ( y(1) *2^17 + tbs*x )/(2^17 + tbs);

%%
y(1) = y(2);
out = y(2);
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
