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
%        2.
%        3. 
%--------------------------------------------------------------------------

close all;
clear all;
addpath('func');
start_pos = 5000;
N = 10000;
filepath = 'data/0916_1337_x/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
% gpxbr,hfcra,lfcrp
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
tmp2 = tmp2(start_pos:start_pos+N-1,:);
% 
tmp3 = textread([filepath , 'gppro.txt']);
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
    az_Gz(i,1) = G(gppl(i)  , TBS(i));
end
az_Gz = quzheng(az_Gz);
amarm = sita_b_dot2.*ararm_par;
pgrav = sita_b.^2.*TBS.^2*pgscale;
plmct = -(az_Gz + pgrav + amarm);
prmct = -(az_Gz + pgrav - amarm);
% plmct = quzheng(plmct);
% prmct = quzheng(prmct);
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

%% 陀螺仪替代
%%一些参量记得单位就可以了，但是这里为什么不对
pitchParameter = 2.0970;
sampleDistance = 0.25;
compf = -1;
pitch = gypitch;
for i = 1:length(pitch)
    pitch(i) = C(pitch(i) , TBS(i));
end
temp = sampleDistance * pitch * compf * pitchParameter;
plmct = - ( temp + amarm );
prmct = - ( temp - amarm );

% 其余的都是一样的
% pmcol = plmct + gpvlo_dot2;
% pmcor = prmct + gpvro_dot2;

%% 长波滤波并积分
zL = longwave_filter(pmcol,281,71,281,491);
zR = longwave_filter(pmcor,281,71,281,491);

%% 结果对比
pmcol_ref = tmp3(:,1);
pmcol_ref70m = tmp3(:,2);
figure;plot(pmcol_ref,'k','LineWidth',0.5);hold on;plot(pmcol,'r--','LineWidth',0.5)
figure;plot(pmcol_ref - pmcol);

%% 为什么会有一个系数呢？感觉很奇怪
% ratio = 1;figure;plot(ratio.*zL);hold on;plot(pmcol_ref70m(177:end));legend matlab gj

%%
% plot_mag(pmcol_ref,'25m')
% plot_mag(pmcol_ref70m,'70m','hold')
% legend 1 2



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