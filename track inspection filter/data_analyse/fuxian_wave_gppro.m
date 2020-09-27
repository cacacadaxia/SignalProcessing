% =========================================================================
%
%                  高低检测项目的复现
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
start_pos = 1;
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

%% 读取高低相关的数据
gpvlo = fmctrl_data(:, 8);
gpvro = fmctrl_data(:,10);
gppl  = fmctrl_data(:,7);   %%z轴向加速度
TBS   = fmctrl_data(:,end);
sita_b = tmp2(:,1);
%% 信号处理，获得差分值，前两个先不管
% ------------1. 参数设定--------------------------
% 增益值的作用是什么？
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
    az_Gz(i,1) = G(gppl(i) , tbs_);%%加速度信号处理
end
amarm = sita_b_dot2.*ararm_par;
pgrav = sita_b.^2.*TBS.^2*pgscale;
plmct = -(az_Gz + pgrav + amarm);
prmct = -(az_Gz + pgrav - amarm);
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

%% 长波滤波并积分
zL = longwave_filter(pmcol);


%% 结果对比
pmcol_ref = tmp3(:,1);
pmcol_ref70m = tmp3(:,2);
figure;plot(pmcol_ref,'k','LineWidth',0.5);hold on;plot(pmcol,'r--','LineWidth',0.5)
figure;plot(pmcol_ref - pmcol);
%% 结果对比2
figure;plot(zL);hold on;plot(pmcol_ref70m);

%% 函数
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




