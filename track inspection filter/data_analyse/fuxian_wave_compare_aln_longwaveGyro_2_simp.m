
% =========================================================================
%
%                  复现轨道检测的算法部分
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月10日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.首先复现最为简单的轨向部分
%        2.加上滤波器的部分，并进行对比。主要是轨向，结果很好
%        3. 修改长波计算的方法，用陀螺仪进行代替
%        4. 简化代码
% 
%--------------------------------------------------------------------------
close all;
clear all;
filepath = 'data/0916_1337_x/';
start_pos = 1; N = 1e4;
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;%%km表示

%%
tmp2 = textread([filepath , 'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(start_pos : start_pos+N-1 , :);
end
% gpxbr,hfcra,lfcrp
% sita_b = sita_b/3276.8/180*pi;
%%
sita_b = tmp2(:,1);
ay = fmctrl_data(:,5);      %%gpan

%% 轨距的对比

% ****************参数设定*********************
delay = 418;%%一直是这个数吗？
%%这两个参数是啥意思？这两个量是准确的
G_par = 3.8259e-14*141500.0;
ht = 3.90398e-05*141500.0*0.268;
tbs = fmctrl_data(:,end);
tbs_s = tbs/1e5;
% ***************step1 模型搭建*************************
rou_l = fmctrl_data(:,9);
rou_r = fmctrl_data(:,11);

for i = 3:length(rou_l)
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
    
end
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
end
camo = ay_Gz -  G_par .* sita_b .* tbs.^2 + ht * sita_b_dot2;

%% 陀螺仪取代
gyroYaw = fmctrl_data(:,6);
for i = 1:length(gyroYaw)
    gyroYaw_ = gyroYaw(i);
    tbs_ = tbs(i);
    yaw_Rz(i,1) = C(gyroYaw_ , tbs_);%%为什么需要这个？
end
sampleDistance = 0.25;
yawParameter = 2.0970;      %% 135751/9.8*6.0135*6.0135/4294.97*pi/180
%%yawParameter这个参数是怎么得到的？
gpyawReviseTemp1 = yaw_Rz * yawParameter * sampleDistance*( -1 );
gpyawRevise = gpyawReviseTemp1 + 0;
gpyawRevise = zeros(size(gpyawRevise));

camo = - gpyawRevise + ht * sita_b_dot2;

%%
camo = -camo;
camo = quzheng(camo);
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;

%%
figure;plot(amcol);hold on;plot(aln(:,3));
figure; plot( amcol - aln(:,3) ); title('amcol之间的对比');

%% 分析倍数
error = amcol ./ aln(:,3);

%%
plot_mag(ay,'ay');%%有问题

 %%
alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
sscal = 0.000825;
sbsci = 0.019802;
fscal = 0.1;
sbsc = 101.000;
Num = 768;
amcol_array = zeros(Num,1);
amcol_arraytmp = zeros(Num,1);
in1 = 533; in2 = 432; in4=382; in6=331; in7=230;%%这些参数都是固定的，不能变
in = 539;%%控制偏移值

%%正常的25m的滤波
for i = 1:length(amcol)           %%简单积分，肯定是不对的
    %%
   amcol_array(in) = amcol(i);
    %%由此可见就是输入出了问题，所以查找输入
    %%-----------------------------
    alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
    elupp = alu;
    elup = elup + elupp;
    elu = elu + elup;
    emco = - amcol_array(in4);
    als = als + amcol_array(in2) - amcol_array(in6);
    
    alss = alss + als + sbsc*emco;
    alsss = alsss + alss;
    xtemp = (sbsci*alsss - sscal*elu)*fscal;
    
    alsss_save(i) = alsss;
    elu_save(i) = elu;
    alss_save(i) = alss;
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
% figure;plot(yL,'k','LineWidth',0.5);hold on;plot(aln(:,4),'r','LineWidth',0.5);legend matlab gj
% figure;plot((yL - aln(:,4)));%%基本完全一致
title('左轨向的结果');
%%

%%  
% guixiang_l = wave_out(:,3);
% figure;plot(guixiang_l(268:end),'LineWidth',1);
% hold on;plot(aln(:,4));
% title('三点滤波的结果');

    
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

function plot_mag(signal_data , tit)
figure;
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

%% 在长波滤波之后进行积分？
function out = integrational(auat)
persistent  euas euasp  auatp eusa;
if isempty(euasp)
    euasp =  0;
    auatp = 0;
    eusa = 0;
    euas = 0;
end
gyroabase = 57;%%这是类似于窗长吗？
%% 计算
euas = euas + (auat + auatp) * 0.5;
eusa = (2 * gyroabase / (2 * gyroabase + 1) * eusa) + (euas + euasp) * 0.5;

%% 更新
out = eusa;
auatp = auat;
euasp = euas;

end



%% 无意义的代码
% %% 重新整理积分
% 
% for i = 1:length(amcol)           %%简单积分，肯定是不对的
%     %%
%    amcol_array(in) = amcol(i);
%     %%由此可见就是输入出了问题，所以查找输入
%     %%-----------------------------
%     
%     alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
%     elupp = alu;
%     elup = elup + elupp;
%     elu = elu + elup;
%     emco = - amcol_array(in4);
% 
%     als = als + amcol_array(in2) - amcol_array(in6);
%     
%     alss = alss + als + sbsc*emco;
%     alsss = alsss + alss;
%     xtemp = (alsss*sbsci - sscal*elu)*fscal;
%     yL(i,1) = xtemp;
%     
%     %% 更新数组
%     in1 = mod(in1,Num)+1;
%     in2 = mod(in2,Num)+1;
%     in4 = mod(in4,Num)+1;
%     in6 = mod(in6,Num)+1;
%     in7 = mod(in7,Num)+1;
%     in = mod(in,Num)+1;
%     %%
%     save(i,1) = alu;
%     save(i,2) = elu;
%     save(i,3) = als;
%     save(i,4) = alsss;
% end

