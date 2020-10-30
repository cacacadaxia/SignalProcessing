
% =========================================================================
%
%                  复现轨道检测的算法部分
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.首先复现最为简单的轨向部分
%        2.加上滤波器的部分，并进行对比。主要是轨向，结果很好
%        3. 修改积分的方法
%       4. 用fdatool设计长波滤波器，观察效果
% fuxian_wave_compare_aln_longwave为模版，然后更换其积分的部分
%       5. 观察频谱，应该算是正确的方法
% 
%--------------------------------------------------------------------------

close all;
clear all;
filepath = 'data/0916_1337_x/';
start_pos = 1;
N = 19969;
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;

%%
tmp5 = textread([filepath,'tmp2.txt']);
tmp5 = tmp5(start_pos:start_pos+N-1,:);
%%
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
tmp2 = tmp2(start_pos:start_pos+N-1,:);
% gpxbr,hfcra,lfcrp
%% 轨距的对比
chaogao = wave_out(:,6);

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

% rou_l = rou_l/(129.01);
% rou_r = rou_r/(129.01);

for i = 3:length(rou_l)
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end

% ******************step2 加速度计的滤波 **********************************
ay = fmctrl_data(:,5);

% ---------------------- 经过滤波器 --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
end
%% 积分
% sita_b = sita_b/3276.8/180*pi;
sita_b = tmp2(:,1);
for i = 3:length(sita_b)
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz - G_par .* sita_b .* tbs.^2 + ht * sita_b_dot2;
camo = -camo;
camo = quzheng(camo);
% camo = floor(camo);%%这里说明了取整带来了问题
% camo = aln(:,1);%%好像是随机取整吧
%%这里说明了camo测的不准，主要是ay有点问题吧
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;

%% 虽然按照他的方法去实现了，但是并没有获得相同的结果，貌似离散了
%% fdatool设计滤波器
    %     load 'filter2.mat'
    %     [h,f] = freqz(Num,1,10000,4);
    %     lamda = 1./f;
    %     figure;semilogx(lamda,(abs(h)));
    %     amcol_fir = conv(Num,amcol);
    %     NN = (length(Num)-1)/2;
    %     amcol_fir = amcol_fir(NN+1:end-NN);
    %     figure;plot(amcol_fir);hold on;plot(amcol)
    %     ydotL = 0;
    %     yL = 0;
    %     for i = 2:length(amcol_fir)
    %         ydotL = ydotL + amcol_fir(i);
    %         ydotL_save(i) = ydotL;
    %         yL(i,1) = yL(i-1,1) + ydotL + 69;
    %     end
%% 窗函数设计滤波器（参考文档）
%% 最简单的滤波器处理（基于窗函数的滤波器）波长大致在20m以下的，没啥区别
%% 对二阶差分进行滤波，如何保证最后两次求和之后，能得到不发散的结果呢？
%% 现在的结果虽然不发散，但是低频的量还是没办法去除掉，所以需要多个窗函数进行及联，生成一个更好的长波滤波器
%% 轨检中的做法
%     M = 30;     %% 窗函数的半窗长，
%     Num = 100;%%队列的长度
%     data = zeros(1,Num);
%     in = 50;
%     in_1 = 49;
%     in_M_1_p = 19;
%     in_M_f = 80;
%     yp = 0;
%     yL = zeros(2,1);
%     dot = 0;sum_x = 0;
%     for i = 1:length(amcol)
%         data(in_M_f) = amcol(i);%%输入数据
%     %     input_project =  data(in) - data(in_1)+ ( -  data(in_M_f) + data(in_M_1_p) )/ (2*M+1);
%     %     yout = yp + input_project;
%         
%         %% 换一种表示的方法
%         sum_x = sum_x + data(in_M_f) - data(in_M_1_p);
%         yout = data(in) - sum_x/(2*M+1);
%             %% 更新数组
%         in_1 = mod(in_1,Num)+1;
%         in_M_f = mod(in_M_f,Num)+1;
%         in_M_1_p = mod(in_M_1_p,Num)+1;
%         in = mod(in,Num)+1;
%         yp = yout;
%         %% 积分
%         yout_save(i,1) = yout;
%         dot = dot + yout;
%         yLdot(i) = dot;
%         yL(i+1) = yL(i) + dot;
%     end

%% 换一种方式
    %         M = 30;     %% 窗函数的半窗长，
    %         Num = 100;%%队列的长度
    %         data = zeros(1,Num);
    %         in = 50;
    %         in_1 = 49;
    %         in_M_1_p = 19;
    %         in_M_f = 80;
    %         yL = zeros(2,1);
    %         dot = 0;sum_x = 0;
    %         center_pos = 50;
    %         for i = 1:length(amcol)
    %             data(find_index(center_pos + M )) = amcol(i);%%输入数据
    %             %% 换一种表示的方法
    %             sum_x = sum_x + data(find_index(center_pos + M )) - data(find_index(center_pos - M - 1));
    %             yout = data(center_pos) - sum_x/(2*M+1);
    %             %% 积分
    %             yout_save(i,1) = yout;
    %             dot = dot + yout;
    %             yLdot(i) = dot;
    %             yL(i+1) = yL(i) + dot;
    %             %%
    %             center_pos = center_pos + 1;
    %             center_pos = find_index(center_pos);
    %         end
%% 实现目前的滤波器
%% 虽然滤波器的实现是没什么问题，但是时延需要调整
    % FSCAL = 0.1
    % -1.036, 0.036, -0.25, 0.25, FSCAL * 2.0
    % 参数设定
    FSCAL = 0.1;
    longWaveFilterBufferSize = 2048;
    mainTriWinCoef = -1.036;
    auxTriWinCoef = 0.036;
    mainRecWinCoef = -0.25;
    auxRecWinCoef = 0.25;
    scaleCoef =  FSCAL * 2.0;
    %% 三角窗和矩形窗并联的方法值得讨论
    %% 参考目前轨检中的做法，参考文档「高速线路轨道长波不平顺检技术研究」
    %% 要具有封装性的代码
    %% 但是这里的处理可能不是很准确或者合理
    %% 可能这里的窗长并不是很合理的
    mainTriM = 160;
    auxTriM = 40;
    mainRecM = 160;
    auxRecM = 280;
    mainTriFront = 0; auxTriFront = 0; mainTriBehind = 0; auxTriBehind = 0;%%三角窗
    mainRecSum = 0; auxRecSum = 0; mainTriSum = 0; auxTriSum = 0;
    MianREC = 0; AuxREC = 0; MainTRI = 0; AuxTRI = 0;
    Num = 2048;     %%循环队列存储数据
    data = zeros(Num,1);
    center_pos = 1000;
    dot = 0; yL = zeros(2,1);
    for i = 1:length(amcol)
        %% 输入数据
        data_input_index = center_pos + auxRecM;
        data(find_index(data_input_index)) = amcol(i);
        %% 矩形窗
        mainRecSum = mainRecSum + data(find_index(center_pos + mainRecM )) - data(find_index(center_pos - mainRecM - 1 ));
        auxRecSum = auxRecSum + data(find_index(center_pos + auxRecM)) -  data(find_index(center_pos - auxRecM - 1 ));
        %% 三角窗
        mainTriFront = mainTriFront + (data(find_index(center_pos + mainTriM-1)) - data(find_index(center_pos-1))) / mainTriM;
        mainTriBehind = mainTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - mainTriM - 1 )) ) / mainTriM;
        
        auxTriFront = auxTriFront + (data(find_index(center_pos + auxTriM-1)) - data(find_index(center_pos-1))) / auxTriM;
        auxTriBehind = auxTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - auxTriM - 1 )) ) / auxTriM;
        
        mainTriSum = mainTriSum + mainTriFront - mainTriBehind;
        auxTriSum = auxTriSum + auxTriFront - auxTriBehind;
        %% 计算更新
        yout = data(find_index(center_pos)) + mainTriWinCoef*mainTriSum/mainTriM + auxTriWinCoef*auxTriSum/auxTriM + ...
            mainRecWinCoef*mainRecSum/(2*mainRecM+1) + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
%         yout = data(center_pos) - 1.036*mainTriSum/mainTriM + 0.036*auxTriSum/auxTriM;
    %     + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
        yout_save(i) = yout;
        dot = dot - yout * 2*auxTriM/(auxTriM*2+1);
        yL(i+1) = yL(i) + dot;
        %% 更新
        center_pos = center_pos + 1;
        center_pos = find_index(center_pos);
    end
yL = yL * scaleCoef;%%除以固定的系数，然后得到轨向的值
%%
% aln(:,4);
figure;plot([yL(126:end)],'LineWidth',1);hold on;plot(aln(:,4));legend matlab70m gj25m;title('matlab计算的70m长波与25m左轨向对比')
% figure;plot((yL - aln(:,4)));%%基本完全一致
%% 这里有个问题，这个延时是怎么确定的？为什么126？是相对于25m波长的126

%% 25m,70m
tmp7 = textread([filepath,'LongWaveResultForAln_L.txt']);
tmp7 = tmp7(start_pos:start_pos+N-1,:);
longwave25m = tmp7(:,1);
longwave70m = tmp7(:,2);
% longwave120m = tmp7(:,3);
for i = 1:length(yL)
    yL(i) = point3filter(yL(i));
end
figure;plot(longwave70m);hold on;plot([zeros(143,1);yL/1.324]);legend gj matlab;title('70m长波对比');set(gca,'Fontname','Times New Roman','fontsize',16);
%%也就是说经过我写的滤波器之后有一个幅值与其不一致


%% 这里长波的问题就是值太大，而且也不是很对齐啊，难道滤波器的实现还是有问题？
yL_2 = longwave_filter(amcol,281,71,281,491);
for i = 1:length(yL_2)
    yL_2(i) = point3filter(yL_2(i));
end
figure;plot(longwave70m);hold on;plot([zeros(143,1);yL/1.324]);legend gj matlab;title('70m长波对比');set(gca,'Fontname','Times New Roman','fontsize',16);


%% 观察频谱(这种方法有点问题)
% plot_mag(longwave25m/103,'gj中给出的轨向 25m');
% plot_mag(longwave70m/103,'gj中给出的轨向 70m','hold');
% plot_mag(yL/103,'matlab中给出的轨向 70m','hold');legend gj25m gj70m 70m/matlab;



%% 函数
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
%% 参考别的方法
function plot_mag1(signal_data , tit , varargin)
%% 波长实际上就是周期
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
%% 在长波滤波之后进行积分？
function out = integrational(auat)

persistent  euas euasp  auatp eusa;
if isempty(euasp)
    euasp =  0;
    auatp = 0;
    eusa = 0;
    euas = 0;
end
gyroabase = 57;

%% 计算
euas = euas + (auat + auatp) * 0.5;
eusa = (2 * gyroabase / (2 * gyroabase + 1) * eusa) + (euas + euasp) * 0.5;

%% 更新
out = eusa;
auatp = auat;
euasp = euas;

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

