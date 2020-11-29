% =========================================================================
%
%                  对比两个车的数据（京津高铁）
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月31日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.这里用的是京津的数据，这里暂时没有波形文件，暂时完成数据的对齐工作，找到
%               对应的数据，以供后面的长波对比分析
%        2.小车数据多，轨检数据少。问题是这个小车的采样间隔是0.25m吗？
% 11.24
%       3. 
%--------------------------------------------------------------------------

%% step1: 读取小车数据
clear;
close all;
data = load('gaodi_xiaoche_jj.txt');
gpproLTrolleyOri = data(:,2);
xSmallCar = data(:,1);
ll = xSmallCar(2:end)-xSmallCar(1:end-1);
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolleyOri);
gpproLTrolley = out(len+1:end-len);

%% step2: 截断数据
N = 5000;
gpproLTrolley = gpproLTrolley(81572 - N : 81572 + 1538 + N );
xSmallCar = xSmallCar(81572 - N : 81572 + 1538 + N);
gpproLTrolleyOriCut = gpproLTrolleyOri(81572 - N : 81572 + 1538 + N);

%% step3: 将手推车的0.65m转换成为0.25m，与轨检同步
gpproLTrolley25m = [];
for i = 1:floor(length(gpproLTrolley)*0.65*4)-4
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolley(p1) + p2/0.65*( gpproLTrolley(p1+1) - gpproLTrolley(p1) );
end
% figure;plot(0:0.25:(length(gpproLTrolley25m)-1)*0.25 , gpproLTrolley25m,'-r');
% hold on;plot(0:0.65:(length(gpproLTrolley)-1)*0.65 , gpproLTrolley,'--k');

%% step4: 读取轨检车数据
 data2 = load('gaodi_guijian_jj.txt');
 gpproGjOri = data2(:,6);
 xGj = data2(:,1);
% figure1 = figure('Color',[1 1 1]);
% plot(xGj,gpproGjOri);
% %%这个数据对上了，按照excel中的里程信息
% hold on;plot(xSmallCar , gpproLTrolley);
% set(gca,'Fontname','Times New Roman','fontsize',14);

%% step5: 利用相关性找到对齐点
out = conv( gpproLTrolley25m , gpproGjOri(end:-1:1));
[~,index] = max(out);
index = 17875;
index = index - length(gpproGjOri);
figure1 = figure('Color',[1 1 1]);
plot( gpproLTrolley25m);hold on;plot([zeros(1,index) , gpproGjOri'])
xlabel('采样点 0.25m');ylabel('高低 mm');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
legend 小车 轨检


%% step6: 分析数据
gpproLTrolleyOriCut = gpproLTrolleyOri;
gpproLTrolley25m = [];
for i = 1:floor(length(gpproLTrolleyOriCut)*0.65*4)-4
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolleyOriCut(p1) + p2/0.65*( gpproLTrolleyOriCut(p1+1) - gpproLTrolleyOriCut(p1) );
end
b = load('filter_250m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolley25m);
gpproLTrolley_1000m = out(len+1:end-len);
% figure1 = figure('Color',[1 1 1]);plot(gpproLTrolleyOriCut+196);hold on ;plot(gpproLTrolley);
% legend 原始数据 滤波后;
% set(gca,'Fontname','Times New Roman','fontsize',14);


%% step7: 分析不同波长的波形对比
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolley25m);
gpproLTrolley_120m = out(len+1:end-len);
figure1 = figure('Color',[1 1 1]);plot(gpproLTrolley_120m);hold on ;plot(gpproLTrolley_1000m);
legend 120m 250m;
set(gca,'Fontname','Times New Roman','fontsize',14);


%% step8：分析频谱
% figure;plot(gpproLTrolleyOri);
% plot_mag(gpproLTrolleyOri,'轨检小车数据');
% plot_mag(gpproLTrolley_1000m,'轨检小车数据');
% plot_mag(gpproLTrolley_120m,'轨检小车数据','hold');
% legend 250m 120m
%% 点头角速度
% wave1 = gpproLTrolleyOri./1e3;
% N = 5*4;%%5m
% pitch_t_1 = atan( (wave1(N+1:end) - wave1(1:end-N))./0.25./N )/pi*180;
% N = 10*4;%%10m
% pitch_t_2 = atan( (wave1(N+1:end) - wave1(1:end-N))./0.25./N )/pi*180;
% figure1 = figure('Color',[1 1 1]);
% plot(pitch_t_1);hold on;plot(pitch_t_2);legend L=5m L=10m
% xlabel('采样点 (0.25m)');ylabel('角度 (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% 
% %%可见轨面的点头角度都挺小的，即便是在波长比较大的情况下

%%
% plot_mag(pitch_t_1,'角度频谱');
% plot_mag2(pitch_t_1,'角度频谱');


%%
function plot_mag(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end

%% 能量归一化
% E = sum(abs(signal_data).^2);
% signal_data = signal_data./sqrt(E);

%%
fs = 4;     %% 0.25m为一个采样间隔
signal_data = signal_data./1e3;
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
title(tit);
grid on;
end

function plot_mag2(signal_data , tit , varargin)
if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure1 = figure('Color',[1 1 1]);
end

%% 能量归一化
% E = sum(abs(signal_data).^2);
% signal_data = signal_data./sqrt(E);

%%
fs = 4;     %% 0.25m为一个采样间隔
signal_data = signal_data./1e3;
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
plot( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
title(tit);
grid on;
end


