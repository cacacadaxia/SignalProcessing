

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
%  功能： 1.
%        2.
%       3. 
%--------------------------------------------------------------------------

%% 1 对左高低数据进行处理
%% 读取数据
%% 左高低数据
% clear all;
d = textread('gplpe.txt');
d = d/129.01;

%% 对采集的信号进行处理，然后再滤波
%% 现在的fs始终都是4Hz，需要注意长波与时间的对应关系
lamda = 100;        %%d
omega = 2*pi/lamda;
l = length(d);
data.time = 0:0.25:((l-1)*0.25);

for i=1:length(d)
    dtp(i) = d(i) + i/1000;
end

    %% 滤波
    % b = Num;
    % sig_out = conv(data.signals.values,b);
    % figure;
    % plot( d,'LineWidth',1 );
    % hold on;
    % plot( dtp,'LineWidth',0.5 );
    % plot(  sig_out(   ((length(b)-1)/2+1):end   ) )
    % 
    % %% 长波滤波器之后
    % figure;
    % fs = 4;     %% 0.25m为一个采样间隔
    % N = length(sig_out);
    % x = (1:N/2+1)/N*fs;
    % x = 1./x;
    % tp = abs(fftshift(fft(sig_out)));
    % tp = 20*log10(tp);
    % semilogx( x , tp((length(tp)/2):end)   );
    % xlabel('\lambda m')
    % ylabel('Mag dB')
    % set(gca,'Fontname','Times New Roman','fontsize',16);
    % %% 那么fs就是4Hz，且 lamda = 1/fs;
    % % d = d-20*[sin(omega*data.time)]';
    % % figure;
    % hold on
    % fs = 4;     %% 0.25m为一个采样间隔
    % N = length(dtp);
    % x = (1:N/2+1)/N*fs;
    % x = 1./x;
    % tp = abs(fftshift(fft(dtp)));
    % tp = 20*log10(tp);
    % semilogx( x , tp((length(tp)/2):end)   );
    % xlabel('\lambda m')
    % ylabel('Mag dB')
    % set(gca,'Fontname','Times New Roman','fontsize',16);


%% 2 激光摄像组件数据
%% 读取数据
%% 激光摄像组件

    % d = textread('jiguang.txt');
    % d = d/129.01;
    % l = length(d);
    % data.time = 0:0.25:((l-1)*0.25);
    % data.time = 1:l;
    % data.signals.values = d;
    % data.signals.dimensions = 2;
    % 
    % sita_bt = (d(:,2)-d(:,1))/1435;
    % figure;plot(sita_bt/pi*180);


%% 3.读取所有的数据
%% 对于所有的数据进行分析，最好是可以在matlab中搭建出整个模型
load_txt
sensor_data = fmctrl_data;
figure;plot(sensor_data);


%% (1)观察所有的频谱
% for i = 1:size(sensor_data,2)
%     
%     signal_data = sensor_data(:,i);
%     figure;
%     fs = 4;     %% 0.25m为一个采样间隔
%     N = length(signal_data);
%     x = (1:N/2+1)/N*fs;
%     x = 1./x;
%     tp = abs(fftshift(fft(signal_data)));
%     tp = 20*log10(tp);
%     semilogx( x , tp((length(tp)/2):end)   );
%     xlabel('\lambda m')
%     ylabel('Mag dB')
%     set(gca,'Fontname','Times New Roman','fontsize',16);
% end

%% 分析信号

% 1。速度
tbs = sensor_data(:,end);
tbs_s = tbs/1e5;
v = 0.9./tbs_s;
distance = 0:0.25:0.25*(length(tbs)-1);
figure;plot(distance , v);
title('速度 km/h')

% 2。加速度
acc_z = sensor_data(:,7);

t = 0;
for i=2:length(tbs_s)
    t(i) = t(i-1) + tbs_s(i);
end
figure;plot(t./3600);
% 大概0.6小时


