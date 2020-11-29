
% =========================================================================
%
%                  惯组数据处理
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月24日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.转换成为标准单位，并进行测试
%        2.为什么x轴向的加速度就是完全错的？这个问题保留
%        3.
%   10.15
%           1. 经过修改之后，应该问题不大
%
%--------------------------------------------------------------------------


%%
close all
clear all
% 读取数据
N = 5e5;
start_pos = 1;
filepath = 'data/0916_1337_x/';
fmctrl_data = textread([filepath,'fmctrl_data_1337.txt']);
if length(fmctrl_data)>N
    fmctrl_data = fmctrl_data(start_pos:start_pos+N-1,:);
end
tmp2 = textread([filepath , 'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(start_pos : start_pos+N-1,:);
end
% gpxbr,hfcra,lfcrp

%% 位移计
gpvlo = fmctrl_data(:, 8);
gpvro = fmctrl_data(:,10);
gpvlo = gpvlo/258.016;
gpvro = gpvro/258.016;%%mm
G = 59.5 * 25.4;
mm2o = 1/G*180/pi;
sita_bt = (gpvro-gpvlo) * mm2o;
figure1 = figure('Color',[1 1 1]);plot(sita_bt);xlabel('里程 /0.25m');ylabel('\theta_{bt} /度');set(gca,'Fontname','Times New Roman','fontsize',16);

%% 惯组单位换算
gyroll = fmctrl_data(:,1)/246083.654;  %%rad/s
gypitch = fmctrl_data(:,2)/246083.654; %%rad/s
gyyaw = fmctrl_data(:,6)/246083.654;   %%rad/s
accx = fmctrl_data(:,3)/19748;%%m/s^2
accy = fmctrl_data(:,4)/13852;%%m/s^2
accz = fmctrl_data(:,7)/13852;%%m/s^2
TBS = fmctrl_data(:,16)/1e5;

figure;plot(accy);hold on;plot(accx*2);legend 1 2;
%%这里有一个问题，为什么这里的ax的数据完全是错的？
%%这一点需要注意。
%%之后需要将这两者的结果进行对比。


%% 加速度计的比较

accx_body = fmctrl_data(:,13)/( 32768/2/9.8 );%%m/s^2
accy_body = fmctrl_data(:,14)/( 32768/2/9.8 );%%m/s^2
accz_body = fmctrl_data(:,15)/( 32768/2/9.8 );%%m/s^2
% figure;plot(accx);hold on; plot(accx_body);legend 检测梁 车体
% figure;plot(accy);hold on; plot(accy_body);legend 检测梁 车体
% figure;plot(accz);hold on; plot(accz_body);legend 检测梁 车体
%%检测梁的数据已经是滤完波的了，所以检测梁的高频分量要少一些

%% 算角度
yaw = 0;pitch = 0; roll = 0;
for i = 1:length(gyyaw)
   gyyaw_ = gyyaw(i);
   gyroll_ = gyroll(i);
   gypitch_ = gypitch(i);
   tbs = TBS(i);
   yaw = yaw + tbs * gyyaw_;
   pitch = pitch + tbs*gypitch_;
   roll = roll + tbs*gyroll_;
   yaw_save(i) = yaw;
   pitch_save(i) = pitch;
   roll_save(i) = roll;
end

figure;plot(yaw_save/pi*180);title('yaw')
figure;plot(pitch_save/pi*180);title('pitch')       %%这显然看起来不太像
figure;plot(roll_save/pi*180);title('roll');        %%roll存在一定漂移
%%yaw角度没啥用？

%% 算速度 
%%暂时不用管，这里明显是有问题的，但是为什么直接对acc进行积分是没有用的呢？
%%因为包含了各种角度吧是，但是也不对，车体就是一直沿着x轴向的方向运行的
for i = 1:length(accx)
    tbs = TBS(i);
    v(i) = 0.25/tbs;
end
vtp = 0;
for i = 2:length(accx)
    accx_ = accx(i);
    tbs = TBS(i);
    vtp(i) = vtp(i-1) + tbs * accx_;
end
% figure;plot(vtp);hold on;plot(v);legend acc true

%% 对于v进行滤波（三点滤波）
b = ones(1,150)/150;
v_filter = conv(b,v);
figure;plot(v);hold on; plot(v_filter);
licheng = 0:0.25:(N-1)*0.25;
licheng = licheng/1e3;
% figure;plot(licheng,v*3.6,'r','LineWidth',1);grid on;xlabel('里程 km');ylabel('速度 km/h');
figure1 = figure('Color',[1 1 1]);plot(v*3.6,'r','LineWidth',1);grid on;xlabel('里程 /0.25m');ylabel('速度 km/h');
set(gca,'Fontname','Times New Roman','fontsize',16);

%% 计算ax，这里的单位是m/s^2

for i = 2:length(TBS)
    accx_tbs(i) = (v_filter(i)-v_filter(i-1)) / TBS(i);
end
figure;plot(accx_tbs);




