
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
%  功能： 1.转换成为标准单位，并进行测试(目前这个文件是有问题的)
%        2.为什么x轴向的加速度就是完全错的？
%        3.
%
%--------------------------------------------------------------------------
close all
clear all

%% 读取数据
N = 5e4;
start_pos = 14;
filepath = 'data/0916_1337_x/';
tmp3 = textread([filepath,'fmctrl_data_1337.txt']);
if length(tmp3)>N
    tmp3 = tmp3(start_pos:start_pos+N-1,:);
end

%%
tmp2 = textread([filepath , 'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(start_pos : start_pos+N-1,:);
end
% gpxbr,hfcra,lfcrp

%% 惯组单位换算
gyroll = tmp3(:,1)/246083.654;  %%rad/s
gypitch = tmp3(:,2)/246083.654; %%rad/s
gyyaw = tmp3(:,6)/246083.654;   %%rad/s
accx = tmp3(:,3)/19748;%%m/s^2
accy = tmp3(:,4)/13852;%%m/s^2
accz = tmp3(:,7)/13852;%%m/s^2
TBS = tmp3(:,16)/1e5;
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
%%
figure;plot(yaw_save/pi*180);title('yaw')
figure;plot(pitch_save/pi*180);title('pitch')       %%这显然看起来不太像
figure;plot(roll_save/pi*180);title('roll');        %%roll存在一定漂移


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
figure;plot(licheng,v*3.6,'r','LineWidth',1);grid on;xlabel('里程 km');ylabel('速度 km/h');
set(gca,'Fontname','Times New Roman','fontsize',16);



%% 计算ax，这里的单位是m/s^2

for i = 2:length(TBS)
    accx_tbs(i) = (v_filter(i)-v_filter(i-1)) / TBS(i);
end
% figure;plot(accx_tbs);


%%

for i = 2:length(TBS)
    gyro_tmp(i) = 0.5*yaw_save(i)*(TBS(i) - TBS(i-1))/( TBS(i) + TBS(i-1) );
end
