% =========================================================================
%
%                  验证滤波器的准确性
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月18日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.Ps这里有问题，数值太小
%        2.不准的问题基本解决，但是在初始化的过程中出现的问题来自于：
%               初始化做的不好
%        3.滤波器的输出的初始化怎么办？
%        4.分析角度大小：发现单位始终是一个问题，不能准确对的上
%        5.
%
%--------------------------------------------------------------------------
clear all;
close all;
N = 10000;
tmp = textread('Ps3_filter_wx.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end
% wx,tbs,hfcra,t4,...

%%
tmp2 = textread('tmp_zhongjian_1337.txt');
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
% gpxbr,hfcra,lfcrp

%% 读取yaw rate
tmp3 = textread('fmctrl_data_1337.txt');
if length(tmp3)>N
    tmp3 = tmp3(1:N,:);
end
gpvlo = tmp3(:,8);
gpvro = tmp3(:,10);
% 见excel表格

%%
ct1 = 278528.0;	 %%/* ct1=2**19*(1/2+1/2**5) */
ct2 = 7.3282879e10;  %%/* ct2=2**37*(1/2+1/2**5+1/2**9) */
ct3 = 1.1962686e15;  %%/* ct3=2**51*(1/2+1/2**5) */
ct4 = 1.8446744e19;  %%/* ct4=2**65*(1/2) */
ct4x = 1.8446744e19; %%/* ct4x=ct4 */
cx0 = 262144.0;	   
cx1 = 409600.0;	   
cx2 = 7.3282879e10;  
cx3 = 3.0;		   
cx4 = 1.1962686e15;  
cx5 = 278528.0;	   
cx6 = 0.5;		  
cx7 = 1.3370061e15; 
cx8 = 8.2409684e10;  
cx9 = 837632.0;	   
cx10 = 6.8719476e10; 
%%  其余参数
degti = 1.037249;


%% 初始化的问题需要避免
rldbp = 0;
x1p3p = 0;
x1p2p = 0;
x2p3p = 0;
frlpp = 0;
fr3pp = 0;
frdp3 = 0;
frlp = 0;
ftrp2 = 0; 
p1flg = 0;
frdp2 = 0;
%% 初始化

rldb_1 = tmp(1,1)*degti;
rldb_2 = tmp(2,1)*degti;
rldb_3 = tmp(3,1)*degti;
rldbp = rldb_3;
x1p3p = rldb_3 - rldb_2;
x1p2p = rldb_3 - 2*rldb_2 + rldb_1;
INILIZED = 0;
for i = 4:length(tmp)
    
    rldd = tmp(i,1);
    ttbs = tmp(i,2);
    xtbs = ttbs;
    xtbsp = tmp(i-1,2);%%这一步就很重要
    
    dnom4 = ( ttbs  + ct1 ) * ttbs;
    dnom4 = ( dnom4 + ct2 ) * ttbs;
    dnom4 = ( dnom4 + ct3 ) * ttbs;
    dnom4 = dnom4 + ct4;
    %%//  dnom4=ct4+ttbs*(ct3+ttbs*(ct2+ttbs*(ct1+ttbs)));
    
    %%/* numerator part X1 */
    rldb = rldd*degti;
    %%
    x1p3=rldb-rldbp;
    %%/* SECOND */
    x1p2=x1p3-x1p3p;
    %%/* THIRD */
    t3 = x1p2-x1p2p;
    
    
    %% /* Multiply by 2 ** 62 */
    t3 = ct4x*t3;
    
    
    x2p3 = x1p3 * xtbs;
    
    %%/* Do another, Multiply by 2 ** 49 * (1/2 + 1/2 ** 3) */
    t3 = t3 + cx7*(x2p3-x2p3p);
    
    
    t3 = t3 + cx8*(xtbs*x2p3);
    
    
    t3 = t3 + cx9*xtbs*xtbs*xtbs*rldb;
    if (frlpp == 0)&&(INILIZED==0)
        frlpp = frlp;
        
        INILIZED = 1;
    end
    frdp1=frlp-frlpp;
    fr3p=cx3*frdp1;
    
    t3 = t3+ct4x*(fr3p-fr3pp+frdp3+frlp);
    ftrp = frlp*xtbs;
    ftrp1 = xtbsp*frdp1;
    t3=t3+cx4*(ftrp1+ftrp1-ftrp2+ftrp);
    
    
    t3=t3+cx2*xtbs*(ftrp1+ftrp);
    t3=t3+cx5*xtbs*xtbs*ftrp;
    %%/* Result */
    t4 = t3/dnom4;
    %%/* ripple recursive computations, filtered roll rate evaluation */
    rldbp=rldb;
    x1p3p=x1p3;
    x1p2p=x1p2;
    x2p3p=x2p3;
%     x2p2p=x2p2;
    frdp3=frdp2;
    frdp2=frdp1;
    fr3pp=fr3p;
    ftrp2=ftrp1;
    
    %%/*Previous result */
    frlpp = frlp;
    %%/*Result */
    frlp = t4;
    frlp_save(i) = t4;
    %%  /* check for roll filter first pass case */
    if(p1flg==0)
        p1flg=-1;
        rollp=frlp;
    end
    
    t3 = frlp - rollp;
    rollp = frlp;
    t3 = t3*cx6;
    hfcra(i,1) = t3;
end

%% 角度分析（from gj）

gpxbr = tmp2(:,1);      %%sita_b
hfcra = tmp2(:,2);      %%高频
lfcrp = tmp2(:,3);      %%低频

dt74 = (gpvro - gpvlo)/2;

scalr = 0.5;
gpxl = gpxbr + dt74*scalr;      %%这就是最终的角度？那么角度是多少？

gpxl = gpxl/32.2520/100;        %%单位与文档中一致
gpxbr = gpxbr/32.2520/100;      %%单位为度
figure;plot(gpxl);

%% 陀螺仪数据分析 wx
wx = tmp(:,1)/3276.8/180*pi;

% 积分
sita = 0;

for i = 1:length(wx)
    
    tbs = tmp(i,2)/1e5;
    
    sita = sita + wx(i)*tbs;
    
    sita_save(i,1) = sita;
end

%%角度对比
figure;plot(gpxbr)
hold on;plot(sita_save/pi*180);
legend gj 直接积分

% figure;plot(gyro/pi*180);legend 1 2 3%%感觉这里的单位是错误的，尤其是角速度的单位

%%

sita_bt = dt74*scalr;
figure;plot(sita_bt/32.2520/100);%%sita_bt只有0.1度？


%% 惯组
gyroll = tmp3(:,1)/3276.8/180*pi;
gypitch = tmp3(:,2)/3276.8/180*pi;
gyyaw = tmp3(:,6)/3276.8/180*pi;
accx = tmp3(:,3)/129.01;%%这单位到底是什么？
accy = tmp3(:,4)/129.01;
accz = tmp3(:,7)/129.01;

yaw = 0;pitch = 0; roll = 0;
for i = 1:length(gyyaw)
   gyyaw_ = gyyaw(i);
   gyroll_ = gyroll(i);
   gypitch_ = gypitch(i);
   tbs = tmp(i,2)/1e5;
   yaw = yaw + tbs*gyyaw_;
   pitch = pitch + tbs*gypitch_;
   roll = roll + tbs*gyroll_;
   yaw_save(i) = yaw;
   pitch_save(i) = pitch;
   roll_save(i) = roll;
end
figure;plot(yaw_save/pi*180);
figure;plot(pitch_save/pi*180);
figure;plot(roll_save/pi*180);
%%这里的问题是直接对yaw的积分是十分不准确的，相当于列车在1km之内转弯为70度，那是非
%%常离谱的，待后续修改

%% 角速度积分
%% 角度积分依然存在问题
gyro = [gyroll,gypitch,gyyaw];
angle = zeros(1,3);
for t = 2:length(gyro)
    tbs = tmp(i,2)/1e5;
    ang_k = angle(t-1,:);
    del = cau_w(ang_k(1),ang_k(2),ang_k(3))*gyro(t,:)';
    angle(t,:) = (ang_k + del.'*tbs);
end
figure;plot(angle(:,1)/pi*180);%%deg，这就是积分带来的问题？
hold on;
plot(gpxbr)
legend 1 2 3


%% 加速度积分
TBS = tmp(:,2);


v = 0.25./(TBS./1e5);
for i = 21:length(v)
    tbs = TBS(i)/1e5;
    a(i) = (v(i) - v(i-20)) / tbs / 20;
end


v_ac = 0;
for i = 2:length(accx)
    accx_ = accx(i);
    tbs = TBS(i)/1e5;
    v_ac(i) = v_ac(i-1)+accx_*tbs;
end





%% 函数
 % figure;plot(gyro);legend 1 2 3
function W = cau_w(fai,sita,psi)
W = [1,sin(fai)*tan(sita) , cos(fai)*tan(sita);
    0, cos(fai) ,       -sin(fai);
    0, sin(fai)/cos(sita) , cos(fai)/cos(sita)];
end





