% =========================================================================
%
%                       超高程序的理解
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月10日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.整合超高的低频和高频的计算，并且根据自己的理解修改
%        2.
%        3. 
% 注意：
%       1.1.31与1.31两个值一致，这一点需要注意，一个是陀螺仪的单位，另一个是omega1的倒数
%       2.int(gj programme) = 量纲*float(true)
%       3.
%--------------------------------------------------------------------------
%%
clear all;
close all;
%%
N = 10000;
tmp2 = textread('Ps3_filter_wx.txt');
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
gpro = tmp2(:,1);
%% 读取yaw rate
tmp3 = textread('fmctrl_data_1337.txt');
if length(tmp3)>N
    tmp3 = tmp3(1:N,:);
end
yaw = tmp3(:,6);
gpin = tmp3(:,4);
%%
tmp4 = textread('Hz_filter_inc.txt');
if length(tmp4)>N
    tmp4 = tmp4(1:N,:);
end
inc_comp = tmp4(:,1);
lfcrp_comp = tmp4(:,3);
%%
tmp5 = textread('Bz_filter.txt');
if length(tmp5)>N
    tmp5 = tmp5(1:N,:);
end
gpin_2 = tmp5(:,1);
infp_2 = tmp5(:,2);
%% 读取数据
tmp = textread('fz_filter_gaodi.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end
gpvlo = tmp(:,1);%%32768/5 inch
gpvro = tmp(:,2);
dt74 = gpvro/2-gpvlo/2;%%为什么是除以2 %%%%32768/10 inch
fd74 = -dt74;       %%不清楚
result = tmp(:,4);
%%

%% 计算
%% 这些参量分别都是什么意思？
% hight = 3.820;
% fim = 0.820210;     %%0.25/0.3048 1feet等于0.3048m
fim = 0.25;%%m
dgyro = 6.56 * 0.3048;%%m       %%陀螺的位置
compf = -1;
scali = 1230.2646;
mytbp = [];

DDISP = 59.5 * 0.3048;           %%DDISP*0.0254 = 1.5113 这个吧？
% sd74 = 126376.6875;     %%这才是正解，这是啥意思？
% sd74 = 100000.0 * ( 75.194127 / DDISP );
% % SD74 IS 10**5 * (5*0.4*3.2768/L*SIN 5)
sd74 =  ( 75.194127 / DDISP );%%sd74 1e5或者1都不影响，对最终结果影响很小
%%sd74 = 246083.654/3276.8，这一点需要注意，就是为了让位移计的测量出来的角速度单位换算到与陀螺仪的角速度单位一致用的

%% 经过打印，确定过后的参数
degti = 1.0375;
%% 都转换成为m
hight = 1.8 * 0.3048;         %%feet   
% scali = 1/246083.654/0.54426*1638.4;%%这个值改变还是有影响的
dincl = 1.32 * 0.3048;       %%倾角计到转向架的距离
mfd12 = -(dgyro*2)*(dgyro*2)/12;%%？  --> feet
yaw_dotp = 0;

for i = 2:length(result)
    % ---------------------------------
    %%其中dot代表差分，diff代表微分
    tbs = tmp(i,3);
    fd74_dot = (fd74(i) - fd74(i-1))*sd74;  %%--> 量纲 246083.654
    fd74_diff = fd74_dot / ( tbs / 1e5 ); %%tbs/1e5代表T
    w_bt(i,1) = B1(fd74_diff,tbs);
    w_bt_dot = w_bt(i,1) - w_bt(i-1,1);
    wt(i,1) = gpro(i) - w_bt(i);        %%他们两个单位统一吗？应该是
    wt_dot = wt(i,1) - wt(i-1,1);
    
    yaw_dot = yaw(i) - yaw(i-1);
    yaw_dot2 = yaw_dot - yaw_dotp;
    % ----------------------------------
    % 这是一个很复杂的量，逻辑就是将所有的量都转换成为feet就好了
    dtemp = fim  * yaw(i) * compf / (tbs/1e5);%% wz -->(fim -- 0.25变成 feet ，fim/tbs代表速度)
    dtemp = dtemp - dgyro * yaw_dot / (tbs/1e5);%% L*s*F(s)-->(yaw_dot / tbs代表微分diff)
    dtemp = dtemp +  mfd12 * yaw_dot2 / fim / (tbs/1e5);%% 不要也罢，二阶的量，对结果影响很小
    dtemp = dtemp + hight  * wt_dot / (tbs/1e5);%% ht*s*F(s) %% wt_dot / tbs-->(wt_dot / tbs 代表微分diff)
%     dtemp = scali * ( dtemp + dincl * w_bt_dot /  (tbs/1e5) );%% -->( w_bt_dot / tbs 代表微分diff)
%     dtemp = (dtemp + dincl * w_bt_dot /  (tbs/1e5)) / 246083.654 /( 0.016932 * 9.8 ) * 3276.8 * 1.0057; %%1.0057是前后两个scali的参数比值
    dtemp = ( dtemp + dincl * w_bt_dot /  (tbs/1e5) ) / 246083.654 / 9.8 * 3276.8; %%1.0057是前后两个scali的参数比值
 
    %% sacli-->/246083是将陀螺仪转换成为deg/s的单位？
    %% 0.54426-->0.016932g/inch
    %%
    dtmp_Fz = F(dtemp , tbs);    %%这个好像有问题？
    out(i,1) = dtmp_Fz;
    %% next step
    %     1/(0.016932*9.8) * 3276.8 = 19748
    infp = B2( gpin(i) / ( 3276.8/(0.016932*9.8) ) * 3276.8, tbs);   %%这个B多次使用所以产生影响  %%3276.8*rad作为量纲
    infp = infp / 9.8;      %%除以g常量
    infp_save(i,1) = infp;
    inc = infp + dtmp_Fz;     
    lfcrp(i,1) = H3z( inc , tbs );
    %%转换量纲
    lfcrp(i,1) = lfcrp(i,1) / 3276.8 * 1638.4 * 180 / pi; %%量纲都换成1638.4*deg 就行,与轨检程序中的结果相符合
    %% update, last step
    yaw_dotp = yaw_dot;
    %%
    %%在dtemp多处理了0.016932，在gpin中少处理了0.016932。同时如果转换成为rad的话，那么又少了180/pi转换成为deg
    %%两个值相抵消，那么就多出了(0.016932*180/pi)==0.97这样的一个倍数
    
end
% lfcrp = lfcrp/1638.4/180*pi;
% lfcrp;%%度数
figure;plot(lfcrp /(0.016932*180/pi) - lfcrp_comp);     %%这一点需要注意一下，与原本程序的区别
figure; plot(lfcrp / 1638.4);
%%直接变成的


%% % wx,tbs,hfcra,t4,...
start_pos = 1;
tmp = textread('Ps3_filter_wx.txt');
if length(tmp)>N
    tmp = tmp(start_pos : start_pos+N-1,:);
end
hfcra_ref = tmp(:,3);
TBS = tmp(:,2);
GYRO = tmp(:,1);

%%
tmp4 = textread('Hz_filter_inc.txt');
if length(tmp4)>N
    tmp4 = tmp4(start_pos : start_pos+N-1,:);
end
lfcrp_ref = tmp4(:,3);


%% 应该除以一个系数，这里为什么没有？
s1 = 246083.654;
GYRO = GYRO  / 3276.8 / 1.31 * 1638.4; %%转换成为 3276.8*deg/s --> 1638.4*deg/s，积分之后正好结果是deg，结果就与最终对的上了
hfcra_func = Pz3_gj( GYRO , TBS );
figure;plot( hfcra_ref -  hfcra_func );figure;plot( hfcra_ref );hold on;plot( hfcra_func );%%不准

%%
reslut = hfcra_func + lfcrp / (0.016932*180/pi);%%多了一个这个
figure;plot(reslut);
reslut_ref = lfcrp_comp + hfcra_ref;
hold on;plot(reslut_ref)
legend matlab gj

%%

function out = B1(x_k,tbs)
persistent y;
if isempty(y)
    y = zeros(2,1);
end
y(2) = ( y(1)*2^17 + tbs*x_k )/(2^17 + tbs );
% tbs*x_k
%%
y(1) = y(2);
out = y(2);
end

function out = B2(x_k,tbs)
persistent y;
if isempty(y)
    y = zeros(2,1);
end
y(2) = ( y(1)*2^17 + tbs*x_k )/(2^17 + tbs );
% tbs*x_k
%%
y(1) = y(2);
out = y(2);
end

function out = F(x,tbs)
%% 滤波器设定
%% 对于x[3]，其与同理
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2

%% 这个滤波器存在延时，几十不等？
%% 从哪里获得的？
persistent y;
if isempty(y)
    y = zeros(3,1);
end
y(3) = ( y(2)*(2*2^28 + 2^14*tbs) - y(1)*2^28 + tbs^2 * x  )/(2^28 + 2^14*tbs + tbs^2);
%% 更新temp量
y(1) = y(2);
y(2) = y(3);
out = y(3);

end

function out = H3z(x_k,tbs)

persistent x y;

if isempty(x)
    x = zeros(2,1);
    y = zeros(3,1);
end
x(2) = x_k;
x_dot = x(2) - x(1);
up = tbs*( 2^18 * (1+1/2+1/2^4) *x_dot + 2^18*y(2) + tbs*x(2)) + 2^36*(2*y(2) - y(1));
down = 2^36 + tbs*(2^18 + tbs);
y(3) = up/down;

%% 更新
y(1) = y(2);
y(2) = y(3);
x(1) = x(2);

out = y(3);
end

function out = P3z(x_k,tbs_k)


persistent x y tbs;

Omega1 = 10^5/2^17;
Omega2 = 10^5/2^14;
Omega3 = 10^5/2^18;

zeta = Omega1+Omega2+Omega3;
nang = Omega2^2 + Omega2*Omega3 + Omega2*Omega1 + Omega1*Omega3 +Omega3^2;
kai = (Omega1 + Omega3)*Omega2^2 + Omega2*Omega1*Omega3 + (Omega1 + Omega2)*Omega3^2;

if isempty(x)
    x = zeros(5,1);
    y = zeros(5,1);
    tbs = zeros(2,1);
end
x(5) = x_k;
tbs(2) = tbs_k;


T = tbs(2)/1e5;         %%这里注意
x_dot = x(5) - x(4);
x_dot2 = x(5) - 2*x(4) + x(3);
x_dot3 = x(5) - 3*x(4) + 3*x(3) - x(4);
x_dot4 = x(5) - 4*x(4) + 6*x(3) - 4*x(2) + x(1);
N1 = kai*T^3*x_dot + nang*T^2*x_dot2 + zeta*T*x_dot3 + x_dot4;

y_k_1_dot = y(4) - y(3);
y_k_2_dot = y(3) - y(2);
y_k_3_dot = y(2) - y(1);
y_k_1 = y(4);       %%n-1时刻的y
N2_1 = 2^19*(1/2+1/2^5)*tbs(2)^3*y_k_1;
N2_2 = 2^37*(1/2+1/2^5+1/2^9) *tbs(2)*( tbs(1) * y_k_1_dot + tbs(2)*y_k_1 );
N2_3 = 2^51*(1/2+1/2^5) * ( tbs(1)* ( y_k_1_dot + y_k_1_dot - y_k_2_dot ) + tbs(2)*y_k_1 );
N2_4 = 2^64*( 3*( y_k_1_dot - y_k_2_dot ) + y_k_3_dot + y_k_1 );
N2 = 2^(-64) * ( N2_1 + N2_2 + N2_3 + N2_4 );

up = 2^64*( N1/Omega1 + N2 );
down = 2^64 + tbs(2) *2^51 *(1/2+1/2^5) + tbs(2)^2*2^51*(1/2+1/2^5+1/2^9) + tbs(2)^3*2^19*(1/2+1/2^5) + tbs(2)^4;

y(5) = up/down;
%% 更新
x(1) = x(2);
x(2) = x(3);
x(3) = x(4);
x(4) = x(5);

y(1) = y(2);
y(2) = y(3);
y(3) = y(4);
y(4) = y(5);

tbs(1) = tbs(2);
out = y(5);
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
x_dot = x(3)-x(2);
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
%%
function out = B(x_k,tbs)
persistent  y;
if isempty(y)
    y = zeros(2,1);
end
y(2) = ( y(1)*2^17 + tbs*x_k )/(2^17 + tbs );

%%
y(1) = y(2);
out = y(2);

end
%%
function out = Pz3_gj(x1 , x2_tbs)

%%
% x1 = tmp(:,1);
% x2_tbs = tmp(:,2);
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
x_k_1 = 0;
x_k_1_dot = 0;
x1p2p = 0;
x2p3p = 0;
y_k_2 = 0;
fr3pp = 0;
y_k_3_dot = 0;
y_k_1 = 0;
ftrp2 = 0; 
p1flg = 0;
y_k_2_dot = 0;
%% 初始化
rldb_1 = x1(1,1)*degti;
rldb_2 = x1(2,1)*degti;
rldb_3 = x1(3,1)*degti;
x_k_1 = rldb_3;
x_k_1_dot = rldb_3 - rldb_2;
x1p2p = rldb_3 - 2*rldb_2 + rldb_1;
INILIZED = 0;
%%
x_k_1 = 0;
x_k_1_dot = 0;
x1p2p = 0;
x_k_2_dot2 = 0;
tbspp = 0;
y_k_1_dot = 0;
for i = 4:length(x1)
    %%
    rldd = x1(i);
    tbs(2) = x2_tbs(i);
    dnom4 = ( tbs(2)  + ct1 ) * tbs(2);
    dnom4 = ( dnom4 + ct2 ) * tbs(2);
    dnom4 = ( dnom4 + ct3 ) * tbs(2);
    dnom4 = dnom4 + ct4;
    x_k = rldd * degti;%%why？
    %%
    x_k_dot = x_k-x_k_1;
    %%/* SECOND */
    x_k_1_dot2= x_k_dot-x_k_1_dot;
    %%/* THIRD */
    N1 = x_k_1_dot2 - x_k_2_dot2;
    %/* Multiply by 2 ** 62 */
    N1 = ct4x*N1;
    x2p3 = x_k_dot * tbs(2);
    %%/* Do another, Multiply by 2 ** 49 * (1/2 + 1/2 ** 3) */
    N1 = N1 + cx7*(x_k_dot * tbs(2)-x_k_1_dot * tbs(1));
    N1 = N1 + cx8*(tbs(2)*x_k_dot * tbs(2));
    N1 = N1 + cx9*tbs(2)*tbs(2)*tbs(2)*x_k;
    if (y_k_2 == 0)&&(INILIZED==0)
        y_k_2 = y_k_1;
        INILIZED = 1;
    end
    %%
    y_k_1_dot = y_k_1 - y_k_2;
    N2 = 0 + ct4x*(3*(y_k_1_dot - y_k_2_dot) + y_k_3_dot + y_k_1);
    N2 = N2 + cx4*(tbs(1) * (y_k_1_dot + y_k_1_dot) - tbspp*y_k_2_dot + y_k_1 * tbs(2));    %%注意这个量ftrp2
    N2 = N2 + cx2*tbs(2)*( tbs(1) * y_k_1_dot + y_k_1 * tbs(2) );
    N2 = N2 + cx5*tbs(2)*tbs(2)*y_k_1*tbs(2);
    N1 = N1 / 0.76;%%注意这一点
    up = N1 + N2;
    t4 = up / dnom4;
    x_k_1 = x_k;
    x_k_1_dot = x_k_dot;
    x_k_2_dot2 = x_k_1_dot2;
%%
    y_k_3_dot = y_k_2_dot;%%y dot i-3
    y_k_2_dot = y_k_1_dot;%%y dot i-2
    
    %% 差分
    hfcra(i,1) = (t4 - y_k_1);
    y_k_2 = y_k_1;%%y i-2
    y_k_1 = t4;%%y i-1
    %% 更新
    tbspp = tbs(1);
    tbs(1) = tbs(2);
end

%%
out = hfcra;
end
