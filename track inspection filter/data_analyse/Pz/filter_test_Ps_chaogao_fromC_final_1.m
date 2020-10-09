% =========================================================================
%
%                  验证超高的高通滤波器的准确性
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
% 0929
%        4.角度的变化
%        5.
%
%--------------------------------------------------------------------------

clear all;
close all;
N = 20000;
start_pos = 1;

%% % wx,tbs,hfcra,t4,...
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
    tmp4 = tmp4(start_pos:start_pos+N-1,:);
end
lfcrp_ref = tmp4(:,3);

%% 应该除以一个系数，这里为什么没有？
s1 = 246083.654;
GYRO = GYRO / s1;
hfcra_func = Pz3_gj( GYRO , TBS );
figure;plot( hfcra_ref -  hfcra_func * s1);figure;plot( hfcra_ref );hold on;plot( hfcra_func * s1);%%不准


%% 对比角度
% hfcra_ref = out;
gpxbr = lfcrp_ref + hfcra_ref;
gpro = tmp(:,1);
s1 = 246083.654;
gpro_ = gpro/s1;
TBS = tmp(:,2)/1e5;
roll = 0;
for i = 1:length(TBS)
    roll = roll + TBS(i).*gpro_(i);
    roll_save(i,1) = roll;
end

sita_c = gpxbr/1638.4/1.03725;
figure;plot(roll_save/pi*180);hold on;plot(sita_c);legend matlab积分 gj;
%%这里可以看到漂移比较严重，但是基本可以确定其变化的范围

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
    up = N1 + N2;
    t4 = up / dnom4;
    x_k_1 = x_k;
    x_k_1_dot = x_k_dot;
    x_k_2_dot2 = x_k_1_dot2;
%%
    y_k_3_dot = y_k_2_dot;%%y dot i-3
    y_k_2_dot = y_k_1_dot;%%y dot i-2
    
    %% 差分
    hfcra(i,1) = (t4 - y_k_1)/2;
    y_k_2 = y_k_1;%%y i-2
    y_k_1 = t4;%%y i-1
    %% 更新
    tbspp = tbs(1);
    tbs(1) = tbs(2);
end

%%
out = hfcra;
end



