




% =========================================================================
%
%                  验证滤波器的准确性
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月28日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.超高的对比，模型的简化
%        2.按照自己的理解修改其中的参量，超高的最终形态
%        3. 角度的换算
%  注意：
%        1.只要是和角速度相关的都需要注意comp
% 
%--------------------------------------------------------------------------



clear all;
close all;

%%
N = 10000;
tmp2 = textread('Ps3_filter_wx.txt');
if length(tmp2)>N
    tmp2 = tmp2(1:N,:);
end
gpro = tmp2(:,1);
%% 读取dtemp
dtemp_comp = textread('tmp_dtemp.txt');
if length(dtemp_comp)>N
    dtemp_comp = dtemp_comp(1:N,:);
end

dt0 = dtemp_comp(:,1);
dt1 = dtemp_comp(:,2);
dt2 = dtemp_comp(:,3);
dt3 = dtemp_comp(:,4);
dt4 = dtemp_comp(:,5);
dt = dtemp_comp(:,6);
frt_ = dtemp_comp(:,7);
frct_ = dtemp_comp(:,8);
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
%% 导入数据
%% 数据格式是什么？
tmp = textread('fz_filter_gaodi.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end
gpvlo = tmp(:,1);
gpvro = tmp(:,2);
dt74 = gpvro/2-gpvlo/2;%%为什么是除以2
fd74 = -dt74;
result = tmp(:,4);

%% 计算
%% 这些参量分别都是什么意思？
hight = 3.820;
fim = 0.820210;     %%0.25/0.3048 1feet等于0.3048m
dgyro = 6.56;       %%陀螺的位置
compf = -1;
scali = 1230.2646;
mytbp = [];

DDISP = 59.5;           %%DDISP*0.0254 = 1.5113 这个吧？
sd74 = 126376.6875;     %%这才是正解，这是啥意思？
sd74 = 100000.0 * ( 75.194127 / DDISP );
% SD74 IS 10**5 * (5*0.4*3.2768/L*SIN 5)
sd74 =  ( 75.194127 / DDISP );%%sd74 1e5或者1都不影响

%% 经过打印，确定过后的参数
%% 这里的参数就是因为很多代码冗余，所以改起来很费劲
hight = 1.8;         %%feet   
scali = 1230.264648 / 1e5;
% scali = 1/246083.654/0.54426*1638.4;%%这个值改变还是有影响的
dincl = 1.32;       %%倾角计到转向架的距离
% mfd12 = -(dgyro*2)*(dgyro*2)/12;%%why？不是很懂，这个值与英寸有关系的  --> feet(变成feet)
yaw_dotp = 0;
for i = 2:length(result)
    % ---------------------------------
    %%其中dot代表差分，diff代表微分
    tbs = tmp(i,3);
    fd74_dot = (fd74(i) - fd74(i-1))*sd74;
    fd74_diff = fd74_dot / ( tbs/1e5 ); %%tbs/1e5代表T
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
%     dtemp = dtemp +  mfd12 * yaw_dot2 / fim / (tbs/1e5);%%不要也罢，二阶的量
    dtemp = dtemp + hight  * wt_dot / (tbs/1e5);%% ht*s*F(s) %% wt_dot / tbs-->(wt_dot / tbs 代表微分diff)
    dtemp = scali * ( dtemp + dincl * w_bt_dot /  (tbs/1e5) );%% -->( w_bt_dot / tbs 代表微分diff)
    % 应该是tbs/1e5=T，计算微分需要注意这个量，在scali中包含
    
    %% scali包含了g
    %% sacli-->/246083是将陀螺仪转换成为弧度/s的单位？
    %% 0.54426-->0.016932g/inch
    %% scali 最终单位变成什么了？
    %%
    dtmp_Fz = F(dtemp , tbs);    %%这个好像有问题？
    out(i,1) = dtmp_Fz;
    %% next step
    s1 = 19748;                     %%系数对最终结果不会有多少影响
    infp = B2(gpin(i)/s1,tbs)*s1;   %%这个B多次使用所以产生影响
    infp_save(i,1) = infp;
    inc = 0.5 * infp + dtmp_Fz;     %%为什么除以2？这一点让人很困惑
    lfcrp(i,1) = H3z(inc,tbs);
    
    %% update
    yaw_dotp = yaw_dot;
end


%%
% gpxlp0 = (int)(gpxbr+(long)dt74*scalr)/2;
% scalr= 1638.4/(3276.8*DDISP/REFL);//两个位移计之间的距离

%% 直接积分

s1 = 246083.654;%%这是那来的？
gpro_ = gpro/s1;
TBS = tmp(:,3)/1e5;
roll = 0;
for i = 1:length(TBS)
    roll = roll + TBS(i).*gpro_(i);
    roll_save(i,1) = roll;
end
figure;plot(roll_save/pi*180);
sita_c = lfcrp_comp/1638.4/1.03725;%%这个系数很重要，相当于统一的单位
hold on;plot(sita_c);
legend matlab 

%% 结果对比
figure;plot(result - out);figure;plot(result);hold on;plot(out);
figure;plot(lfcrp);hold on;plot(lfcrp_comp);legend matlab gj;
figure;plot(lfcrp_comp - lfcrp);


%% 超高相关函数
function out = filter_1_unknow(x_k,tbs)
%% 差分后经过Bz
persistent  y;
if isempty(y)
    y = zeros(2,1);
end
%%
y(2) = ( x_k + 2^17 * y(1) )/(2^17 + tbs);
% x_k
%% update
y(1) = y(2);
out = y(2);
end


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


function out = F_xiuzheng(x,tbs)
% 这个F滤波器少了点什么？就是少了x*tbs^2，这一个项
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

y(3) = ( y(2)*(2*2^28 + 2^14*tbs) - y(1)*2^28 +  x  )/(2^28 + 2^14*tbs + tbs^2);

%% 更新temp量
y(1) = y(2);
y(2) = y(3);
out = y(3);

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



