




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
%        2.
%        3. 
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
sd74 = 126376.6875;     %%这才是正解，这是啥意思？

%% 经过打印，确定过后的参数
%% 这里的参数就是因为很多代码冗余，所以改起来很费劲
hight = 1.8;         %%feet   
scali = 1230.264648;
dincl = 1.32;       %%倾角计到转向架的距离
mfd12 = -(dgyro*2)*(dgyro*2)/12;%%why？不是很懂
yaw_dotp = 0;
for i = 2:length(result)
    % ---------------------------------
    tbs = tmp(i,3);
    fd74_dot = (fd74(i) - fd74(i-1))*sd74;
    fd74_diff = fd74_dot/tbs;
    w_bt(i,1) = filter_1_unknow(fd74_dot , tbs);
%     w_bt(i,1) = B(fd74_diff*sd74,tbs);
    w_bt_dot = w_bt(i,1) - w_bt(i-1,1);
    wt(i,1) = gpro(i) - w_bt(i);
    wt_dot = wt(i,1) - wt(i-1,1);
    yaw_dot = yaw(i) - yaw(i-1);
    yaw_dot2 = yaw_dot - yaw_dotp;
    % ----------------------------------
    % 这是一个很复杂的量，逻辑就是将所有的量都转换成为feet就好了
    dtemp = fim  * yaw(i) * compf / tbs;%% wz -->(fim -- 0.25变成 feet )
    dtemp = dtemp - dgyro * yaw_dot / tbs;%% L*s*F(s)-->(yaw_dot / tbs代表微分diff)
%     dtemp = dtemp +  mfd12 * yaw_dot2 / fim / tbs;%%不要也罢
    dtemp = dtemp + hight  * wt_dot / tbs;%% ht*s*F(s) %% wt_dot / tbs-->(wt_dot / tbs 代表微分diff)
    dtemp = scali * ( dtemp + dincl * w_bt_dot / tbs);%% -->( w_bt_dot / tbs 代表微分diff)
    
    %%
    dtmp_Fz = F(dtemp,tbs);    %%这个好像有问题？
    
    %% next step
    
    infp = B(gpin(i),tbs);
    infp_save(i,1) = infp;
    inc = 0.5 * infp + dtmp_Fz;
    lfcrp(i,1) = H3z(inc,tbs);
    
    %% update
    yaw_dotp = yaw_dot;
end

%%
% %% 新的超高
% for i = 2:length(result)
% 
%     tbs = tmp(i,3)/1e5;
%     
% end



%% 结果对比

% figure;plot(out);hold on;plot(result);legend matlab gj;
figure;plot(lfcrp_comp - lfcrp);


%% 超高相关函数
function out = filter_1_unknow(x_k,tbs)
%% 差分后经过Bz
persistent x y;
if isempty(x)
    x = zeros(2,1);
    y = zeros(2,1);
end
sd74 = 82281.94545;
sd74 = 126376.6875;     %%这才是正解，这是啥意思？
% 100000.0 * ( 75.194127 / DDISP );     DDISP位移计之间的距离
x(2) = x_k;
x_dot = x(2) - x(1);
%%
y(2) = ( x_k + 2^17 * y(1) )/(2^17 + tbs);

%% update
x(1) = x(2);
y(1) = y(2);
out = y(2);
end


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



