




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
%  功能： 1.Hz没啥问题
%        2.Fz没问题
%        3. 顺便验证了超高的第一部分
%        4. 超高的第二部分
%        5. 低通滤波器已经差不多完成，虽然有一些不一致的地方（尤其是刚开始的那一段）
%               初步认为是与初始化相关
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
dt74 = gpvro/2-gpvlo/2;
fd74 = -dt74;
result = tmp(:,4);

%% 计算
%% 这些参量分别都是什么意思？
hight = 3.820;
fim = 0.820210;
dgyro = 6.56;
fim1 = 1.2192;
compf = -1;
mfd12 = -374.08;
scali = 1230.2646;
mytbp = [];
dincl = 2.320;

%% 经过打印，确定过后的参数
%% 这里的参数就是因为很多代码冗余，所以改起来很费劲
hight = 1.8;
scali = 1230.264648;
dincl = 1.32;
mfd12 = -14.344533;
for i = 2:length(result)
    % ---------------------------------
    tbs = tmp(i,3);
    frct(i,1) = filter_1_unknow(fd74(i), tbs);
    frt(i,1) = gpro(i) - frct(i);
    % ----------------------------------
    dtemp = fim * tbs * yaw(i) * compf;
    mytb = tbs * (yaw(i) - yaw(i-1));
    dtemp = dtemp - dgyro * tbs * (yaw(i) - yaw(i-1));
    if isempty(mytbp)
        mytbp = mytb;
    end
    
    dtemp = dtemp + fim1 * mfd12 *(mytb - mytbp );
    dtemp = dtemp + hight*tbs*(frt(i) - frt(i-1));
    dtemp = scali*(dtemp + dincl * tbs * (  frct(i) - frct(i-1) ));
    
    %%
    out2(i,1) = dtemp;
    dtmp_Fz = F_xiuzheng(dtemp,tbs);    %%这个好像有问题？
    out(i,1) = dtmp_Fz;
    
    %% next step
    
    infp = B(gpin(i),tbs);
    infp_save(i,1) = infp;
    inc = 0.5 * infp + dtmp_Fz;
    inc_save(i,1) = inc;
    lfcrp(i,1) = H3z(inc,tbs);
    
    %% update
    mytbp = mytb;%%就这一个数比较不好弄
end

%% 结果对比

figure;plot(out);hold on;plot(result);
legend myresult gjresult;

figure;plot(out - result);
figure;plot(lfcrp_comp - lfcrp);
figure;plot(inc_comp - inc_save);%%这两个有区别，所以后面的就有区别

figure;plot(dt0-out2);

%% 超高相关函数
function out = filter_1_unknow(x_k,tbs)
persistent x y;
if isempty(x)
    x = zeros(2,1);
    y = zeros(2,1);
end
sd74 = 82281.94545;
sd74 = 126376.6875;     %%这才是正解
x(2) = x_k;
x_dot = x(2) - x(1);
%%
y(2) = ( x_dot*sd74 + 2^17 * y(1) )/(2^17 + tbs);

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



