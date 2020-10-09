




% =========================================================================
%
%                  ��֤�˲�����׼ȷ��
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��28��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ߵĶԱȣ�ģ�͵ļ�
%        2.�����Լ�������޸����еĲ��������ߵ�������̬
%        3. �ǶȵĻ���
%  ע�⣺
%        1.ֻҪ�Ǻͽ��ٶ���صĶ���Ҫע��comp
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
%% ��ȡdtemp
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
%% ��ȡyaw rate
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
%% ��������
%% ���ݸ�ʽ��ʲô��
tmp = textread('fz_filter_gaodi.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end
gpvlo = tmp(:,1);
gpvro = tmp(:,2);
dt74 = gpvro/2-gpvlo/2;%%Ϊʲô�ǳ���2
fd74 = -dt74;
result = tmp(:,4);

%% ����
%% ��Щ�����ֱ���ʲô��˼��
hight = 3.820;
fim = 0.820210;     %%0.25/0.3048 1feet����0.3048m
dgyro = 6.56;       %%���ݵ�λ��
compf = -1;
scali = 1230.2646;
mytbp = [];

DDISP = 59.5;           %%DDISP*0.0254 = 1.5113 ����ɣ�
sd74 = 126376.6875;     %%��������⣬����ɶ��˼��
sd74 = 100000.0 * ( 75.194127 / DDISP );
% SD74 IS 10**5 * (5*0.4*3.2768/L*SIN 5)
sd74 =  ( 75.194127 / DDISP );%%sd74 1e5����1����Ӱ��

%% ������ӡ��ȷ������Ĳ���
%% ����Ĳ���������Ϊ�ܶ�������࣬���Ը������ܷѾ�
hight = 1.8;         %%feet   
scali = 1230.264648 / 1e5;
% scali = 1/246083.654/0.54426*1638.4;%%���ֵ�ı仹����Ӱ���
dincl = 1.32;       %%��ǼƵ�ת��ܵľ���
% mfd12 = -(dgyro*2)*(dgyro*2)/12;%%why�����Ǻܶ������ֵ��Ӣ���й�ϵ��  --> feet(���feet)
yaw_dotp = 0;
for i = 2:length(result)
    % ---------------------------------
    %%����dot�����֣�diff����΢��
    tbs = tmp(i,3);
    fd74_dot = (fd74(i) - fd74(i-1))*sd74;
    fd74_diff = fd74_dot / ( tbs/1e5 ); %%tbs/1e5����T
    w_bt(i,1) = B1(fd74_diff,tbs);
    w_bt_dot = w_bt(i,1) - w_bt(i-1,1);
    wt(i,1) = gpro(i) - w_bt(i);        %%����������λͳһ��Ӧ����
    wt_dot = wt(i,1) - wt(i-1,1);
    yaw_dot = yaw(i) - yaw(i-1);
    yaw_dot2 = yaw_dot - yaw_dotp;
    % ----------------------------------
    % ����һ���ܸ��ӵ������߼����ǽ����е�����ת����Ϊfeet�ͺ���
    dtemp = fim  * yaw(i) * compf / (tbs/1e5);%% wz -->(fim -- 0.25��� feet ��fim/tbs�����ٶ�)
    dtemp = dtemp - dgyro * yaw_dot / (tbs/1e5);%% L*s*F(s)-->(yaw_dot / tbs����΢��diff)
%     dtemp = dtemp +  mfd12 * yaw_dot2 / fim / (tbs/1e5);%%��ҪҲ�գ����׵���
    dtemp = dtemp + hight  * wt_dot / (tbs/1e5);%% ht*s*F(s) %% wt_dot / tbs-->(wt_dot / tbs ����΢��diff)
    dtemp = scali * ( dtemp + dincl * w_bt_dot /  (tbs/1e5) );%% -->( w_bt_dot / tbs ����΢��diff)
    % Ӧ����tbs/1e5=T������΢����Ҫע�����������scali�а���
    
    %% scali������g
    %% sacli-->/246083�ǽ�������ת����Ϊ����/s�ĵ�λ��
    %% 0.54426-->0.016932g/inch
    %% scali ���յ�λ���ʲô�ˣ�
    %%
    dtmp_Fz = F(dtemp , tbs);    %%������������⣿
    out(i,1) = dtmp_Fz;
    %% next step
    s1 = 19748;                     %%ϵ�������ս�������ж���Ӱ��
    infp = B2(gpin(i)/s1,tbs)*s1;   %%���B���ʹ�����Բ���Ӱ��
    infp_save(i,1) = infp;
    inc = 0.5 * infp + dtmp_Fz;     %%Ϊʲô����2����һ�����˺�����
    lfcrp(i,1) = H3z(inc,tbs);
    
    %% update
    yaw_dotp = yaw_dot;
end


%%
% gpxlp0 = (int)(gpxbr+(long)dt74*scalr)/2;
% scalr= 1638.4/(3276.8*DDISP/REFL);//����λ�Ƽ�֮��ľ���

%% ֱ�ӻ���

s1 = 246083.654;%%���������ģ�
gpro_ = gpro/s1;
TBS = tmp(:,3)/1e5;
roll = 0;
for i = 1:length(TBS)
    roll = roll + TBS(i).*gpro_(i);
    roll_save(i,1) = roll;
end
figure;plot(roll_save/pi*180);
sita_c = lfcrp_comp/1638.4/1.03725;%%���ϵ������Ҫ���൱��ͳһ�ĵ�λ
hold on;plot(sita_c);
legend matlab 

%% ����Ա�
figure;plot(result - out);figure;plot(result);hold on;plot(out);
figure;plot(lfcrp);hold on;plot(lfcrp_comp);legend matlab gj;
figure;plot(lfcrp_comp - lfcrp);


%% ������غ���
function out = filter_1_unknow(x_k,tbs)
%% ��ֺ󾭹�Bz
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
% ���F�˲������˵�ʲô����������x*tbs^2����һ����
%% �˲����趨
%% ����x[3]������ͬ��
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2

%% ����˲���������ʱ����ʮ���ȣ�
%% �������õģ�
persistent y;
if isempty(y)
    y = zeros(3,1);
end

y(3) = ( y(2)*(2*2^28 + 2^14*tbs) - y(1)*2^28 +  x  )/(2^28 + 2^14*tbs + tbs^2);

%% ����temp��
y(1) = y(2);
y(2) = y(3);
out = y(3);

end

function out = F(x,tbs)
%% �˲����趨
%% ����x[3]������ͬ��
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2

%% ����˲���������ʱ����ʮ���ȣ�
%% �������õģ�
persistent y;
if isempty(y)
    y = zeros(3,1);
end
y(3) = ( y(2)*(2*2^28 + 2^14*tbs) - y(1)*2^28 + tbs^2 * x  )/(2^28 + 2^14*tbs + tbs^2);
%% ����temp��
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

%% ����
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


T = tbs(2)/1e5;         %%����ע��
x_dot = x(5) - x(4);
x_dot2 = x(5) - 2*x(4) + x(3);
x_dot3 = x(5) - 3*x(4) + 3*x(3) - x(4);
x_dot4 = x(5) - 4*x(4) + 6*x(3) - 4*x(2) + x(1);
N1 = kai*T^3*x_dot + nang*T^2*x_dot2 + zeta*T*x_dot3 + x_dot4;

y_k_1_dot = y(4) - y(3);
y_k_2_dot = y(3) - y(2);
y_k_3_dot = y(2) - y(1);
y_k_1 = y(4);       %%n-1ʱ�̵�y
N2_1 = 2^19*(1/2+1/2^5)*tbs(2)^3*y_k_1;
N2_2 = 2^37*(1/2+1/2^5+1/2^9) *tbs(2)*( tbs(1) * y_k_1_dot + tbs(2)*y_k_1 );
N2_3 = 2^51*(1/2+1/2^5) * ( tbs(1)* ( y_k_1_dot + y_k_1_dot - y_k_2_dot ) + tbs(2)*y_k_1 );
N2_4 = 2^64*( 3*( y_k_1_dot - y_k_2_dot ) + y_k_3_dot + y_k_1 );
N2 = 2^(-64) * ( N2_1 + N2_2 + N2_3 + N2_4 );

up = 2^64*( N1/Omega1 + N2 );
down = 2^64 + tbs(2) *2^51 *(1/2+1/2^5) + tbs(2)^2*2^51*(1/2+1/2^5+1/2^9) + tbs(2)^3*2^19*(1/2+1/2^5) + tbs(2)^4;

y(5) = up/down;
%% ����
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
%% ����kʱ��
x(3) = x_k;
tbs(2) = tbs_k;

%% ��ɢ�˲�
x_dot = x(3)-x(2);
x_dot_2 = x(3) - 2*x(2) + x(1);
y1(2) = (x(3) + x(2)) *tbs(2)/2^15 + x_dot;
y = (y1(2)+y1(1))* (tbs(2) + tbs(1))/2^16 + x_dot_2;

%% ����temp��
y1(1) = y1(2);
x(1) = x(2);
x(2) = x(3);
tbs(1) = tbs(2);
out = y;

end



