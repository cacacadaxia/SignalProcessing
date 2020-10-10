% =========================================================================
%
%                       ���߳�������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��10��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ϳ��ߵĵ�Ƶ�͸�Ƶ�ļ��㣬���Ҹ����Լ�������޸�
%        2.
%        3. 
% ע�⣺
%       1.1.31��1.31����ֵһ�£���һ����Ҫע�⣬һ���������ǵĵ�λ����һ����omega1�ĵ���
%       2.int(gj programme) = ����*float(true)
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
%% ��ȡ����
tmp = textread('fz_filter_gaodi.txt');
if length(tmp)>N
    tmp = tmp(1:N,:);
end
gpvlo = tmp(:,1);%%32768/5 inch
gpvro = tmp(:,2);
dt74 = gpvro/2-gpvlo/2;%%Ϊʲô�ǳ���2 %%%%32768/10 inch
fd74 = -dt74;       %%�����
result = tmp(:,4);
%%

%% ����
%% ��Щ�����ֱ���ʲô��˼��
% hight = 3.820;
% fim = 0.820210;     %%0.25/0.3048 1feet����0.3048m
fim = 0.25;%%m
dgyro = 6.56 * 0.3048;%%m       %%���ݵ�λ��
compf = -1;
scali = 1230.2646;
mytbp = [];

DDISP = 59.5 * 0.3048;           %%DDISP*0.0254 = 1.5113 ����ɣ�
% sd74 = 126376.6875;     %%��������⣬����ɶ��˼��
% sd74 = 100000.0 * ( 75.194127 / DDISP );
% % SD74 IS 10**5 * (5*0.4*3.2768/L*SIN 5)
sd74 =  ( 75.194127 / DDISP );%%sd74 1e5����1����Ӱ�죬�����ս��Ӱ���С
%%sd74 = 246083.654/3276.8����һ����Ҫע�⣬����Ϊ����λ�ƼƵĲ��������Ľ��ٶȵ�λ���㵽�������ǵĽ��ٶȵ�λһ���õ�

%% ������ӡ��ȷ������Ĳ���
degti = 1.0375;
%% ��ת����Ϊm
hight = 1.8 * 0.3048;         %%feet   
% scali = 1/246083.654/0.54426*1638.4;%%���ֵ�ı仹����Ӱ���
dincl = 1.32 * 0.3048;       %%��ǼƵ�ת��ܵľ���
mfd12 = -(dgyro*2)*(dgyro*2)/12;%%��  --> feet
yaw_dotp = 0;

for i = 2:length(result)
    % ---------------------------------
    %%����dot�����֣�diff����΢��
    tbs = tmp(i,3);
    fd74_dot = (fd74(i) - fd74(i-1))*sd74;  %%--> ���� 246083.654
    fd74_diff = fd74_dot / ( tbs / 1e5 ); %%tbs/1e5����T
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
    dtemp = dtemp +  mfd12 * yaw_dot2 / fim / (tbs/1e5);%% ��ҪҲ�գ����׵������Խ��Ӱ���С
    dtemp = dtemp + hight  * wt_dot / (tbs/1e5);%% ht*s*F(s) %% wt_dot / tbs-->(wt_dot / tbs ����΢��diff)
%     dtemp = scali * ( dtemp + dincl * w_bt_dot /  (tbs/1e5) );%% -->( w_bt_dot / tbs ����΢��diff)
%     dtemp = (dtemp + dincl * w_bt_dot /  (tbs/1e5)) / 246083.654 /( 0.016932 * 9.8 ) * 3276.8 * 1.0057; %%1.0057��ǰ������scali�Ĳ�����ֵ
    dtemp = ( dtemp + dincl * w_bt_dot /  (tbs/1e5) ) / 246083.654 / 9.8 * 3276.8; %%1.0057��ǰ������scali�Ĳ�����ֵ
 
    %% sacli-->/246083�ǽ�������ת����Ϊdeg/s�ĵ�λ��
    %% 0.54426-->0.016932g/inch
    %%
    dtmp_Fz = F(dtemp , tbs);    %%������������⣿
    out(i,1) = dtmp_Fz;
    %% next step
    %     1/(0.016932*9.8) * 3276.8 = 19748
    infp = B2( gpin(i) / ( 3276.8/(0.016932*9.8) ) * 3276.8, tbs);   %%���B���ʹ�����Բ���Ӱ��  %%3276.8*rad��Ϊ����
    infp = infp / 9.8;      %%����g����
    infp_save(i,1) = infp;
    inc = infp + dtmp_Fz;     
    lfcrp(i,1) = H3z( inc , tbs );
    %%ת������
    lfcrp(i,1) = lfcrp(i,1) / 3276.8 * 1638.4 * 180 / pi; %%���ٶ�����1638.4*deg ����,��������еĽ�������
    %% update, last step
    yaw_dotp = yaw_dot;
    %%
    %%��dtemp�ദ����0.016932����gpin���ٴ�����0.016932��ͬʱ���ת����Ϊrad�Ļ�����ô������180/piת����Ϊdeg
    %%����ֵ���������ô�Ͷ����(0.016932*180/pi)==0.97������һ������
    
end
% lfcrp = lfcrp/1638.4/180*pi;
% lfcrp;%%����
figure;plot(lfcrp /(0.016932*180/pi) - lfcrp_comp);     %%��һ����Ҫע��һ�£���ԭ�����������
figure; plot(lfcrp / 1638.4);
%%ֱ�ӱ�ɵ�


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


%% Ӧ�ó���һ��ϵ��������Ϊʲôû�У�
s1 = 246083.654;
GYRO = GYRO  / 3276.8 / 1.31 * 1638.4; %%ת����Ϊ 3276.8*deg/s --> 1638.4*deg/s������֮�����ý����deg������������նԵ�����
hfcra_func = Pz3_gj( GYRO , TBS );
figure;plot( hfcra_ref -  hfcra_func );figure;plot( hfcra_ref );hold on;plot( hfcra_func );%%��׼

%%
reslut = hfcra_func + lfcrp / (0.016932*180/pi);%%����һ�����
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
%%  �������
degti = 1.037249;

%% ��ʼ����������Ҫ����
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
%% ��ʼ��
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
    x_k = rldd * degti;%%why��
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
    N2 = N2 + cx4*(tbs(1) * (y_k_1_dot + y_k_1_dot) - tbspp*y_k_2_dot + y_k_1 * tbs(2));    %%ע�������ftrp2
    N2 = N2 + cx2*tbs(2)*( tbs(1) * y_k_1_dot + y_k_1 * tbs(2) );
    N2 = N2 + cx5*tbs(2)*tbs(2)*y_k_1*tbs(2);
    N1 = N1 / 0.76;%%ע����һ��
    up = N1 + N2;
    t4 = up / dnom4;
    x_k_1 = x_k;
    x_k_1_dot = x_k_dot;
    x_k_2_dot2 = x_k_1_dot2;
%%
    y_k_3_dot = y_k_2_dot;%%y dot i-3
    y_k_2_dot = y_k_1_dot;%%y dot i-2
    
    %% ���
    hfcra(i,1) = (t4 - y_k_1);
    y_k_2 = y_k_1;%%y i-2
    y_k_1 = t4;%%y i-1
    %% ����
    tbspp = tbs(1);
    tbs(1) = tbs(2);
end

%%
out = hfcra;
end
