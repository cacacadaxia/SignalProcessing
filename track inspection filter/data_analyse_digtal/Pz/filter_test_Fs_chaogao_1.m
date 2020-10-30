




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
%        2.
%        3. 
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
sd74 = 126376.6875;     %%��������⣬����ɶ��˼��

%% ������ӡ��ȷ������Ĳ���
%% ����Ĳ���������Ϊ�ܶ�������࣬���Ը������ܷѾ�
hight = 1.8;         %%feet   
scali = 1230.264648;
dincl = 1.32;       %%��ǼƵ�ת��ܵľ���
mfd12 = -(dgyro*2)*(dgyro*2)/12;%%why�����Ǻܶ�
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
    % ����һ���ܸ��ӵ������߼����ǽ����е�����ת����Ϊfeet�ͺ���
    dtemp = fim  * yaw(i) * compf / tbs;%% wz -->(fim -- 0.25��� feet )
    dtemp = dtemp - dgyro * yaw_dot / tbs;%% L*s*F(s)-->(yaw_dot / tbs����΢��diff)
    dtemp = dtemp +  mfd12 * yaw_dot2 / fim / tbs;%%��ҪҲ��
    dtemp = dtemp + hight  * wt_dot / tbs;%% ht*s*F(s) %% wt_dot / tbs-->(wt_dot / tbs ����΢��diff)
    dtemp = scali * ( dtemp + dincl * w_bt_dot / tbs);%% -->( w_bt_dot / tbs ����΢��diff)
    %% scali������g
    %% scali���յ�λ���ʲô�ˣ�
    %%
    dtmp_Fz = F(dtemp,tbs);    %%������������⣿
    out(i,1) = dtmp_Fz;
    %% next step
    
    infp = B(gpin(i),tbs);
    infp_save(i,1) = infp;
    inc = 0.5 * infp + dtmp_Fz;
    lfcrp(i,1) = H3z(inc,tbs);
    
    %% update
    yaw_dotp = yaw_dot;
end



%% ����Ա�

figure;plot(out - result);
figure;plot(lfcrp_comp - lfcrp);


%% ������غ���
function out = filter_1_unknow(x_k,tbs)
%% ��ֺ󾭹�Bz
persistent x y;
if isempty(x)
    x = zeros(2,1);
    y = zeros(2,1);
end
sd74 = 82281.94545;
sd74 = 126376.6875;     %%��������⣬����ɶ��˼��
% 100000.0 * ( 75.194127 / DDISP );     DDISPλ�Ƽ�֮��ľ���
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



