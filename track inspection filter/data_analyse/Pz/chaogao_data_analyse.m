% =========================================================================
%
%                  ��֤�˲�����׼ȷ��
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��18��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.Ps���������⣬��ֵ̫С
%        2.��׼�������������������ڳ�ʼ���Ĺ����г��ֵ����������ڣ�
%               ��ʼ�����Ĳ���
%        3.�˲���������ĳ�ʼ����ô�죿
%        4.�����Ƕȴ�С�����ֵ�λʼ����һ�����⣬����׼ȷ�Ե���
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

%% ��ȡyaw rate
tmp3 = textread('fmctrl_data_1337.txt');
if length(tmp3)>N
    tmp3 = tmp3(1:N,:);
end
gpvlo = tmp3(:,8);
gpvro = tmp3(:,10);
% ��excel���

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
%% ��ʼ��

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
    xtbsp = tmp(i-1,2);%%��һ���ͺ���Ҫ
    
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

%% �Ƕȷ�����from gj��

gpxbr = tmp2(:,1);      %%sita_b
hfcra = tmp2(:,2);      %%��Ƶ
lfcrp = tmp2(:,3);      %%��Ƶ

dt74 = (gpvro - gpvlo)/2;

scalr = 0.5;
gpxl = gpxbr + dt74*scalr;      %%��������յĽǶȣ���ô�Ƕ��Ƕ��٣�

gpxl = gpxl/32.2520/100;        %%��λ���ĵ���һ��
gpxbr = gpxbr/32.2520/100;      %%��λΪ��
figure;plot(gpxl);

%% ���������ݷ��� wx
wx = tmp(:,1)/3276.8/180*pi;

% ����
sita = 0;

for i = 1:length(wx)
    
    tbs = tmp(i,2)/1e5;
    
    sita = sita + wx(i)*tbs;
    
    sita_save(i,1) = sita;
end

%%�ǶȶԱ�
figure;plot(gpxbr)
hold on;plot(sita_save/pi*180);
legend gj ֱ�ӻ���

% figure;plot(gyro/pi*180);legend 1 2 3%%�о�����ĵ�λ�Ǵ���ģ������ǽ��ٶȵĵ�λ

%%

sita_bt = dt74*scalr;
figure;plot(sita_bt/32.2520/100);%%sita_btֻ��0.1�ȣ�


%% ����
gyroll = tmp3(:,1)/3276.8/180*pi;
gypitch = tmp3(:,2)/3276.8/180*pi;
gyyaw = tmp3(:,6)/3276.8/180*pi;
accx = tmp3(:,3)/129.01;%%�ⵥλ������ʲô��
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
%%�����������ֱ�Ӷ�yaw�Ļ�����ʮ�ֲ�׼ȷ�ģ��൱���г���1km֮��ת��Ϊ70�ȣ����Ƿ�
%%�����׵ģ��������޸�

%% ���ٶȻ���
%% �ǶȻ�����Ȼ��������
gyro = [gyroll,gypitch,gyyaw];
angle = zeros(1,3);
for t = 2:length(gyro)
    tbs = tmp(i,2)/1e5;
    ang_k = angle(t-1,:);
    del = cau_w(ang_k(1),ang_k(2),ang_k(3))*gyro(t,:)';
    angle(t,:) = (ang_k + del.'*tbs);
end
figure;plot(angle(:,1)/pi*180);%%deg������ǻ��ִ��������⣿
hold on;
plot(gpxbr)
legend 1 2 3


%% ���ٶȻ���
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





%% ����
 % figure;plot(gyro);legend 1 2 3
function W = cau_w(fai,sita,psi)
W = [1,sin(fai)*tan(sita) , cos(fai)*tan(sita);
    0, cos(fai) ,       -sin(fai);
    0, sin(fai)/cos(sita) , cos(fai)/cos(sita)];
end





