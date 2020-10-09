% =========================================================================
%
%                  ��֤���ߵĸ�ͨ�˲�����׼ȷ��
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
% 0929
%        4.����ʹ����д���˲������ڹ���һ��ʱ��֮��ͷ�ɢ�ˣ�����ԭ���������
%           �Ҳ�ȥ���ǣ�������Ŀǰ���Ĵ���
%        5.
%
%--------------------------------------------------------------------------

clear all;
close all;
N = 20000;
start_pos = 20000;
tmp = textread('Ps3_filter_wx.txt');
if length(tmp)>N
    tmp = tmp(start_pos : start_pos+N-1,:);
end
% wx,tbs,hfcra,t4,...

%%
tmp4 = textread('Hz_filter_inc.txt');
if length(tmp4)>N
    tmp4 = tmp4(start_pos:start_pos+N-1,:);
end
lfcrp_ref = tmp4(:,3);

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
rldb_1 = tmp(1,1)*degti;
rldb_2 = tmp(2,1)*degti;
rldb_3 = tmp(3,1)*degti;
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
for i = 4:length(tmp)
    %%
    rldd = tmp(i,1);
    tbs(2) = tmp(i,2);
    dnom4 = ( tbs(2)  + ct1 ) * tbs(2);
    dnom4 = ( dnom4 + ct2 ) * tbs(2);
    dnom4 = ( dnom4 + ct3 ) * tbs(2);
    dnom4 = dnom4 + ct4;
    x_k = rldd*degti;
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
    up = N1 + N2;
    t4 = up / dnom4;
    x_k_1 = x_k;
    x_k_1_dot = x_k_dot;
    x_k_2_dot2 = x_k_1_dot2;
    x2p3p = x2p3;
%%
    y_k_3_dot = y_k_2_dot;%%y dot i-3
    y_k_2_dot = y_k_1_dot;%%y dot i-2
    hfcra(i,1) = (t4 - y_k_1)/2;
    y_k_2 = y_k_1;%%y i-2
    y_k_1 = t4;%%y i-1
    frlp_save(i) = t4;
    
    %% ����
    tbspp = tbs(1);
    tbs(1) = tbs(2);
end

%%  �ҵ��˲����Ա�
for i =1:length(tmp)
    gpro = tmp(i,1)*degti;
    tbs = tmp(i,2);
    [gpro_Pz, tp_N1] = P3z(gpro,tbs);
    [gpro_Pz1, tp_N1] = P3z1(gpro,tbs);
    out(i,1) = gpro_Pz;
    out1(i,1) = gpro_Pz1;
end
hfcra_ref = tmp(:,3);
figure;plot( hfcra_ref -  hfcra );figure;plot(hfcra_ref);hold on;plot(hfcra);%%ֻ�������׼��
figure;plot( hfcra_ref -  out1  );figure;plot(hfcra_ref);hold on;plot(out1);%%��׼
figure;plot( hfcra_ref -  out   );figure;plot(hfcra_ref);hold on;plot(out);%%��׼

%% �ԱȽǶ�
hfcra_ref = tmp(:,3);
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
figure;plot(roll_save/pi*180);hold on;plot(sita_c);legend matlab���� gj;
%%������Կ���Ư�ƱȽ����أ����ǻ�������ȷ����仯�ķ�Χ



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

function [out,tmp] = P3z(x_k,tbs_k)

persistent x y tbs tbsp ;

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
    tbsp = 0;

end
x(5) = x_k;
tbs(2) = tbs_k;


T = tbs(2)/1e5;         %%����ע��
x_dot = x(5) - x(4);
x_dot2 = x(5) - 2*x(4) + x(3);
x_dot3 = x(5) - 3*x(4) + 3*x(3) - x(2);
x_dot4 = x(5) - 4*x(4) + 6*x(3) - 4*x(2) + x(1);
% N1 = kai*T^3*x_dot + nang*T^2*x_dot2 + zeta*T*x_dot3 + x_dot4;
N1 = kai*T^3*x_k + nang*T^2*x_dot + zeta*T*x_dot2 + x_dot3;
N1 = N1 * 2^64;

y_k_1_dot = y(4) - y(3);
y_k_2_dot = y(3) - y(2);
y_k_3_dot = y(2) - y(1);
y_k_1 = y(4);       %%n-1ʱ�̵�y
N2_1 = 2^19*(1/2+1/2^5)*tbs(2)^3*y_k_1;
N2_2 = 2^37*(1/2+1/2^5+1/2^9) * tbs(2)*( tbs(1) * y_k_1_dot + tbs(2)*y_k_1 );
N2_3 = 2^51*(1/2+1/2^5) * ( tbs(1)* ( y_k_1_dot + y_k_1_dot ) - y_k_2_dot * tbsp + tbs(2)*y_k_1 );
N2_4 = 2^64*( 3*( y_k_1_dot - y_k_2_dot ) + y_k_3_dot + y_k_1 );
N2 = ( N2_1 + N2_2 + N2_3 + N2_4 );

up = ( N1 + N2 );

down = 2^64 + tbs(2) *2^51 *(1/2+1/2^5) + tbs(2)^2*2^37*(1/2+1/2^5+1/2^9) + tbs(2)^3*2^19*(1/2+1/2^5) + tbs(2)^4;
tmp = down;
%%

ct1 = 278528.0;	 %%/* ct1=2**19*(1/2+1/2**5) */
ct2 = 7.3282879e10;  %%/* ct2=2**37*(1/2+1/2**5+1/2**9) */
ct3 = 1.1962686e15;  %%/* ct3=2**51*(1/2+1/2**5) */
ct4 = 1.8446744e19;  %%/* ct4=2**65*(1/2) */
    ttbs = tbs(2);
    dnom4 = ( ttbs  + ct1 ) * ttbs;
    dnom4 = ( dnom4 + ct2 ) * ttbs;
    dnom4 = ( dnom4 + ct3 ) * ttbs;
    dnom4 = dnom4 + ct4;
    down = dnom4;

y(5) = up/down;
out = y(5) - y(4);
out = out/2;
%% ����
x(1) = x(2);
x(2) = x(3);
x(3) = x(4);
x(4) = x(5);

y(1) = y(2);
y(2) = y(3);
y(3) = y(4);
y(4) = y(5);

tbsp = tbs(1);
tbs(1) = tbs(2);
end


function [out,tmp] = P3z1(x_k,tbs_k)

Omega1 = 10^5/2^17;
Omega2 = 10^5/2^14;
Omega3 = 10^5/2^18;

zeta = Omega1+Omega2+Omega3;
nang = Omega2^2 + Omega2*Omega3 + Omega2*Omega1 + Omega1*Omega3 +Omega3^2;
kai = (Omega1 + Omega3)*Omega2^2 + Omega2*Omega1*Omega3 + (Omega1 + Omega2)*Omega3^2;

persistent x y tbs y_k_1 y_k_1_dot y_k_2_dot y_k_3_dot y_k_2 yy tbsp; 
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
if isempty(x)
    x = zeros(5,1);
    y = zeros(5,1);
    tbs = zeros(2,1);
    y_k_1_dot = 0;
    y_k_2_dot = 0;
    y_k_3_dot = 0;
    y_k_1 = 0;
    y_k_2 = 0;
    yy = 0;
    tbsp = 0;
end
%%
tbs(2) = tbs_k;
x(5) = x_k;

%%
T = tbs(2)/1e5;         %%����ע��
x_dot = x(5) - x(4);
x_dot2 = x(5) - 2*x(4) + x(3);
x_dot3 = x(5) - 3*x(4) + 3*x(3) - x(2);
x_dot4 = x(5) - 4*x(4) + 6*x(3) - 4*x(2) + x(1);
% N1 = kai*T^3*x_dot + nang*T^2*x_dot2 + zeta*T*x_dot3 + x_dot4;
N1 = kai*T^3*x_k + nang*T^2*x_dot + zeta*T*x_dot2 + x_dot3;
N1 = N1 * 2^64;


%%
y_k_1_dot = y_k_1 - y_k_2;
N2 = 0 + ct4x*(3*(y_k_1_dot - y_k_2_dot) + y_k_3_dot + y_k_1);
N2 = N2 + cx4*(tbs(1) * (y_k_1_dot + y_k_1_dot) - tbsp*y_k_2_dot + y_k_1 * tbs(2));    %%ע�������ftrp2
N2 = N2 + cx2*tbs(2)*( tbs(1) * y_k_1_dot + y_k_1 * tbs(2) );
N2 = N2 + cx5*tbs(2)*tbs(2)*y_k_1*tbs(2);
up = N1 + N2;

down = 2^64 + tbs(2) *2^51 *(1/2+1/2^5) + tbs(2)^2*2^37*(1/2+1/2^5+1/2^9) + tbs(2)^3*2^19*(1/2+1/2^5) + tbs(2)^4;
t4 = up / down;
out = (t4 - y_k_1)/2;
tmp = 0;
%%

x(1) = x(2);
x(2) = x(3);
x(3) = x(4);
x(4) = x(5);
tbsp = tbs(1);
tbs(1) = tbs(2);
y_k_3_dot = y_k_2_dot;%%y dot i-3
y_k_2_dot = y_k_1_dot;%%y dot i-2

y_k_2 = y_k_1;%%y i-2
y_k_1 = t4;%%y i-1
end




