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
%        4.
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


%%  �ҵ��˲����Ա�
for i =1:length(tmp)
    gpro = tmp(i,1);
    tbs = tmp(i,2);
%     gpro = B(gpro,tbs);
    gpro_Pz = P3z(gpro,tbs);
    out(i) = gpro_Pz;
end
figure;plot(out);hold on;figure;plot(frlp_save);




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
% T = tbs(2);
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
N2_2 = 2^37*(1/2+1/2^5+1/2^9) * tbs(2)*( tbs(1) * y_k_1_dot + tbs(2)*y_k_1 );
N2_3 = 2^51*(1/2+1/2^5) * ( tbs(1)* ( y_k_1_dot + y_k_1_dot - y_k_2_dot ) + tbs(2)*y_k_1 );
N2_4 = 2^64*( 3*( y_k_1_dot - y_k_2_dot ) + y_k_3_dot + y_k_1 );
N2 = 2^(-64) * ( N2_1 + N2_2 + N2_3 + N2_4 );

up = 2^64*( N1/Omega1 + N2 );

down = 2^64 + tbs(2) *2^51 *(1/2+1/2^5) + tbs(2)^2*2^51*(1/2+1/2^5+1/2^9) + tbs(2)^3*2^19*(1/2+1/2^5) + tbs(2)^4;
    

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








