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
%  功能： 1.Ps这里有问题，数值太小
%        2.不准的问题基本解决，但是在初始化的过程中出现的问题来自于：
%               初始化做的不好
%        3.滤波器的输出的初始化怎么办？不好弄
%        4.高频的部分，是对陀螺仪进行滤波加积分
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

tp1 = tmp(:,5);
tp2 = tmp(:,6);
tp3 = tmp(:,7);
tp4 = tmp(:,8);

middle = textread('Ps3_filter_wx_tmp4_err.txt');
if length(middle)>N
    middle = middle(1:N,:);
end
% not important

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
%%
degti = 1.037249;


%% 初始化的问题需要避免
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
%% 初始化

rldb_1 = tmp(1,1)*degti;
rldb_2 = tmp(2,1)*degti;
rldb_3 = tmp(3,1)*degti;
rldbp = rldb_3;
x1p3p = rldb_3 - rldb_2;
x1p2p = rldb_3 - 2*rldb_2 + rldb_1;

INILIZED=0;
for i = 4:length(tmp)
    
    rldd = tmp(i,1);
    ttbs = tmp(i,2);
    xtbs = ttbs;
    xtbsp = tmp(i-1,2);%%这一步就很重要
    
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
    
    tmp1(i) = dnom4;
    
    %% /* Multiply by 2 ** 62 */
    t3 = ct4x*t3;
    tmp2(i) = t3;
    
    
    x2p3 = x1p3*xtbs;
    
    %%/* Do another, Multiply by 2 ** 49 * (1/2 + 1/2 ** 3) */
    t3 = t3+cx7*(x2p3-x2p3p);
    
    
    t3 = t3+cx8*(xtbs*x2p3);
    
    
    t3 = t3+cx9*xtbs*xtbs*xtbs*rldb;
    tmp3(i) = t3;
    if (frlpp == 0)&&(INILIZED==0)
        frlpp = frlp;
 
    end
    frdp1=frlp-frlpp;
    fr3p = cx3*frdp1;
    
    if (fr3pp == 0)&&(INILIZED==0)
        fr3pp = fr3p;
    end
    t3 = t3+ct4x*(fr3p-fr3pp+frdp3+frlp);
    if (frdp3 == 0)&&(INILIZED==0)
        frdp2 = frdp1;
        frdp3 = frdp2;
    end
    ftrp=frlp*xtbs;
    ftrp1=xtbsp*frdp1;
    if (ftrp2 == 0)&&(INILIZED==0)
        ftrp2 = ftrp1;
        INILIZED = 1;
    end
    t3 = t3 + cx4*(ftrp1+ftrp1-ftrp2+ftrp);
    
    tmp4(i) = t3;
    tmp4_1(i,1) = ftrp1;
    tmp4_2(i,1) = ftrp2;
    tmp4_3(i,1) = ftrp;
    
    t3=t3+cx2*xtbs*(ftrp1+ftrp);
    
    t3=t3+cx5*xtbs*xtbs*ftrp;
    %%/* Result */
    t4 = t3/dnom4;
    t4_save(i,1) = t4;%%save
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
    
    %%  /* check for roll filter first pass case */
    if(p1flg==0)
        p1flg=-1;
        rollp=frlp;
    end
    
    t3 = frlp - rollp;
    rollp = frlp;
    
    t3 = t3*cx6;
   
    hfcra(i) = t3;
end


% 
figure;plot(tmp(:,3) - ceil(hfcra).');%%误差很小
figure;plot(tmp(:,3),'LineWidth',1);hold on;plot(hfcra)
figure;plot(tp1'-tmp1);axis([1,10000,-10,10]);
figure;plot(tp2'-tmp2);
figure;plot(tp3'-tmp3);
figure;plot(tp4'-tmp4);
figure;plot(t4_save-tmp(:,4));


figure;plot(tmp4_1 - middle(:,1));
figure;plot(tmp4_2 - middle(:,2));
figure;plot(tmp4_3 - middle(:,3));



