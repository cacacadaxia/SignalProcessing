


% =========================================================================
%
%                   长波滤波器设计
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：
%   作者：du
%--------------------------------------------------------------------------
%  功能： 1.三角窗的表示可能会有点问题
%        2.
%        3. 
%--------------------------------------------------------------------------



close all
clear all
%************重要**************du
% 三角窗+矩形窗并联；梯形窗并联；矩形窗级联 比较

% load para.txt
% load para2.txt
% para = para2;
%
N=128^2;
w=linspace(0,pi,N);
lamda = 2*pi./w./4;
lam = [];


% 65是什么？

%三角窗+矩形窗并联
[win25L,win25M,win25N] = winlen(65);
[win42L,win42M,win42N] = winlen(111);
[win70L,win70M,win70N] = winlen(189);
[win120L,win120M,win120N] = winlen(319);
% 都是FIR，所以参数就直接组成数组就可以了
para=[65, [win25L,win25M,win25N];
      111, [win42L,win42M,win42N] ;
      189,[win70L,win70M,win70N];
      319, [win120L,win120M,win120N]];

% 25、42、70分别是波长
[Acc_TriRec_25, TriDB_TriRec_25, maxV_TriRec_25, minV_TriRec_25, BW_TriRec_25, slope_TriRec_25] = tri_rec_paral_my(25, para);
[Acc_TriRec_42, TriDB_TriRec_42, maxV_TriRec_42, minV_TriRec_42, BW_TriRec_42, slope_TriRec_42] = tri_rec_paral_my(42, para);
[Acc_TriRec_70, TriDB_TriRec_70, maxV_TriRec_70, minV_TriRec_70, BW_TriRec_70, slope_TriRec_70] = tri_rec_paral_my(70, para);
[Acc_TriRec_120, TriDB_TriRec_120, maxV_TriRec_120, minV_TriRec_120, BW_TriRec_120, slope_TriRec_120] = tri_rec_paral_my(120, para);

%梯形窗并联
[Acc_Trap_25, TriDB_Trap_25, maxV_Trap_25, minV_Trap_25, BW_Trap_25, slope_Trap_25] = tri_paral(25);
[Acc_Trap_42, TriDB_Trap_42, maxV_Trap_42, minV_Trap_42, BW_Trap_42, slope_Trap_42] = tri_paral(42);
[Acc_Trap_70, TriDB_Trap_70, maxV_Trap_70, minV_Trap_70, BW_Trap_70, slope_Trap_70] = tri_paral(70);
[Acc_Trap_120, TriDB_Trap_120, maxV_Trap_120, minV_Trap_120, BW_Trap_120, slope_Trap_120] = tri_paral(120);

%4矩形窗级联
cas = 4;
[Acc_4Rec_25, TriDB_4Rec_25, maxV_4Rec_25, minV_4Rec_25, BW_4Rec_25, slope_4Rec_25] = rec_casc(25,cas);
[Acc_4Rec_42, TriDB_4Rec_42, maxV_4Rec_42, minV_4Rec_42, BW_4Rec_42, slope_4Rec_42] = rec_casc(42,cas);
[Acc_4Rec_70, TriDB_4Rec_70, maxV_4Rec_70, minV_4Rec_70, BW_4Rec_70, slope_4Rec_70] = rec_casc(70,cas);
[Acc_4Rec_120, TriDB_4Rec_120, maxV_4Rec_120, minV_4Rec_120, BW_4Rec_120, slope_4Rec_120] = rec_casc(120,cas);

cas = 3;
[Acc_3Rec_25, TriDB_3Rec_25, maxV_3Rec_25, minV_3Rec_25, BW_3Rec_25, slope_3Rec_25] = rec_casc(25,cas);
[Acc_3Rec_42, TriDB_3Rec_42, maxV_3Rec_42, minV_3Rec_42, BW_3Rec_42, slope_3Rec_42] = rec_casc(42,cas);
[Acc_3Rec_70, TriDB_3Rec_70, maxV_3Rec_70, minV_3Rec_70, BW_3Rec_70, slope_3Rec_70] = rec_casc(70,cas);
[Acc_3Rec_120, TriDB_3Rec_120, maxV_3Rec_120, minV_3Rec_120, BW_3Rec_120, slope_3Rec_120] = rec_casc(120,cas);

cas = 2;
[Acc_2Rec_25, TriDB_2Rec_25, maxV_2Rec_25, minV_2Rec_25, BW_2Rec_25, slope_2Rec_25] = rec_casc(25,cas);
[Acc_2Rec_42, TriDB_2Rec_42, maxV_2Rec_42, minV_2Rec_42, BW_2Rec_42, slope_2Rec_42] = rec_casc(42,cas);
[Acc_2Rec_70, TriDB_2Rec_70, maxV_2Rec_70, minV_2Rec_70, BW_2Rec_70, slope_2Rec_70] = rec_casc(70,cas);
[Acc_2Rec_120, TriDB_2Rec_120, maxV_2Rec_120, minV_2Rec_120, BW_2Rec_120, slope_2Rec_120] = rec_casc(120,cas);

cas = 1;
[Acc_1Rec_25, TriDB_1Rec_25, maxV_1Rec_25, minV_1Rec_25, BW_1Rec_25, slope_1Rec_25] = rec_casc(25,cas);
[Acc_1Rec_42, TriDB_1Rec_42, maxV_1Rec_42, minV_1Rec_42, BW_1Rec_42, slope_1Rec_42] = rec_casc(42,cas);
[Acc_1Rec_70, TriDB_1Rec_70, maxV_1Rec_70, minV_1Rec_70, BW_1Rec_70, slope_1Rec_70] = rec_casc(70,cas);
[Acc_1Rec_120, TriDB_1Rec_120, maxV_1Rec_120, minV_1Rec_120, BW_1Rec_120, slope_1Rec_120] = rec_casc(120,cas);

Title = strcat('Amplitude-Frequency Curve--- Wavelength: ',num2str(lam),' m ');
figure;
semilogx((lamda), ((Acc_TriRec_25)), 'b','LineWidth',2);
hold on
semilogx((lamda), ((Acc_TriRec_42)), 'g','LineWidth',2);
semilogx((lamda), ((Acc_TriRec_70)), 'r','LineWidth',2);
semilogx((lamda), ((Acc_TriRec_120)), 'k','LineWidth',2);

semilogx((lamda), (abs(Acc_Trap_25)), '--b','LineWidth',2);
semilogx((lamda), (abs(Acc_Trap_42)), '--g','LineWidth',2);
% semilogx((lamda), (abs(Acc_Trap_70)), '--r','LineWidth',2);
semilogx((lamda), (abs(Acc_Trap_120)), '--k','LineWidth',2);

% semilogx((lamda), ((Acc_4Rec_25)), ':b','LineWidth',2);
% semilogx((lamda), ((Acc_4Rec_42)), ':g','LineWidth',2);
% semilogx((lamda), ((Acc_4Rec_70)), ':r','LineWidth',2);
semilogx((lamda), ((Acc_4Rec_120)), ':k','LineWidth',2);

% legend 三角窗+矩形窗 梯形窗并联 矩形窗并联