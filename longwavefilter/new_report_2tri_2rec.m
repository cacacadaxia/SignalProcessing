% -----------------------
% 
%   1. 对修改的系数进行检验
% 
% ------------------------

clear all;
close all;
% 横坐标配置
N=128^2;
w=linspace(0,pi,N);
lamda = 2*pi./w./4;
% 参数配置
para = winLen_para_set();
% para(3,:) = [160*2,40*2,169,280];%%word文档中的滤波器参数
% 获取滤波器参数
[Acc_TriRec_25, TriDB_TriRec_25, maxV_TriRec_25, minV_TriRec_25, BW_TriRec_25, slope_TriRec_25] = tri_rec_paral_repair(25, para);
[Acc_TriRec_42, TriDB_TriRec_42, maxV_TriRec_42, minV_TriRec_42, BW_TriRec_42, slope_TriRec_42] = tri_rec_paral_repair(42, para);
[Acc_TriRec_70, TriDB_TriRec_70, maxV_TriRec_70, minV_TriRec_70, BW_TriRec_70, slope_TriRec_70] = tri_rec_paral_repair(70, para);
[Acc_TriRec_120, TriDB_TriRec_120, maxV_TriRec_120, minV_TriRec_120, BW_TriRec_120, slope_TriRec_120] = tri_rec_paral_repair(120, para);

% plot 25
fig = figure;
% semilogx((lamda), ((Acc_TriRec_25)), 'b','LineWidth',2);
% plot_line_func(Acc_TriRec_25,lamda,25);
% xlabel('Wavelength /m')
% ylabel('Magnitude /dB')
% set(gca,'Fontname','Times New Roman','fontsize',14);
% % plot 42
% % figure;
% semilogx((lamda), ((Acc_TriRec_42)), 'b','LineWidth',2);
% plot_line_func(Acc_TriRec_42,lamda,42);
% xlabel('Wavelength /m')
% ylabel('Magnitude /dB')
% set(gca,'Fontname','Times New Roman','fontsize',14);
% plot 70
% figure;
semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,70);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
% % plot 120
% % figure;
% semilogx((lamda), ((Acc_TriRec_120)), 'b','LineWidth',2);
% plot_line_func(Acc_TriRec_120,lamda,120);
% xlabel('Wavelength /m')
% ylabel('Magnitude /dB')
% set(gca,'Fontname','Times New Roman','fontsize',14);

function [ensco1, tridb, maxV, minV, BW, slope] = tri_rec_paral_repair(lam, para)
%三角窗+矩形窗并联幅频特性分析
N=128^2;
w=linspace(0,pi,N);

Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);

switch (lam)
    case 25
        M1 = para(1,1);
        M2 = para(1,2);
        M3 = para(1,3);
        M4 = para(1,4);
    case 42
        M1 = para(2,1);
        M2 = para(2,2);
        M3 = para(2,3);
        M4 = para(2,4);
    case 70
        M1 = para(3,1);
        M2 = para(3,2);
        M3 = para(3,3);
        M4 = para(3,4);
    case 120
        M1 = para(4,1);
        M2 = para(4,2);
        M3 = para(4,3);
        M4 = para(4,4);
end
        
%% way 2 (程序中原有的做法)

% m1 = floor(lam*2.4+1.9+0.5);
% M1 = m1 + 1 - mod(m1,2);
% m2 = floor(M1*0.3216 + 0.5);
% M2 = m2 + 1 - mod(m2,2);
% M3 = round(M1*1.1228)*2+1;
% M4 = round(M1*1.5439)*2+1;

% 矩形窗的固有形式
    M1 = M1/2;
    M2 = M2/2;
for i=1:N
%     三角窗
    Tri1(i)= ((sin(w(i)*M1/2)/sin(w(i)/2))/M1)*((sin(w(i)*M1/2)/sin(w(i)/2))/M1);
    Tri2(i)= ((sin(w(i)*M2/2)/sin(w(i)/2))/M2)*((sin(w(i)*M2/2)/sin(w(i)/2))/M2);
    
%     矩形窗
    Rec1(i)= ((sin(w(i)*M3/2)/sin(w(i)/2))/M3);
    Rec2(i)= ((sin(w(i)*M4/2)/sin(w(i)/2))/M4);
end
% 并联的系数
MulFilter = 1.036*Tri1-0.036*Tri2+0.25*(Rec1-Rec2);
% a = 1.0;
% b = 1.0;
% MulFilter = Tri1+a*0.036*(Tri1-Tri2)+b*0.25*(Rec1-Rec2);

% test
% MulFilter = Tri1-Tri2+Rec1-Rec2;
ensco1 = 1 - MulFilter;         %%这个为什么是这样的

ensco1 = MulFilter;

Acc = ensco1;
[tridb, maxV, minV, BW, slope] = filter_para(Acc);
end

function para = winLen_para_set()
list = [25,42,70,120];
for i = 1:4
    Wavelen = list(i);
    [Tri1,Tri2,Rec1,Rec2] = winlen_lamda(Wavelen);
    para(i,:) = [Tri1,Tri2,Rec1,Rec2];
end
end

function [Tri1,Tri2,Rec1,Rec2] = winlen_lamda(lamda)

Tri1 = 4.4940*lamda-0.2418; Tri1 = CalWinLen(Tri1);
Tri2 = 1.2341*lamda-0.7893; Tri2 = CalWinLen(Tri2);
Rec1 = 4.4940*lamda-0.2418; Rec1 = CalWinLen(Rec1);
Rec2 = 7.8231*lamda-1.6313; Rec2 = CalWinLen(Rec2);

out = [Tri1,Tri2,Rec1,Rec2];
end

function Len = CalWinLen(tmp)
if mod(floor(tmp),2) == 1
    Len = floor(tmp);
else
    Len = floor(tmp)-1;
end
end