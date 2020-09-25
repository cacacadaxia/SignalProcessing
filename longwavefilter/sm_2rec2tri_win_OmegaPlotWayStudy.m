
% =========================================================================
%
%                   长波滤波器原理探究
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.旧版本的长波滤波器的（70m波长）副频曲线对比
%        2.
%        3. 
%--------------------------------------------------------------------------


clear all;
close all;
% 横坐标配置
N=128^2;
w=linspace(0,pi,N);
lamda = 2*pi./w./4;%%这样表示波长很奇怪

%%  参数配置
%三角窗+矩形窗并联幅频特性分析
N=128^2;
Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);  
%% way 2 (程序中原有的做法)
M1 = 160;
M2 = 40;
M3 = 160;
M4 = 280;
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
ensco1 = 1 - MulFilter;         %%这个为什么是这样的
Acc_TriRec_70 = ensco1;

%%
fig = figure;

semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,70);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);


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