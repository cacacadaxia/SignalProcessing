
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
%        2.不同滤波器对比，首先运行new_report_2tri_2rec_final.m
%        3. 注意三角窗用M，矩形窗用N，这是正确方法
%        4.
%        5.
% 1120
%       1.加上两个滤波器及联和并联的对比
%--------------------------------------------------------------------------


%%
clear all;
close all;
lamda = 1:0.1:10000;
pesi  = 1./lamda;           %空间域频率
omega = 2*pi*pesi*0.25;


%%  参数配置
%三角窗+矩形窗并联幅频特性分析
N = length(lamda);
Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);  

%% 调整窗长
N1 = 539;
N2 = 147;
N3 = 539;
N4 = 937;
% 			1571	431	1571	2735	350m
%% 计算半窗长
M1 = (N1-1)/2;%%三角窗的半窗长
M2 = (N2-1)/2;%%三角窗的半窗长

%%
for i=1:N
    % 三角窗
    Tri1(i)= ((sin(omega(i)*M1/2)/sin(omega(i)/2))/M1)*((sin(omega(i)*M1/2)/sin(omega(i)/2))/M1);
    Tri2(i)= ((sin(omega(i)*M2/2)/sin(omega(i)/2))/M2)*((sin(omega(i)*M2/2)/sin(omega(i)/2))/M2);
    % 矩形窗
    Rec1(i)= ((sin(omega(i)*N3/2)/sin(omega(i)/2))/N3);
    Rec2(i)= ((sin(omega(i)*N4/2)/sin(omega(i)/2))/N4);
end
% 并联的系数
MulFilter = 1.036 * Tri1 - 0.036*Tri2 + 0.25 * (Rec1-Rec2);

% 级联的滤波器
MulFilter2 = Tri1.*Tri2.*Rec1.*Rec2;

ensco1 = 1 - MulFilter;
ensco2 = 1 - MulFilter2;
Acc_TriRec_70 = ensco1;


%%
figure1 = figure('Color',[1 1 1]);
% semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,120);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
hold on;
plot_line_func(abs(ensco2),lamda,120);
%% 对比fdatool
% b = load('filter_120m.mat');
% [h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],10000,4);
% % figure1 = figure('Color',[1 1 1]);
% semilogx(1./f, abs(h) ,'LineWidth','1');grid on;xlabel('波长 /m');
% ylabel('幅值 /dB')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% legend 窗函数 FIR高阶滤波器
