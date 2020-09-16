clear all
clc
close all
%step1 分析GJ-4型轨检车（三窗级联）函数幅频特性
%step2 矩形窗 三窗级联 最优化设计函数幅频特性
%step3 三窗级联 最优化设计 对比 原GJ-4型窗函数幅频特性
%step4 确认窗函数选择
lamda = 1:0.001:1000;    %空间波长
pesi  = 1 ./ lamda;      %空间频率

%% 矩形窗1
L = 38;
N = 2*L + 1;
V_z_1  = ( (exp(-1j.*pesi)).^L - (exp(-1j.*pesi)).^(-L-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / N ;

%%矩形窗2
L = 58;
N = 2*L + 1;
V_z_2  = ( (exp(-1j.*pesi)).^L - (exp(-1j.*pesi)).^(-L-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / N ;

V_z_3 = V_z_1.*V_z_2;

figure
semilogx( pesi, 20*log10(V_z_3) );

%% 三阶不等长矩形窗级联  55 69 83   
L_base = 97;

L1   = ( L_base*0.6505 - 1 )/2 ;
% L1   = 28;
R_54 = ( (exp(-1j.*pesi)).^L1 - (exp(-1j.*pesi)).^(-L1-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / (2*L1) ;

L2   = ( L_base*0.8313 - 1 )/2 ;
% L2   = 40;
R_69 = ( (exp(-1j.*pesi)).^L2 - (exp(-1j.*pesi)).^(-L2-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / (2*L2);

L3   = ( L_base- 1 ) /2 ;
% L3   = 48;
R_83 = ( (exp(-1j.*pesi)).^L3 - (exp(-1j.*pesi)).^(-L3-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / (2*L3);
R_z_3 = R_54 .* R_69 .*R_83 ;

figure
semilogx( pesi, 20*log10(V_z_3) ,'r', pesi, 20*log10(R_z_3),'g');
legend  ('二窗级联','三窗级联最优化');
