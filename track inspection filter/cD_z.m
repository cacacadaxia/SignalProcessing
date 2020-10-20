% =========================================================================
%
%                  曲率低通滤波器
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月16日
%   作者：
%--------------------------------------------------------------------------
%  功能： 1.看起来是两个矩形窗及联的
%        2.留下100m波长以上的，因为要计算曲率，就是要选择更长波长的
%        3. 
%--------------------------------------------------------------------------

clear all;
close all;
%曲率数字低通滤波器D(z)
%方法1
K = 38;
M = 2*K+1;
L = 58;
N = 2*L+1;
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %空间域频率
length = 0.25;
temp   = exp(1j*2*pi.*pesi*length);
% DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
%            .*( temp.^L - temp.^(-L-1) )...
%            ./( ( 1- exp(-1j*2*pi.*pesi*length) ) .^2) ...
%             /M /N ;
%%换一种表示，就是两个矩形窗
DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
           .*( temp.^L - temp.^(-L-1) )...
           ./( ( 1- temp.^(-1) ) .^2) ...
            /M /N ;

%% 画图

figure
semilogx( 1./pesi , 20*log10(DjpesiDz) );
xlabel('波长');
ylabel('dB');
title ('曲率低通滤波器');

figure;
semilogx( pesi , 20*log10(DjpesiDz) );
xlabel('Ψ（1/m）');
ylabel('dB');
title ('曲率低通滤波器');


%%

% %方法2
% M = 77;
% N = 117;
% lamda = 1:0.1:10000;
% pesi  = 1./lamda;    %空间域频率
% length = 0.25;
% Dpesi  =    sin(M*pi.*pesi*length)      .* sin(N*pi.*pesi*length) ...
%           ./   (sin(pi.*pesi*length)) ./    (sin(pi.*pesi*length))...
%            /M/N;
% semilogx( lamda ,20*log10(Dpesi) );

% %曲率低通滤波器
% K = 38;
% M = 2*K+1;
% L = 58;
% N = 2*L+1;
% lamda = 1:0.1:1000;
% pesi  = 1./lamda;    %空间域频率
% length = 0.25;
% temp   = exp(1j*2*pi.*pesi*length);
% DjpesiEz =     ( temp.^K - temp.^(-K-1) )...
%            .*( temp.^L - temp.^(-L-1) )...
%            ./ ( 1- exp(-1j*2*pi.*pesi*length) ) ...
%             /M /N ;
% figure
% semilogx( pesi , 20*log10(DjpesiEz) );
% xlabel('波长（m）');
% ylabel('dB');
% title ('曲率低通滤波器');
