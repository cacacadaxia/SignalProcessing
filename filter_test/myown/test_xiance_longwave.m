

% =========================================================================
%
%                  弦测法的讨论
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月23日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.
%        2.
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;
%% 1m弦测
lambda = 0.5:0.01:500;%%最小波长为0.5m
L = 1;
HlambdaMag = 1 - cos(pi./lambda*L);
Hfunc = @(lambda)(1 *(1 - cos(pi./lambda*L)));
figure1 = figure('Color',[1 1 1]);loglog((1./lambda),(HlambdaMag));xlabel('频率 ');ylabel('M(x)/f(x)')
figure1 = figure('Color',[1 1 1]);semilogx((lambda),(HlambdaMag));xlabel('波长 /m ');ylabel('M(x)/f(x)')

%% 逆滤波器 I(z)
% lambdaI = 1:1:100;      %%FIR初定是500这么长
% psiI = 1./lambdaI;

Ifunc = @(lambda)(1 ./ (1 - cos(pi./lambda*L)));
% MagI = Ifunc(lambdaI);

% hold on;semilogx(lambdaI ,MagI );
%%不能直接进行ifft变换，应该需要在频率上进行插值
% figure;semilogx(1./lambdaI ,MagI ,'-xr');
% timedomain = ifft(MagI);
%%如果变成500级的滤波器，那么显然没办法运算

%%
% 1:10
psiI = 0.1:0.01:1;
N = length(psiI);
lambdaI = 1./psiI;


%% 逆傅立叶变换
hn = 0;
for i=1:N
    hn(i) = 
end



