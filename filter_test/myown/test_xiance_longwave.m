

% =========================================================================
%
%                  �Ҳⷨ������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��23��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;
%% 1m�Ҳ�
lambda = 0.5:0.01:500;%%��С����Ϊ0.5m
L = 1;
HlambdaMag = 1 - cos(pi./lambda*L);
Hfunc = @(lambda)(1 *(1 - cos(pi./lambda*L)));
figure1 = figure('Color',[1 1 1]);loglog((1./lambda),(HlambdaMag));xlabel('Ƶ�� ');ylabel('M(x)/f(x)')
figure1 = figure('Color',[1 1 1]);semilogx((lambda),(HlambdaMag));xlabel('���� /m ');ylabel('M(x)/f(x)')

%% ���˲��� I(z)
% lambdaI = 1:1:100;      %%FIR������500��ô��
% psiI = 1./lambdaI;

Ifunc = @(lambda)(1 ./ (1 - cos(pi./lambda*L)));
% MagI = Ifunc(lambdaI);

% hold on;semilogx(lambdaI ,MagI );
%%����ֱ�ӽ���ifft�任��Ӧ����Ҫ��Ƶ���Ͻ��в�ֵ
% figure;semilogx(1./lambdaI ,MagI ,'-xr');
% timedomain = ifft(MagI);
%%������500�����˲�������ô��Ȼû�취����

%%
% 1:10
psiI = 0.1:0.01:1;
N = length(psiI);
lambdaI = 1./psiI;


%% �渵��Ҷ�任
hn = 0;
for i=1:N
    hn(i) = 
end



