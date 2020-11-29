% =========================================================================
%
%                  ���ʵ�ͨ�˲���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��16��
%   ���ߣ�
%--------------------------------------------------------------------------
%  ���ܣ� 1.���������������δ�������
%        2.����100m�������ϵģ���ΪҪ�������ʣ�����Ҫѡ�����������
%        3. ���ǵ�ͨ�˲�����ֱ���ô������Ϳ���ʵ����
%--------------------------------------------------------------------------

clear all;
close all;
%�������ֵ�ͨ�˲���D(z)
%����1
K = 38;
M = 2*K+1;
L = 58;
N = 2*L+1;
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
length = 0.25;
temp   = exp(1j*2*pi.*pesi*length);
% DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
%            .*( temp.^L - temp.^(-L-1) )...
%            ./( ( 1- exp(-1j*2*pi.*pesi*length) ) .^2) ...
%             /M /N ;
%%��һ�ֱ�ʾ�������������δ�
DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
           .*( temp.^L - temp.^(-L-1) )...
           ./( ( 1- temp.^(-1) ) .^2) ...
            /M /N ;


%% ��ͼ

figure1 = figure('Color',[1 1 1]);
semilogx( 1./pesi , 20*log10(DjpesiDz) ,'LineWidth',1);
xlabel('\lambda m');
ylabel('Mag dB');
title ('���ʵ�ͨ�˲���');
set(gca,'Fontname','Times New Roman','fontsize',14);
grid on;

figure1 = figure('Color',[1 1 1]);
semilogx( pesi , 20*log10(DjpesiDz) ,'LineWidth',1 );
xlabel('\psi��1/m��');
ylabel('Mag dB');
title ('���ʵ�ͨ�˲���');
set(gca,'Fontname','Times New Roman','fontsize',14);
grid on;

%% �µķ���

% %����2
% M = 77;
% N = 117;
% lamda = 1:0.1:10000;
% pesi  = 1./lamda;    %�ռ���Ƶ��
% length = 0.25;
% Dpesi  =    sin(M*pi.*pesi*length)      .* sin(N*pi.*pesi*length) ...
%           ./   (sin(pi.*pesi*length)) ./    (sin(pi.*pesi*length))...
%            /M/N;
% semilogx( lamda ,20*log10(Dpesi) );

% %���ʵ�ͨ�˲���
% K = 38;
% M = 2*K+1;
% L = 58;
% N = 2*L+1;
% lamda = 1:0.1:1000;
% pesi  = 1./lamda;    %�ռ���Ƶ��
% length = 0.25;
% temp   = exp(1j*2*pi.*pesi*length);
% DjpesiEz =     ( temp.^K - temp.^(-K-1) )...
%            .*( temp.^L - temp.^(-L-1) )...
%            ./ ( 1- exp(-1j*2*pi.*pesi*length) ) ...
%             /M /N ;
% figure
% semilogx( pesi , 20*log10(DjpesiEz) );
% xlabel('������m��');
% ylabel('dB');
% title ('���ʵ�ͨ�˲���');
