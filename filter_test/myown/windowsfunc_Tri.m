
% =========================================================================
%
%                  ���δ�
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��23��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�������Ļ�ͼ���ο��ĵ���������·���������ƽ˳��⼼���о�.docx��
%        2.�������Ҫ�������˸��ַ���֮��Ĺ�ϵ
%        3. ���ֻ�ͼ�ķ������ǳ���Ҫ
%--------------------------------------------------------------------------

close all;
clear all;

%% ���δ�
%% 
M = 30;
N = 2*M+1;
psi = linspace(0,pi,10000);
Omega = 2.*pi./psi.*0.25;
Wn = sin(N.*psi/2)./sin(psi./2)./N;
figure1 = figure('Color',[1 1 1]);semilogx(Omega,20*log10(abs(Wn)));
grid on;


%%
M = 30;
N = 2*M+1;
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
Omega = 2*pi*pesi*0.25;
Wn = sin(N.*Omega/2)./sin(Omega./2)./N;
hold on;semilogx(lamda,20*log10(abs(Wn)));

%% ���Ǵ�
% M = 30;
% N = 2*M+1;
% psi = linspace(0,pi,10000);
% Omega = 2.*pi./psi.*0.25;
% Wn = (sin(M*psi/2)./sin(psi/2)).^2/M^2;
% figure;semilogx(Omega,20*log10(abs(Wn)));
% grid on;

%% freqz
b = ones(1,N)/N;
[h,f] = freqz(b,[zeros(1,N-1),1],10000,4);
hold on;semilogx(1./f,20*log10(h));grid on;

%% z
K = 30;
M = 2*K+1;
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
temp = exp(1j*2*pi.*pesi*0.25);
%%��һ�ֱ�ʾ�������������δ�
DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
           ./( ( 1- temp.^(-1) ) ) ...
            /M ;
hold on;semilogx( lamda , 20*log10(DjpesiDz) );

%%

pahse = angle(DjpesiDz)/pi*180;
figure;semilogx( lamda , pahse );
%% ����˲����������ٶ��ˣ�
%%��ʵ�ֵ�������
% lamda = 1:0.001:900;
% psi = 1./lamda;
% fs = 0.25;
% z = exp(1j*2*pi*psi/fs);
% Hz1 = 1/83*(z.^41 - z.^(-42))./( 1-z.^(-1) );
% Hz2 = 1/249*( z.^124 - z.^(-125) )./(1 - z.^(-1));
% Hz = 1 - Hz1 - 1/8*(Hz1 - Hz2);
% figure;semilogx(lamda , Hz);






