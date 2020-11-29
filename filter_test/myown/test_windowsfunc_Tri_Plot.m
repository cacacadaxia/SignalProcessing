
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
% 
%   1102��
%   1.�����˲���+����
%    11.28:
%   1.���������ĳ����˲�����δ���
%--------------------------------------------------------------------------

close all;
clear all;

%%
figure1 = figure('Color',[1 1 1]);
M = 30;
N = 2*M+1;
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
Omega = 2*pi*pesi*0.25;
Wn = sin(N.*Omega/2)./sin(Omega./2)./N;
semilogx(lamda,20*log10(abs(Wn)));

%% freqz
b = ones(1,N)/N;
[h,f] = freqz(b,[zeros(1,N-1),1],10000,4);
hold on;semilogx(1./f,20*log10(h));
grid on;

%% z
K = 30;
M = 2*K+1;
% lamda = 1:0.1:1000;
% pesi  = 1./lamda;           %�ռ���Ƶ��
pesi = 0.001:0.0001:1;
lamda = 1./pesi;
temp = exp(1j*2*pi.*pesi*0.25);
%%��һ�ֱ�ʾ�������������δ�
DjpesiDz =     ( temp.^K - temp.^(-K-1) )...
           ./( ( 1- temp.^(-1) ) ) ...
            /M ;
% hold on;
semilogx( lamda , 20*log10(DjpesiDz) ,'LineWidth',1);
xlabel('\lambda /m');ylabel('Mag /dB')
set(gca,'Fontname','Times New Roman','fontsize',16);

%% �����˲�����������(�е�����)
k = 70;
m = floor(0.829*k);
n = floor(0.646*k);
z_1 = temp.^(-1);
Tz = 1-1/k/m/n.*(1-z_1.^k).*(1-z_1.^m).*(1-z_1.^n)./(1-z_1).^3;
figure1 = figure('Color',[1 1 1]);
semilogx( lamda , 20*log10(Tz) ,'LineWidth',1);
xlabel('\lambda /m');ylabel('Mag /dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
grid on;




%% ���ײ�����
% H_integral = 1./( 1 - 2*temp.^(-1) + temp.^(-2) );
% Hz = DjpesiDz.*H_integral;
% hold on;grid on;
% % semilogx(lamda , 20*log10(Hz),'LineWidth',1);
% legend ;
%% fdatool
% b = load('filter1.mat');
% [h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],length(lamda)*2,4);
% figure1 = figure('Color',[1 1 1]);semilogx(1./f(1:length(f)/2),20*log10(h(1:length(f)/2)));grid on;xlabel('���� /m');
% % figure1 = figure('Color',[1 1 1]);semilogx(1./f(length(f)/2+1:end),20*log10(h(length(f)/2+1:end)));grid on;xlabel('���� /m');
% ylabel('��ֵ /dB')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% Hz2 = h(1:length(f)/2).*H_integral.';
% hold on;
% semilogx(lamda , 20*log10(abs(Hz2)));
% %%����õ����ۣ��źŶ��ײ�ֵĽ�ֹ�������źű����ֹ������һ���¡��Զ��ײ�ֽ����˲����൱�ڶ��źŽ����˲���
% %%�����ڷ����źž��������˲�����Ƶ��ʱ�������������Ľ����




