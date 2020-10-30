% =========================================================================
%
%                  ����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��22��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------
clear all
close all

Omega2 = 10^5/2^14;
b = [Omega2^2];
a = [1,Omega2,Omega2^2];
[H,Om] = freqs(b,a);

figure1 = figure('Color',[1 1 1]);
semilogx(Om/2/pi , 20*log10(abs(H)),'LineWidth',2);
hold on
%%
%%�޸�T��ֵ������ı��˲���������
T = 1/500;
b = [Omega2^2*T^2];
a = [(1+Omega2*T+(Omega2*T)^2) ,  -(2+Omega2*T) , 1];
[h,f] = freqz(b,a,800000,1/T);
% figure;
semilogx( f , 20*log10(abs(h)),'--r','LineWidth',1);
xlabel('f /Hz')
ylabel('Mag /dB');
set(gca,'Fontname','Times New Roman','fontsize',16);

%%Ϊʲô��ֹƵ����1Hz��

%% �ռ���
%%�������ⲻһ���ԣ�
v = 100;v = v/3.6;lamda = v./f;
figure;
semilogx( lamda , 20*log10(abs(h)),'--r','LineWidth',1);
hold on;
v = 300;v = v/3.6;lamda = v./f;
semilogx( lamda , 20*log10(abs(h)),'--r','LineWidth',1);
xlabel('���� /m')
ylabel('Mag /dB');
set(gca,'Fontname','Times New Roman','fontsize',16);




