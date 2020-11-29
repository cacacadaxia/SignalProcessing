% =========================================================================
%
%                  �˲������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��27��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.��֤�˲�����׼ȷ�ԣ�����������˲����ķ���
%        2.
%        3. 
%--------------------------------------------------------------------------

clear all;
close all;
psi = 0.001:0.001:1/4;%%���fs/2
% psi = 0.001:0.001:0.01;
Omega = 2*pi*psi;
z_1 = exp(-1j*2*pi*psi/4);

%%
% z_1 = exp(-1j*2*pi*psi*0.25);
% FIR1 = (2-z_1+3*z_1.^2)./( 1-z_1+z_1.^2 );
% temp  = 5/0.25;
% FIR1 = (1-z_1.^temp)./5.*(0.25./(1-z_1) + 1/1.29);
% 
% figure1 = figure('Color',[1 1 1]);
% semilogx(1./psi,20*log10(abs(FIR1)) ,'r','LineWidth',1);
% xlabel('\lambda m');ylabel('Mag dB');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% 
% % figure1 = figure('Color',[1 1 1]);
% % semilogx(psi,20*log10(abs(FIR1)) ,'r','LineWidth',2);
% % xlabel('\psi Hz');ylabel('Mag dB');
% % set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% 
% 
% IIR1 = 1./FIR1;
% hold on;
% semilogx(1./psi,20*log10(abs(IIR1)) ,'k','LineWidth',1);
% semilogx(1./psi,20*log10(abs(IIR1.*FIR1)) ,'b--','LineWidth',1);
%% 
L = 5;
Hdelz_japan = (1 - exp(-1j*Omega*L));

Res_s = Hdelz_japan./L.*0.25*1./(1 - z_1);
H_final = Hdelz_japan/4 + Res_s;%%~~H3z - Hdelz_japan/2

H_final_fir = 1./H_final;

figure1 = figure('Color',[1 1 1]);
semilogx(1./psi,20*log10(abs(H_final)) ,'r','LineWidth',2);
hold on;
semilogx(1./psi,20*log10(abs(H_final_fir)) ,'r','LineWidth',2);
xlabel('\lambda m');ylabel('Mag dB');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;




%%
compFir = (1+L*(1-z_1))./( 4*L*(1-z_1) );
figure1 = figure('Color',[1 1 1]);
semilogx(1./psi,20*log10(abs(Hdelz_japan)) ,'r','LineWidth',2);
hold on;
semilogx(1./psi,20*log10(abs(Hdelz_japan.*compFir)) ,'b','LineWidth',2);
xlabel('\lambda m');ylabel('Mag dB');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

legend �����ַ� �����˲���


figure1 = figure('Color',[1 1 1]);
plot(psi,angle(Hdelz_japan.*compFir)/pi*180)
xlabel('\psi Hz');ylabel('�Ƕ�');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

