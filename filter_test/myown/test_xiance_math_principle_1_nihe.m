% =========================================================================
%
%                  ��ַ�������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��24��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���������Ҳⷨ���ַ����Աȣ��Ӵ��ݺ�������ֵ���ʽ
%        2.���߶���������λ�ģ�������̫�࿼���ⷽ��Ķ���
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;
psi = 0.001:0.001:1/4;%%���fs/2
% psi = 0.001:0.001:0.01;
Omega = 2*pi*psi;
L = 5;
Hdelz_japan = (1 - exp(-1j*Omega*L));
Hdelz_psi = Hdelz_japan./L;
figure1 = figure('Color',[1 1 1]);
% semilogx(1./psi,abs(Hdelz_japan),'y','LineWidth',2);
xlabel('\lambda m');ylabel('Mag');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle(Hdelz_japan)/pi*180)
% xlabel('\psi Hz');ylabel('�Ƕ�');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%% �ҡ�(pitch*0.25)   ?  z2
z_tmp = exp(1j*2*pi*psi/4);
Res_s = Hdelz_psi.*0.25*1./(1-z_tmp.^(-1));

% ��ֵ
% figure1 = figure('Color',[1 1 1]);

% hold on;
% semilogx(1./psi,abs(Res_s),'r-','LineWidth',1);
xlabel('\lambda m');ylabel('Mag');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% ��λ(ע�������λ�����Ե�)
% ��ͺ���֣���������5�أ�
% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle(Res_s)/pi*180,'r-','LineWidth',1);
% xlabel('\psi Hz');ylabel('Angle (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%% z1+z2/2 ? z2
H3z = ( 1 + exp(-1j * Omega * L) )/2;
% hold on;
% semilogx(1./psi,abs(H3z),'b-','LineWidth',1);


%% �ҵ���Ϲ�ϵ
% semilogx(1./psi,20*log10(abs(Res_s + Hdelz_japan./1.1)),'r-','LineWidth',1);
% hold on
% semilogx(1./psi,20*log10(abs(H3z + Hdelz_japan)),'g-','LineWidth',1);
H8z = (Res_s + Hdelz_japan/1.1)./(H3z + Hdelz_japan);
% semilogx(1./psi,20*log10(abs(H8z)),'m-','LineWidth',1);
grid on;
% semilogx(1./psi,20*log10(abs(Res_s + Hdelz_japan/4)),'k-','LineWidth',1);


%%
H_final = Hdelz_japan*0.4090 + Res_s;%%~~H3z - Hdelz_japan/2
semilogx(1./psi,20*log10(abs(H_final)),'m-','LineWidth',1);
grid on;
% title('�����໥�غϵ���������')
legend H_3(s)+H_1(s)/1.12 H_0(s)+H_1(s);
xlabel('\lambda m');ylabel('Mag');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%%
figure1 = figure('Color',[1 1 1]);
plot(psi,angle(H_final)/pi*180,'r-','LineWidth',1);
hold on;
plot([0,0.1,0.2],[-36.2+73.799999999999997,-36.2,-110],'b-','LineWidth',1);
% plot(psi,angle(Res_s)/pi*180,'k-','LineWidth',1);
% plot(psi,angle(Hdelz_japan)/pi*180,'b-','LineWidth',1);
xlabel('\psi Hz');ylabel('Angle (deg)');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%%
% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle( H8z )/pi*180,'k-','LineWidth',1);
% hold on;
% plot(psi,angle( Res_s + Hdelz_japan/1.12 )/pi*180,'b-','LineWidth',1);
% plot(psi,angle( H3z + Hdelz_japan )/pi*180,'r-','LineWidth',1);
% xlabel('\psi Hz');ylabel('Angle (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% figure1 = figure('Color',[1 1 1]);
% % H_tmp = (Res_s + Hdelz_japan/1.12) ./ (H3z + Hdelz_japan);
% semilogx(1./psi,abs(Res_s + Hdelz_japan/4),'k-','LineWidth',1);
% hold on;
% % semilogx(1./psi,abs(H_tmp),'r-','LineWidth',1);
% xlabel('\lambda m');ylabel('Mag');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;



%%  ���յķ���
% H_res = Res_s + Hdelz_japan/1.29;
% % H_res = Res_s + Hdelz_japan/1.40;
% H_res_final = ( H_res*2 + exp(-1j*Omega*L) )./3;
% 
% semilogx(1./psi,abs(H_res_final),'g-','LineWidth',2);
% % 
% semilogx(1./psi,abs(Hdelz_japan),'y-','LineWidth',2);
% legend H_4(s) H_6(s) H_1(s)
% title('���ַ����Ա�')

% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle( H_res_final )/pi*180,'b-','LineWidth',1);
% hold on;
% plot(psi,angle(Res_s+Hdelz_japan/4 )/pi*180,'r-','LineWidth',1);
% xlabel('\psi Hz');ylabel('Angle (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% hold on;
% legend H_6(s) H_4(s)
% title('��Ƶ����')
% % plot(psi,angle(Res_s )/pi*180,'k-','LineWidth',1);
% % plot(psi,angle(Hdelz_japan/4  )/pi*180,'g-','LineWidth',1);


