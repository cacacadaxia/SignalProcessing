% =========================================================================
%
%                  差分方法测试
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 11月24日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.两点差分与弦测法两种方法对比，从传递函数到幅值表达式
%        2.两者都是线性相位的，不用做太多考虑这方面的东西
%       3. 
%--------------------------------------------------------------------------

clear all;
close all;
psi = 0.001:0.001:1/6;%%最大到fs/2
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
% xlabel('\psi Hz');ylabel('角度');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%% 弦测法
% % hold on;
% figure1 = figure('Color',[1 1 1]);
% % psi = 1/200:0.001:1/10;
% lambda = 1./psi;
% L = 10;
% HlambdaMag = 1 - cos(pi./lambda*L);
% semilogx(1./psi,abs(HlambdaMag));
% xlabel('\lambda m');ylabel('Mag');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% % legend 两点差分 弦测;
% 
% %% 弦测法的理论
% hold on;
% Hz1 = (1 + exp(-1j*Omega*L))/2 - exp(-1j*Omega*L/2);
% semilogx(1./psi,abs(Hz1),'--');
% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle(Hz1)/pi*180)
% xlabel('\psi Hz');ylabel('角度');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
%% 新方案
%% ∫∫(pitch*0.25)   ?  z2
z_tmp = exp(1j*2*pi*psi/4);
Res_s = Hdelz_psi.*0.25*1./(1-z_tmp.^(-1));

% 幅值
% figure1 = figure('Color',[1 1 1]);

% hold on;
% semilogx(1./psi,abs(Res_s),'r-','LineWidth',1);
xlabel('\lambda m');ylabel('Mag');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% 相位(注意这个相位事线性的)
% 这就很奇怪，哪里来的5呢？
% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle(Res_s)/pi*180,'r-','LineWidth',1);
% xlabel('\psi Hz');ylabel('Angle (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;

%% z1+z2/2 ? z2
H3z = ( 1 + exp(-1j * Omega * L) )/2;
% hold on;
% semilogx(1./psi,abs(H3z),'b-','LineWidth',1);


%% 找到拟合关系
semilogx(1./psi,20*log10(abs(Res_s - 0.1/1.1* Hdelz_japan + Hdelz_japan/2)),'r-','LineWidth',1);
hold on
semilogx(1./psi,20*log10(abs(H3z + Hdelz_japan/2)),'g-','LineWidth',1);
H8z = (Res_s + Hdelz_japan/1.12)./(H3z + Hdelz_japan);
% semilogx(1./psi,abs(H8z),'r-','LineWidth',1);
grid on;
semilogx(1./psi,20*log10(abs(Res_s + Hdelz_japan/4)),'k-','LineWidth',1);

title('集合相互重合的两个曲线')
legend H_3(s)+H_1(s)/1.12 H_0(s)+H_1(s);
xlabel('\lambda m');ylabel('Mag');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;


%%
figure1 = figure('Color',[1 1 1]);
plot(psi,angle( H8z )/pi*180,'k-','LineWidth',1);
hold on;
plot(psi,angle( Res_s + Hdelz_japan/1.12 )/pi*180,'b-','LineWidth',1);
plot(psi,angle( H3z + Hdelz_japan )/pi*180,'r-','LineWidth',1);
xlabel('\psi Hz');ylabel('Angle (deg)');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
figure1 = figure('Color',[1 1 1]);
% H_tmp = (Res_s + Hdelz_japan/1.12) ./ (H3z + Hdelz_japan);
semilogx(1./psi,abs(Res_s + Hdelz_japan/4),'k-','LineWidth',1);
hold on;
% semilogx(1./psi,abs(H_tmp),'r-','LineWidth',1);
xlabel('\lambda m');ylabel('Mag');
set(gca,'Fontname','Times New Roman','fontsize',14);grid on;



%%  最终的方案
H_res = Res_s + Hdelz_japan/1.29;
% H_res = Res_s + Hdelz_japan/1.40;
H_res_final = ( H_res*2 + exp(-1j*Omega*L) )./3;

semilogx(1./psi,abs(H_res_final),'g-','LineWidth',2);
% 
semilogx(1./psi,abs(Hdelz_japan),'y-','LineWidth',2);
legend H_4(s) H_6(s) H_1(s)
title('三种方案对比')

% figure1 = figure('Color',[1 1 1]);
% plot(psi,angle( H_res_final )/pi*180,'b-','LineWidth',1);
% hold on;
% plot(psi,angle(Res_s+Hdelz_japan/4 )/pi*180,'r-','LineWidth',1);
% xlabel('\psi Hz');ylabel('Angle (deg)');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% hold on;
% legend H_6(s) H_4(s)
% title('相频曲线')
% % plot(psi,angle(Res_s )/pi*180,'k-','LineWidth',1);
% % plot(psi,angle(Hdelz_japan/4  )/pi*180,'g-','LineWidth',1);

%% 设计逆滤波器


% figure1 = figure('Color',[1 1 1]);
% % semilogx(1./psi,abs(Pz.*H_res_final),'g-','LineWidth',2);
% % hold on;
% % semilogx(1./psi,abs(Pz),'r-','LineWidth',2);
% 
% H5z = H3z + Hdelz_japan;
% H7z = (H5z + exp(-1j*Omega*L)./2)/3*2;
% semilogx(1./psi,abs(H8z),'k-','LineWidth',2);
% 
% xlabel('\lambda m');ylabel('Mag');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;






%% z1+z2/2 == ? ∫∫(pitch*0.25)
% Res_s2 = Res_s.*( 2./( 1 + exp(-1j * Omega * L) ) );
% % 幅值
% figure1 = figure('Color',[1 1 1]);
% semilogx(1./psi,abs(Res_s2),'r-','LineWidth',1);
% xlabel('\lambda m');ylabel('Mag');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% 
% 
% %% (z1+z2)/2 ? z2
% Hs4 = ( 1 + exp(-1j * Omega * L))./2;
% figure1 = figure('Color',[1 1 1]);
% semilogx(1./psi,abs(Hs4),'r-','LineWidth',1);
% xlabel('\lambda m');ylabel('Mag');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
%% 测试两者是否相等(论文中的2米弦测法)
%%这种方法太局限，并不能实用
%%只有是2的时候才相等
% psi = 0.001:0.001:2;
% lambda = 1./psi;
% lambda = 200;
% L = 2;
% x = 1:0.01:200;
% tmp0 = sin(2*pi/lambda*x) - sin(2*pi/lambda*(x-L));
% tmp1 = L/2*2*sin(2*pi/lambda)*cos(2*pi/lambda*x - 2*pi/lambda*L/2);
% 
% figure1 = figure('Color',[1 1 1]);
% plot(x,tmp0,'r','LineWidth',1);hold on;plot(x,tmp1,'k--','LineWidth',1)
% % xlabel('');ylabel('');
% set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% legend 1 2
%% 两点差分法设计


