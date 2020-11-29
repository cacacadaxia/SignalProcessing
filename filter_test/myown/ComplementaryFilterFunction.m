% =========================================================================
%
%                  �����˲����Ĺ��ܲ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��18��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ֲ�û�����Ϊ1�Ľ�����֣��������ڽ���������⣬Ҳ�������ĵ�д����
%        2.��Ϊ��Ҫ����ǰ����˲���������Fs����Bz��
%        3. 
%--------------------------------------------------------------------------


clear all;
close all;
v = 100;
T = 0.25/v;
fs = 1./T;
f = 0.0001:0.0001:4*v/2;

%% H3s ��ͨ�˲���
Omega1 = 0.76;
Omega2 = 6.10;
Omega3 = 0.38;

tp1 = Omega1*Omega2^2*Omega3^3;
tp2 = Omega2*Omega3*(Omega1*Omega2 + Omega2*Omega3 + Omega1 * Omega3);
tp3 = Omega2^2*Omega1;
b = [tp2,tp1];
a = tp3*[1,Omega3,Omega3^2];

s = 2*pi*f*1j;
temp1 = s.^2 + Omega3.*s+Omega3^2;
H3s = (tp1 + tp2.*s) ./(tp3*( temp1 ));

figure1 = figure('Color',[1 1 1]);
semilogx(f , 20*log10(abs(H3s)),'LineWidth',1);



%%
z_1 = exp(-1j*2*pi.*f/fs);
H3z = (Omega3*T)^2 + (Omega3*T)*(1+Omega3/Omega1+Omega3/Omega2)*(1-z_1);
H3z = H3z./( 1 + Omega3*T +(Omega3*T)^2  - ( (2+Omega3*T)*z_1 - z_1.^2 ) );
hold on;
semilogx(f , 20*log10(abs(H3z)),'LineWidth',1);
% legend ģ���˲��� �����˲���

%%
% figure;plot(f , 180*(angle(H3z))/pi,'LineWidth',1);

%% ��ͨ�˲���
% Omega1 = 10^5/2^17;
% Omega2 = 10^5/2^14;
% Omega3 = 10^5/2^18;

zeta = Omega1+Omega2+Omega3;
nang = Omega2^2 + Omega2*Omega3 + Omega2*Omega1 + Omega1*Omega3 +Omega3^2;
kai = (Omega1 + Omega3)*Omega2^2 + Omega2*Omega1*Omega3 + (Omega1 + Omega2)*Omega3^2;

temp2 = s.^2 + Omega2.*s+Omega2^2;
P2s = s.^5 + zeta*s.^4 + nang*s.^3 + kai*s.^2;
P2s = P2s./( Omega1.*temp1.*temp2 )./s;

%% plot
hold on;
semilogx(f , 20*log10(abs(P2s)),'LineWidth',1);
semilogx(f , 20*log10(abs(-P2s + H3z*2)),'LineWidth',1);
legend ��ͨ�˲��� ��ͨ�˲��� �����˲���
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('f /Hz');
ylabel('Mag /dB');

%%
figure;plot(f , 180*(angle(P2s))/pi,'LineWidth',1);

