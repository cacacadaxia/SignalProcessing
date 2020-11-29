  %
%                  ������˲����ķ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 11��9��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���׵�ͨ�˲�����ȥ�Ʊ��˲�����Ϊ�о��Ķ���
%        2.�Ѿ�����˴���ı�д
%       3. 
%--------------------------------------------------------------------------

%% �����˲�����ȥ�Ʊ��˲�����������
clear all;
close all;
f = 0.01:0.01:250;%%����500�ͺܹ֡�
fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz
% ��������
Omega_1 = (10^5)/(2^17);%%��ֹƵ��
Omega_2 = 10^5/2^14;
omiga =  2*pi.*f;           % ģ���Ƶ�� rad/sec
% Fs = Omega_1 ./ (1j.*(omiga) + Omega_1);%%����һ�׵�
Fs = Omega_2^2./( (1j.*omiga).^2 + (1j.*omiga).*Omega_2 + Omega_2^2 );

%% 50km/h
figure1 = figure('Color',[1 1 1]);
v1 = 50/3.6;
semilogx(v1./f , 20*log10(abs(Fs)),'LineWidth',1);hold on;
Ti = 0.25/v1;
z_1 = exp( - 1j*2*pi* f/v1 *0.25);
Cz1 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
a = 1/2*Omega_2; b = sqrt(1-0.5^2)*Omega_2;
a2b2 = (10^5/2^14)^2;
w2 = (10^5)/(2^14);
a = w2/2; aa = a^2;
bb = 3*(w2^2)/4;
biir = [1+a*Ti+aa*Ti*Ti ,-2*(1-aa*Ti*Ti), 1-a*Ti+aa*Ti*Ti];
% Gz1 = ( (1-z_1).^2 + (1 - z_1.^2)*a*Ti + (a*Ti + 2*a*Ti*z_1 + a*Ti*z_1.^2) ) ./ (Ti^2*a2b2);
Gz1 = (biir(1) + biir(2).*z_1 + biir(3).*z_1.^2) ./ ((aa+bb)*Ti*Ti);
hold on;semilogx(v1./f , 20*log10(abs(Gz1)));

%% 120km/h
v2 = 120/3.6;
semilogx(v2./f , 20*log10(abs(Fs)));
Ti = 0.25/v2;
z_1 = exp( - 1j*2*pi* f/v2 *0.25);
Cz2 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
% Gz2 = ( (1-z_1).^2 + (1 - z_1.^2)*a*Ti + (a*Ti + 2*a*Ti*z_1 + a*Ti*z_1.^2) ) ./ (Ti^2*a2b2);
biir = [1+a*Ti+aa*Ti*Ti ,-2*(1-aa*Ti*Ti), 1-a*Ti+aa*Ti*Ti];
Gz2 = (biir(1) + biir(2).*z_1 + biir(3).*z_1.^2) ./ ((aa+bb)*Ti*Ti);
hold on;semilogx(v2./f , 20*log10(abs(Gz2)));

%% 350km/h
v3 = 200/3.6;
v3 = 62.5;
semilogx(v3./f , 20*log10(abs(Fs)));
Ti = 0.25/v3;
z_1 = exp( - 1j*2*pi* f/v3 * 0.25);
% Cz3 = ( exp(Omega_1*Ti/2) - exp(- Omega_1*Ti/2)*z_1 )/(Omega_1*Ti);
% Gz3 = ( (1-z_1).^2 + (1 - z_1.^2)*a*Ti + (a*Ti + 2*a*Ti*z_1 + a*Ti*z_1.^2) ) ./ (Ti^2*a2b2);
biir = [1+a*Ti+aa*Ti*Ti ,-2*(1-aa*Ti*Ti), 1-a*Ti+aa*Ti*Ti];
Gz3 = (biir(1) + biir(2).*z_1 + biir(3).*z_1.^2) ./ ((aa+bb)*Ti*Ti);
hold on;semilogx(v3./f , 20*log10(abs(Gz3)));
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m');
ylabel('Mag /dB');
legend 50km/h 50km/h 120km/h 120km/h 350km/h 350km/h;
% legend 50km/h 120km/h 350km/h
grid on;

%% �����Ļ�һ��ȥ�Ʊ��˲���
figure1 = figure('Color',[1 1 1]);
semilogx( v1./f , 20*log10(abs(Gz1)));hold on;
semilogx( v2./f , 20*log10(abs(Gz2)));
semilogx( v3./f , 20*log10(abs(Gz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('������˲���+ȥ�Ʊ��˲���')
legend 50km/h 120km/h 200km/h

%% �����+ȥ�Ʊ�
figure1 = figure('Color',[1 1 1]);
semilogx(v1./f , 20*log10(abs(Fs.*Gz1)));
hold on;
semilogx(v2./f , 20*log10(abs(Fs.*Gz2)));
semilogx(v3./f , 20*log10(abs(Fs.*Gz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('������˲���+ȥ�Ʊ��˲���')
legend 50km/h 120km/h 200km/h

%% �µ��о��� ������
%%��Ҫע����ǣ������ֻ��2Hz��λ�þ�Ӧ�ò�����
%%��Ϊ�ռ���Ĳ���Ƶ�ʾ���4Hz�����賵���ٶ�Ϊ100m/s����ô������Ϊ100*4=400Hz
%%psi_max = f./v = 400Hz./v = 4Hz�����Ի���4Hz
figure1 = figure('Color',[1 1 1]);
semilogx(f./v1 , 20*log10(abs(Fs.*Gz1)));
hold on;
semilogx(f./v2 , 20*log10(abs(Fs.*Gz2)));
semilogx(f./v3 , 20*log10(abs(Fs.*Gz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\psi /Hz')
ylabel('Mag /dB')
title('������˲���+ȥ�Ʊ��˲���')
legend 50km/h 120km/h 200km/h
%%ʵ������ѡ���൱�ڵ�ͨ�˲����������źŶ�ʧ�㲻�����ˣ���Ϊ��Ƶ���ź��ڵ�Ƶ����ʱ
%%(����Ĳ�����ָ�ռ����ѡ��)�������Ϊ���Ǵ����źŵĻ���������ͨ�˲����൱�ڿ���
%%��ɿ�����Ĺ��ܡ�

%%��ȻҲ���ܲ�û�п�����Ĺ��ܡ�

%%
%%���������Ĺ��ܻ��Ǵ��ڵİɣ�

%% ��֤һ�����������Ĺ����ǲ��ǳ�����
% figure;semilogx(f,20*log10(abs(Fs)));
%%Ϊʲô��ֹƵ����1Hz���ң�����÷��Ƚ�����

%%
% figure1 = figure('Color',[1 1 1]);
% semilogx(f, 20*log10(abs(Fs.*Gz1)));
% hold on;
% semilogx(f , 20*log10(abs(Fs.*Gz2)));
% semilogx(f , 20*log10(abs(Fs.*Gz3)));
% grid on;
% set(gca,'Fontname','Times New Roman','fontsize',14);
% xlabel('\psi /Hz')
% ylabel('Mag /dB')
% title('������˲���+ȥ�Ʊ��˲���')
% legend 50km/h 120km/h 200km/h


