

% =========================================================================
%
%                  ������˲����ķ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.һ�׵�ͨ�˲�����ȥ�Ʊ��˲�����Ϊ�о��Ķ���
%        2.
%       3. 
%--------------------------------------------------------------------------

%% ģ���˲����������˲���


clear;
close all;
%һ��ģ�⿹���˲���
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;   % ģ��Ƶ��Hz
omiga =  2*pi.*f;           % ģ���Ƶ�� rad/sec
Bs = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
semilogx( f , 20*log10(abs(Bs)) ,'LineWidth',1);
xlabel('ģ��Ƶ�� (Hz)');
ylabel('dB');            %��ֹƵ��Ϊ 0.76/2/pi = 0.121Hz
title('һ�׿����˲���');
grid on
f = 0.01:0.01:250;%%����500�ͺܹ֡�
fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz
z_1 = exp( - 1j*2*pi*f/fs);%%f����

% ��������
Omega_1 = (10^5)/(2^17);
v = 50/3.6;
T = 1/fs;
b = [ Omega_1 * T 0];
a = [ 1 + Omega_1*T ,  (-1)];
% ��ͼ
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx(f , 20*log10(abs(Bz)),'LineWidth',1);
legend ģ���˲��� �����˲���
set(gca,'Fontname','Times New Roman','fontsize',14);

 
%% �����ٶȶ��٣����չ̶�ʱ������������˲�����ģ���˲����ڿռ����ϵı���һ��
clear;
% close all;
%һ��ģ�⿹���˲���
omiga1 = 10^5 / 2^17;
f  = 0.001 : 0.0001 :100;   % ģ��Ƶ��Hz
omiga =  2*pi.*f;           % ģ���Ƶ�� rad/sec
Bs = omiga1 ./ (1j.*(omiga) + omiga1);
figure1 = figure('Color',[1 1 1]);
v = 100/3.6;
semilogx( v./f , 20*log10(abs(Bs)) ,'LineWidth',1);
xlabel('���� /m');
ylabel('dB');            %��ֹƵ��Ϊ 0.76/2/pi = 0.121Hz
title('һ�׿����˲���');
grid on
f = 0.01:0.01:250;%%����500�ͺܹ֡�
fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz
z_1 = exp( - 1j*2*pi*f/fs);%%f����

% ��������
Omega_1 = (10^5)/(2^17);
T = 1/fs;
b = [ Omega_1 * T 0];
a = [ 1 + Omega_1*T ,  (-1)];
% ��ͼ
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx(v./f , 20*log10(abs(Bz)),'LineWidth',1);
legend ģ���˲��� �����˲���
set(gca,'Fontname','Times New Roman','fontsize',14);


%% �����˲�����ȥ�Ʊ��˲�����������
clear all;
f = 0.01:0.01:250;%%����500�ͺܹ֡�
fs = 500;	% ����Ĳ���Ƶ��Ϊ500Hz
z_1 = exp( - 1j*2*pi*f/fs);%%f����
% ��������
Omega_1 = (10^5)/(2^17);%%��ֹƵ��
figure1 = figure('Color',[1 1 1]);

%% 50km/h
v = 50/3.6;
T = 1/fs;
b = [ Omega_1 * T , 0];
a = [ 1 + Omega_1*T ,  (-1)];
% ��ͼ
z_1 = exp( - 1j*2*pi*f/fs);%%f����
Bz = Omega_1 * T./ ( 1 + Omega_1*T - z_1 );
semilogx(v./f , 20*log10(abs(Bz)),'LineWidth',1);
T = 0.25/v;
z_1 = exp( - 1j*2*pi* f/v *0.25);
Cz = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v./f , 20*log10(abs(Cz)));
Bz1 = Bz;
Cz1 = Cz;


%% 120km/h
v = 120/3.6;
T = 1/fs;
b = [ Omega_1 * T , 0];
a = [ 1 + Omega_1*T ,  (-1)];
% ��ͼ
z_1 = exp( - 1j*2*pi*f/fs);%%f����
Bz = Omega_1*T./ ( 1 + Omega_1*T - z_1 );
hold on
semilogx(v./f , 20*log10(abs(Bz)));
T = 0.25/v;
z_1 = exp( - 1j*2*pi* f/v *0.25);
Cz = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v./f , 20*log10(abs(Cz)));

Bz2 = Bz;
Cz2 = Cz;

%% 350km/h
v = 350/3.6;
T = 1/fs;
b = [ Omega_1 * T , 0];
a = [ 1 + Omega_1*T ,  (-1)];
% ��ͼ
z_1 = exp( - 1j*2*pi*f/fs);%%f����
Bz = Omega_1*T./ ( 1 + Omega_1*T - z_1 );
semilogx(v./f , 20*log10(abs(Bz)));
T = 0.25/v;
z_1 = exp( - 1j*2*pi* f/v * 0.25);
Cz = ( exp(Omega_1*T/2) - exp(- Omega_1*T/2)*z_1 )/(Omega_1*T);
hold on;semilogx(v./f , 20*log10(abs(Cz)));
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
legend 50km/h 50km/h 120km/h 120km/h 350km/h 350km/h
% legend 50km/h 120km/h 350km/h
grid on;

Bz3 = Bz;
Cz3 = Cz;

%%
figure1 = figure('Color',[1 1 1]);
semilogx(50/3.6./f , 20*log10(abs(Bz1.*Cz1)));
hold on;
semilogx(120/3.6./f , 20*log10(abs(Bz2.*Cz2)));
semilogx(350/3.6./f , 20*log10(abs(Bz3.*Cz3)));
grid on;
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('\lambda /m')
ylabel('Mag /dB')
title('������˲���+ȥ�Ʊ��˲���')
legend 50km/h 120km/h 350km/h


%% ģ���˲���+ȥ�Ʊ��˲���
