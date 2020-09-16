clear all
clc
close all
%step1 ����GJ-4�͹�쳵������������������Ƶ����
%step2 ���δ� �������� ���Ż���ƺ�����Ƶ����
%step3 �������� ���Ż���� �Ա� ԭGJ-4�ʹ�������Ƶ����
%step4 ȷ�ϴ�����ѡ��
lamda = 1:0.001:1000;    %�ռ䲨��
pesi  = 1 ./ lamda;      %�ռ�Ƶ��

%% ���δ�1
L = 38;
N = 2*L + 1;
V_z_1  = ( (exp(-1j.*pesi)).^L - (exp(-1j.*pesi)).^(-L-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / N ;

%%���δ�2
L = 58;
N = 2*L + 1;
V_z_2  = ( (exp(-1j.*pesi)).^L - (exp(-1j.*pesi)).^(-L-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / N ;

V_z_3 = V_z_1.*V_z_2;

figure
semilogx( pesi, 20*log10(V_z_3) );

%% ���ײ��ȳ����δ�����  55 69 83   
L_base = 97;

L1   = ( L_base*0.6505 - 1 )/2 ;
% L1   = 28;
R_54 = ( (exp(-1j.*pesi)).^L1 - (exp(-1j.*pesi)).^(-L1-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / (2*L1) ;

L2   = ( L_base*0.8313 - 1 )/2 ;
% L2   = 40;
R_69 = ( (exp(-1j.*pesi)).^L2 - (exp(-1j.*pesi)).^(-L2-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / (2*L2);

L3   = ( L_base- 1 ) /2 ;
% L3   = 48;
R_83 = ( (exp(-1j.*pesi)).^L3 - (exp(-1j.*pesi)).^(-L3-1) ) ...
    ./ (  1 - exp(-1j.*pesi) ) / (2*L3);
R_z_3 = R_54 .* R_69 .*R_83 ;

figure
semilogx( pesi, 20*log10(V_z_3) ,'r', pesi, 20*log10(R_z_3),'g');
legend  ('��������','�����������Ż�');
