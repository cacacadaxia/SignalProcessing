
% =========================================================================
%
%                   �����˲���ԭ��̽��
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��23��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�ɰ汾�ĳ����˲����ģ�70m��������Ƶ���߶Ա�
%        2.��ͬ�˲����Աȣ���������new_report_2tri_2rec_final.m
%        3. ע�����Ǵ���M�����δ���N��������ȷ����
%        4.
%        5.
% 1120
%       1.���������˲��������Ͳ����ĶԱ�
%--------------------------------------------------------------------------


%%
clear all;
close all;
lamda = 1:0.1:10000;
pesi  = 1./lamda;           %�ռ���Ƶ��
omega = 2*pi*pesi*0.25;


%%  ��������
%���Ǵ�+���δ�������Ƶ���Է���
N = length(lamda);
Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);  

%% ��������
N1 = 539;
N2 = 147;
N3 = 539;
N4 = 937;
% 			1571	431	1571	2735	350m
%% ����봰��
M1 = (N1-1)/2;%%���Ǵ��İ봰��
M2 = (N2-1)/2;%%���Ǵ��İ봰��

%%
for i=1:N
    % ���Ǵ�
    Tri1(i)= ((sin(omega(i)*M1/2)/sin(omega(i)/2))/M1)*((sin(omega(i)*M1/2)/sin(omega(i)/2))/M1);
    Tri2(i)= ((sin(omega(i)*M2/2)/sin(omega(i)/2))/M2)*((sin(omega(i)*M2/2)/sin(omega(i)/2))/M2);
    % ���δ�
    Rec1(i)= ((sin(omega(i)*N3/2)/sin(omega(i)/2))/N3);
    Rec2(i)= ((sin(omega(i)*N4/2)/sin(omega(i)/2))/N4);
end
% ������ϵ��
MulFilter = 1.036 * Tri1 - 0.036*Tri2 + 0.25 * (Rec1-Rec2);

% �������˲���
MulFilter2 = Tri1.*Tri2.*Rec1.*Rec2;

ensco1 = 1 - MulFilter;
ensco2 = 1 - MulFilter2;
Acc_TriRec_70 = ensco1;


%%
figure1 = figure('Color',[1 1 1]);
% semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,120);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
hold on;
plot_line_func(abs(ensco2),lamda,120);
%% �Ա�fdatool
% b = load('filter_120m.mat');
% [h,f] = freqz(b.Num,[zeros(1,length(b.Num)-1),1],10000,4);
% % figure1 = figure('Color',[1 1 1]);
% semilogx(1./f, abs(h) ,'LineWidth','1');grid on;xlabel('���� /m');
% ylabel('��ֵ /dB')
% set(gca,'Fontname','Times New Roman','fontsize',16);
% legend ������ FIR�߽��˲���
