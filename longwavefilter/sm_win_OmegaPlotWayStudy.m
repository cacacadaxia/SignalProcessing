
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
%--------------------------------------------------------------------------


%%
clear all;
close all;
lamda = 1:0.1:1000;
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
N1 = 111;
N2 = 29;
N3 = 111;
N4 = 193;
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
ensco1 = 1 - MulFilter;         %%���Ϊʲô��������
Acc_TriRec_70 = ensco1;


%%

figure1 = figure('Color',[1 1 1]);
semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,350);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);

