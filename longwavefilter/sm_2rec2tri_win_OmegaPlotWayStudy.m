
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
%        2.
%        3. 
%--------------------------------------------------------------------------


clear all;
close all;
% ����������
N=128^2;
w=linspace(0,pi,N);
lamda = 2*pi./w./4;%%������ʾ���������

%%  ��������
%���Ǵ�+���δ�������Ƶ���Է���
N=128^2;
Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);  
%% way 2 (������ԭ�е�����)
M1 = 160;
M2 = 40;
M3 = 160;
M4 = 280;
for i=1:N
%     ���Ǵ�
    Tri1(i)= ((sin(w(i)*M1/2)/sin(w(i)/2))/M1)*((sin(w(i)*M1/2)/sin(w(i)/2))/M1);
    Tri2(i)= ((sin(w(i)*M2/2)/sin(w(i)/2))/M2)*((sin(w(i)*M2/2)/sin(w(i)/2))/M2);
%     ���δ�
    Rec1(i)= ((sin(w(i)*M3/2)/sin(w(i)/2))/M3);
    Rec2(i)= ((sin(w(i)*M4/2)/sin(w(i)/2))/M4);
end
% ������ϵ��
MulFilter = 1.036*Tri1-0.036*Tri2+0.25*(Rec1-Rec2);
ensco1 = 1 - MulFilter;         %%���Ϊʲô��������
Acc_TriRec_70 = ensco1;

%%
fig = figure;

semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,70);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);


function [Tri1,Tri2,Rec1,Rec2] = winlen_lamda(lamda)

Tri1 = 4.4940*lamda-0.2418; Tri1 = CalWinLen(Tri1);
Tri2 = 1.2341*lamda-0.7893; Tri2 = CalWinLen(Tri2);
Rec1 = 4.4940*lamda-0.2418; Rec1 = CalWinLen(Rec1);
Rec2 = 7.8231*lamda-1.6313; Rec2 = CalWinLen(Rec2);

out = [Tri1,Tri2,Rec1,Rec2];
end

function Len = CalWinLen(tmp)
if mod(floor(tmp),2) == 1
    Len = floor(tmp);
else
    Len = floor(tmp)-1;
end
end