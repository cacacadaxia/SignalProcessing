
% =========================================================================
%
%                  �����˲�����ͼ
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��27��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------
% -----------------------
% 
%   1. ���޸ĵ�ϵ�����м���
% �����˸��ֲ�ͬ�Ĳ���
% 
% ------------------------

clear all;
close all;
% ����������
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
omega = 2*pi*pesi*0.25;
N = length(lamda);

% ��������
para = winLen_para_set();%%�����趨
% para(3,:) = [160*2,40*2,169,280];%%word�ĵ��е��˲�������
% ��ȡ�˲�������
[Acc_TriRec_25, TriDB_TriRec_25, maxV_TriRec_25, minV_TriRec_25, BW_TriRec_25, slope_TriRec_25] = tri_rec_paral_repair(25, para);
[Acc_TriRec_42, TriDB_TriRec_42, maxV_TriRec_42, minV_TriRec_42, BW_TriRec_42, slope_TriRec_42] = tri_rec_paral_repair(42, para);
[Acc_TriRec_70, TriDB_TriRec_70, maxV_TriRec_70, minV_TriRec_70, BW_TriRec_70, slope_TriRec_70] = tri_rec_paral_repair(70, para);
[Acc_TriRec_120, TriDB_TriRec_120, maxV_TriRec_120, minV_TriRec_120, BW_TriRec_120, slope_TriRec_120] = tri_rec_paral_repair(120, para);

% plot 25

figure1 = figure('Color',[1 1 1]);
semilogx((lamda), ((Acc_TriRec_25)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_25,lamda,25);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
% plot 42
figure1 = figure('Color',[1 1 1]);
semilogx((lamda), ((Acc_TriRec_42)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_42,lamda,42);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
% plot 70
figure1 = figure('Color',[1 1 1]);
semilogx((lamda), ((Acc_TriRec_70)), 'b','LineWidth',2);
plot_line_func(Acc_TriRec_70,lamda,70);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);
% plot 120
figure1 = figure('Color',[1 1 1]);
semilogx((lamda), (20*log10(Acc_TriRec_120)), 'b','LineWidth',2);
% plot_line_func(Acc_TriRec_120,lamda,120);
xlabel('Wavelength /m')
ylabel('Magnitude /dB')
set(gca,'Fontname','Times New Roman','fontsize',14);


fprintf('ȫ���������\n')
fprintf('\tN_Tri1\tN_Tri2\tN_Rec1\tN_Rec2\n');
fprintf("-------------------------------------------\n")
list = {'25m','42m','70m','120m'};
for i = 1:4
    fprintf("%s\t%d\t%d\t%d\t%d\t\n",list{i},para(i,:));
end
function [ensco1, tridb, maxV, minV, BW, slope] = tri_rec_paral_repair(lam, para)
%���Ǵ�+���δ�������Ƶ���Է���
% N=128^2;
% omega=linspace(0,pi,N);
lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
omega = 2*pi*pesi*0.25;
N = length(lamda);

Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);

switch (lam)
    case 25
        N1 = para(1,1);
        N2 = para(1,2);
        N3 = para(1,3);
        N4 = para(1,4);
    case 42
        N1 = para(2,1);
        N2 = para(2,2);
        N3 = para(2,3);
        N4 = para(2,4);
    case 70
        N1 = para(3,1);
        N2 = para(3,2);
        N3 = para(3,3);
        N4 = para(3,4);
    case 120
        N1 = para(4,1);
        N2 = para(4,2);
        N3 = para(4,3);
        N4 = para(4,4);
end
        
%% way 2 (������ԭ�е�����)

% m1 = floor(lam*2.4+1.9+0.5);
% M1 = m1 + 1 - mod(m1,2);
% m2 = floor(M1*0.3216 + 0.5);
% M2 = m2 + 1 - mod(m2,2);
% M3 = round(M1*1.1228)*2+1;
% M4 = round(M1*1.5439)*2+1;

% ���δ��Ĺ�����ʽ
    M1 = N1/2;
    M2 = N2/2;
for i=1:N
%     ���Ǵ�
    Tri1(i)= ((sin(omega(i)*M1/2)/sin(omega(i)/2))/M1)*((sin(omega(i)*M1/2)/sin(omega(i)/2))/M1);
    Tri2(i)= ((sin(omega(i)*M2/2)/sin(omega(i)/2))/M2)*((sin(omega(i)*M2/2)/sin(omega(i)/2))/M2);
    
%     ���δ�
    Rec1(i)= ((sin(omega(i)*N3/2)/sin(omega(i)/2))/N3);
    Rec2(i)= ((sin(omega(i)*N4/2)/sin(omega(i)/2))/N4);
end
% ������ϵ��
MulFilter = 1.036*Tri1-0.036*Tri2+0.25*(Rec1-Rec2);
% a = 1.0;
% b = 1.0;
% MulFilter = Tri1+a*0.036*(Tri1-Tri2)+b*0.25*(Rec1-Rec2);

% test
% MulFilter = Tri1-Tri2+Rec1-Rec2;
ensco1 = 1 - MulFilter;         %%���Ϊʲô��������

% ensco1 = MulFilter;

Acc = ensco1;
[tridb, maxV, minV, BW, slope] = filter_para(Acc);
end

function para = winLen_para_set()
list = [25,42,70,120];
% list = [25,42,70,250];
% list = [25,42,70,3];
for i = 1:4
    Wavelen = list(i);
    [Tri1,Tri2,Rec1,Rec2] = winlen_lamda(Wavelen);
    para(i,:) = [Tri1,Tri2,Rec1,Rec2];
end
end

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

function [triDB, maxV,minV, BW, slope] = filter_para(F)
%triDB: -3db
%maxV: ͨ����󲨶�
%BW������
%slope��3�ֱ���б��
ensco1 = F;

lamda = 1:0.1:1000;
pesi  = 1./lamda;           %�ռ���Ƶ��
omega = 2*pi*pesi*0.25;
N = length(lamda);

cutn = 3;
ensco = ensco1([cutn:end]);
lamda1 = lamda([cutn:end]);
[maxV1,index] = max(ensco);

maxL = 0;
triDB = 0;
minL = 0;
tol = 0.02;
index = abs(ensco-0.707) < tol;
indexNumeric = find(index); 
indexM = floor(median(indexNumeric));

if ~isnan(indexM) 
    triDB = (lamda1(indexM));
    dx=diff(lamda1);
    dy=diff(ensco);
    dydx = dy./dx;
    slope = dydx(indexM);
else
    slope = NaN;
end

% index = abs(ensco-maxV1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));

index = abs(ensco - 1.0) < tol;
indexM = find(index,1,'first');
if ~isnan(indexM) 
    maxL = (lamda1(indexM));
    
end
tmp1 = indexM;

ensco2 = ensco([indexM:end]);
maxV = max(abs(ensco2-1));

% index = abs(ensco-0.1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));

index = abs(ensco - 0.05) < tol;
indexM = find(index,1,'last');
if ~isnan(indexM) 
    minL = (lamda1(indexM));
end
ensco3 = ensco([1:indexM]);
tmp2 = indexM;
minV = max(abs(ensco3-0));


PLOT = 0;
if PLOT
    fprintf('Results of transition zone:\n')
    disp([lamda1(tmp1) , lamda1(tmp2)]);
    disp([(tmp1) , (tmp2)]);
end

BW = minL - maxL;
end