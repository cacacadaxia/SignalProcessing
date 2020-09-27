% =========================================================================
%
%                  �ߵͼ����Ŀ�ĸ���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��26��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�ߵͽ��������������Թ��̽��м�
%        2.
%        3. 
%--------------------------------------------------------------------------

close all;
clear all;
addpath('func');
start_pos = 1;
N = 10000;
filepath = 'data/0916_1337_x/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
% gpxbr,hfcra,lfcrp
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
tmp2 = tmp2(start_pos:start_pos+N-1,:);
% 
tmp3 = textread([filepath , 'gppro.txt']);
tmp3 = tmp3(start_pos:start_pos+N-1,:);

%% ��ȡ�ߵ���ص�����
gpvlo = fmctrl_data(:, 8);
gpvro = fmctrl_data(:,10);
gppl  = fmctrl_data(:,7);   %%z������ٶ�
TBS   = fmctrl_data(:,end);
sita_b = tmp2(:,1);
%% �źŴ�����ò��ֵ��ǰ�����Ȳ���
% ------------1. �����趨--------------------------
% ����ֵ��������ʲô��
canMagpar = -159900;
ararm_par = 2.9280e-5 * abs(canMagpar);
pgscale = 1.9646e-19 * abs(canMagpar);
% ------------2. ������------------------------
for i = 3:length(gpvlo)
    tbs_ = TBS(i);
    gpvlo(i) = gpvlo(i)*2;
    gpvro(i) = gpvro(i)*2;
    gpvlo_dot2(i,1) = gpvlo(i) - 2*gpvlo(i-1) + gpvlo(i-2);
    gpvro_dot2(i,1) = gpvro(i) - 2*gpvro(i-1) + gpvro(i-2);
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
    az_Gz(i,1) = G(gppl(i) , tbs_);%%���ٶ��źŴ���
end
amarm = sita_b_dot2.*ararm_par;
pgrav = sita_b.^2.*TBS.^2*pgscale;
plmct = -(az_Gz + pgrav + amarm);
prmct = -(az_Gz + pgrav - amarm);
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

%% �����˲�������
zL = longwave_filter(pmcol);


%% ����Ա�
pmcol_ref = tmp3(:,1);
pmcol_ref70m = tmp3(:,2);
figure;plot(pmcol_ref,'k','LineWidth',0.5);hold on;plot(pmcol,'r--','LineWidth',0.5)
figure;plot(pmcol_ref - pmcol);
%% ����Ա�2
figure;plot(zL);hold on;plot(pmcol_ref70m);

%% ����
function out = G(x_k,tbs_k)
persistent x y1 tbs;
if isempty(x)
%     y = zeros(2,1);
    x = zeros(3,1);
    y1 = zeros(2,1);
    tbs = zeros(2,1);
end
%% ����kʱ��
x(3) = x_k;
tbs(2) = tbs_k;

%% ��ɢ�˲�
x_dot = x(3) - x(2);
x_dot_2 = x(3) - 2*x(2) + x(1);
y1(2) = (x(3) + x(2)) *tbs(2)/2^15 + x_dot;
y = (y1(2)+y1(1))* (tbs(2) + tbs(1))/2^16 + x_dot_2;


%% ����temp��
y1(1) = y1(2);
x(1) = x(2);
x(2) = x(3);
tbs(1) = tbs(2);
out = y;

end

function out = find_index(in)
%% �����˲����Ķ�����
Num = 2048;
if mod(in,Num) == 0
    out = Num;
else
    out = mod(in,Num);
end
end

function out = point3filter(in)
%% �򵥵������˲�
persistent x;
if isempty(x)
    x = zeros(3,1);
end
x(3) = in;
out = sum(x)/3;
%%
x(1) = x(2);
x(2) = x(3);

end




