% =========================================================================
%
%                  ��������ȡ�����ٶȼƲ�������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��26��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�ߵͽ��������������Թ��̽��м�
%        2.�ߵͻ�����������ֵûɶ����
% 10��
%        3. ����ƽ������
%--------------------------------------------------------------------------

%%
close all;
clear all;
addpath('func');
start_pos = 1;
N = 2e4;
filepath = 'data/0916_1337_x/';
load_txt;
size(wave_out);
N = length(fmctrl_data);
x = 0:0.25:0.25*(N-1);
x = x/1000;
% gpxbr,hfcra,lfcrp
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
tmp2 = tmp2(start_pos:start_pos+N-1,:);
% pmcol\gplpe\gplpe1(70m)
tmp3 = textread([filepath , 'LongWaveResultForGppro_L.txt']);
tmp3 = tmp3(start_pos:start_pos+N-1,:);
% 
tmp4 = textread([filepath , 'fmctrl_data_1337.txt']);
tmp4 = tmp4(start_pos:start_pos+N-1,:);

%% ��ȡ�ߵ���ص�����
gypitch = fmctrl_data(:,2);
gpvlo = fmctrl_data(:, 8);
gpvro = fmctrl_data(:,10);
gppl  = fmctrl_data(:,7);   %%z������ٶ�
TBS   = fmctrl_data(:,end);
sita_b = tmp2(:,1);

%% �źŴ�����ò��ֵ��ǰ�����Ȳ���
% ------------1. �����趨--------------------------
% ����ֵ��������ʲô��������˺�����
% ������Ҫ��ô�죿
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
    
end
for i = 1:length(gppl)
    az_Gz(i,1) = G( gppl(i)  , TBS(i) );
end
az_Gz = quzheng( az_Gz );
amarm = sita_b_dot2.*ararm_par;
pgrav = sita_b.^2.*TBS.^2*pgscale;
% plmct = -(az_Gz + pgrav + amarm);
% prmct = -(az_Gz + pgrav - amarm);
disp(mean(abs(amarm)))
disp(mean(abs(az_Gz)))
plmct = -(az_Gz - pgrav + amarm);%%�޸�֮��Ҳ���仯
prmct = -(az_Gz - pgrav - amarm);

% plmct = quzheng(plmct);
% prmct = quzheng(prmct);
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

pmcol = pmcol(300:end);
pmcor = pmcor(300:end);

%% ���������
pitchParameter = 2.0970;
sampleDistance = 0.25;
compf = -1;
pitch = gypitch;
for i = 1:length(pitch)
    pitch(i) = C(pitch(i) , TBS(i));%%Tbs��������C��
end
temp = sampleDistance * pitch * compf * pitchParameter;%%�����Ȼ�ǶԵ�
% temp = 0;
plmct = - ( temp + amarm );
prmct = - ( temp - amarm );

% ����Ķ���һ����
pmcol = plmct + gpvlo_dot2;
pmcor = prmct + gpvro_dot2;

%% �̲��˲�
zL_ref30m = tmp3(:,2);
zL_ref70m = tmp3(:,3);
zL_30m = shortwave_filter(pmcol);

%% �����˲�������
zL_70m = longwave_filter(pmcol,281,71,281,491);
zR_70m = longwave_filter(pmcor,281,71,281,491);
zL_250m = longwave_filter(pmcol,1123,307,1123,1953);
zR_250m = longwave_filter(pmcor,1123,307,1123,1953);

% figure1 = figure('Color',[1 1 1]);plot(zL_250m);hold on;plot(zR_250m);set(gca,'Fontname','Times New Roman','fontsize',16);

%% fdatool ��Ƹ�ͨ�˲���
b = load('filter1.mat');%%70m�ĳ���
% b = load('filter_250m.mat');%%250m����
% b = load('filter_300m.mat');%%300m����
% b = load('filter_400m.mat');%%400m����
% b = load('filter_1000m.mat');%%1000m����


tmpN = (length(b.Num)-1)/2;
pmcol_70m_tmp = conv(b.Num,pmcol);
pmcol_70m_tmp(1:tmpN) = [];pmcol_70m_tmp(end-tmpN+1:end) = [];%%ȥ����̬��֮��
z_dot = -9.6;   %%70m��������
% z_dot = -19.7;%%250m
% z_dot = -4.7;%%300m
% z_dot = -52.7;%%400m
% z_dot = -261;%%250m ���ٶȼ�
% z_dot = -132;
zL_70m_fdatool = 0;
for i = 1:length(pmcol_70m_tmp)
    z_dot = z_dot + pmcol_70m_tmp(i);
    zL_70m_fdatool = zL_70m_fdatool + z_dot;
    z_dot_save(i) = z_dot;
    zL_70m_fdatool_save(i) = zL_70m_fdatool;
end

% figure1 = figure('Color',[1 1 1]);plot( -zL_70m_fdatool_save/5);hold on;plot(zL_70m(247:end));
% legend fdatool��Ƶ��˲��� �������˲���(70m);set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / (32768/10 inch)');
figure1 = figure('Color',[1 1 1]);plot( -zL_70m_fdatool_save/5/129.01);hold on;plot(zL_70m(247:end)/129.01);
legend fdatool��Ƶ��˲��� �������˲���(70m);set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / mm');
%%����˲����Ǵ���Ư�Ƶģ��ô���������������С�����Ҳ������źŵĻ����Ư�ƣ�������
%%����ķ�ֵ�����仯�ͺ���֣�Ϊʲô��Է�ֵ����Ӱ���أ�
%%����Ҳ�Ƕ�������ײ�ֽ����˲�
%%������Ҫע��ĵض����ǣ�zL_70m(247:end)���ǲο��ĵ��λ�ã��������Ҫ

%% ���ݷ�����fdatool�Ա�
% figure;plot(zL_250m(1000:end));hold on;plot(-zL_70m_fdatool_save/5);legend 1 2
% figure;plot(zL_ref30m/129.01);

%%
% figure1 = figure('Color',[1 1 1]);plot(zL_30m);hold on;plot(zL_70m(89:end));legend ������ȡ������30m ������ȡ������70m;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / (32768/10 inch)');grid on;

%% ����Ա�
pmcol_ref = tmp3(:,1);
zL_ref30m = tmp3(:,2);
zL_ref70m = tmp3(:,3);

% figure;plot(pmcol_ref,'k','LineWidth',0.5);hold on;plot(pmcol,'r--','LineWidth',0.5)
% figure;plot(pmcol_ref - pmcol);

%% Ϊʲô����һ��ϵ���أ��о�����֣����ϵ������ʲô�ģ�
%%��Ϊ����������������������
% ratio = 1; figure1 = figure('Color',[1 1 1]);plot(ratio.*zL_30m);hold on;plot(zL_ref30m(1:end));legend ������ȡ������ ԭ����;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / (32768/10 inch)');title('30m')
% ratio = 1; figure1 = figure('Color',[1 1 1]);plot(ratio.*zL_70m);hold on;plot(zL_ref70m(177:end));legend ������ȡ������ ԭ����;set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / (32768/10 inch)');title('70m')

%% �鿴Ƶ�ף���ô�о��е�����
plot_mag(yL_ref30m,'25m')
plot_mag(yL_ref70m,'70m','hold')
legend 1 2

%% 300m����������1200����(����������һ��300m���ţ�����ʵ�����������ǲ����Ǻܺ�����Ϊ
% �����ܱ�֤һ��ʼ�ľ����Ŷյ�λ�ã�ֻ�ܼ���)
% % 300m-->1200
% start = 4000;
% z_dot = 0;
% z = 0;
% len = 300*4;
% pmcolBridge = pmcol(start+1:start+len);
% pmcolBridge = pmcolBridge - mean(pmcolBridge);
% for i = 1:len
%     z_dot = z_dot + pmcolBridge(i);
%     z = z + z_dot;
%     z_save(i) = z;
%     z_dot_save2(i) = z_dot;
% end
% 
% figure1 = figure('Color',[1 1 1]);grid on;
% subplot(2,1,1);plot(z_dot_save2);title('һ�ײ��');
% xlabel('������ /0.25m');ylabel('int')
% set(gca,'Fontname','Times New Roman','fontsize',14);
% % subplot(2,1,2);
% figure;plot(z_save);title('�ߵͲ�ƽ˳');xlabel('������ /0.25m');ylabel('int');set(gca,'Fontname','Times New Roman','fontsize',14);
% % figure1 = figure('Color',[1 1 1]);plot(zL_ref70m(start:start+len-1));hold on;plot(zL_ref30m(start:start+len-1));set(gca,'Fontname','Times New Roman','fontsize',14);grid on;
% % xlabel('������ /0.25m');ylabel('��ֵ /32768/10 inch');legend 70m 30m
% 
% % max_ = z_save(end);
% % del = max_/len;
% % for i = 1:len
% %     z_save(i) = z_save(i) - del*i;
% % end
% % figure1 = figure('Color',[1 1 1]);
% % plot(z_save);hold on;plot(yL_ref70m(start:start+len-1));
% % xlabel('������ /0.25m');ylabel('��ֵ /32768/10 inch');set(gca,'Fontname','Times New Roman','fontsize',16);





% plot_mag()
%% ��ͼ�ĺ���
function plot_for_paper()
% figure1 = figure('Color',[1 1 1]);plot( -zL_70m_fdatool_save/5/129.01);hold on;plot(zL_70m(247:end)/129.01);
% legend fdatool��Ƶ��˲��� �������˲���(70m);set(gca,'Fontname','Times New Roman','fontsize',16);xlabel('��� /0.25m');ylabel('�ߵ� / mm');
end
%% ����
function out = C(x_k,tbs)
omega1 = 0.76;

persistent x;
if isempty(x)
    x = zeros(2,1);
end
%%
x(2) = x_k;
y = x(2)-x(1) + tbs/2^18*(x(2)+x(1));
y = y/omega1;
%%
x(1) = x(2);
out = y;
end


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

% out = y/tbs(1)/tbs(2);
out = y;
tbs(1) = tbs(2);


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


function out = quzheng(in)
out = zeros(2,1);
for i = 1:length(in)
    if in(i)>=0
        out(i) = floor(in(i));
    elseif in(i)<0
        out(i) = floor(in(i)) + 1;
    end
end
end


function plot_mag(signal_data , tit , varargin)

if (nargin == 3)
    mode = varargin{1};
    if mode == 'hold'
        hold on;
    end
else
    figure;
end

fs = 4;     %% 0.25mΪһ���������
N = length(signal_data);
x = (1:N/2+1)/N*fs;
x = 1./x;
tp = abs(fftshift(fft(signal_data)));
tp = 20*log10(tp);
semilogx( x , tp((length(tp)/2):end)   );
xlabel('\lambda m')
ylabel('Mag dB')
set(gca,'Fontname','Times New Roman','fontsize',16);
title(tit);
grid on;
end


function out = B1(x,tbs)
persistent y;
if isempty(y)
    y = zeros(2,1);
end
Omega1 = 10^5/2^17;
y(2) = ( y(1) *2^17 + tbs*x )/(2^17 + tbs);

%%
y(1) = y(2);
out = y(2);
end
%% �̲��˲���
function out = shortwave_filter(in)
amcol = in;
%% 25m�����˲��ӻ���
alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
sscal = 0.000825;
sbsci = 0.019802;
fscal = 0.1;
sbsc = 101.000;
% �����趨
Num = 768;
amcol_array = zeros(Num,1);
amcol_arraytmp = zeros(Num,1);
in1 = 533; in2 = 432;in4=382;in6=331;in7=230;
in = 539;

for i = 1:length(amcol) 
    amcol_array(in) = amcol(i);
    alu = alu + amcol_array(in1) - 3*amcol_array(in2) + 3*amcol_array(in6) - amcol_array(in7);
    elupp = alu;
    elup = elup + elupp;
    elu = elu + elup;
    emco = - amcol_array(in4);
    als = als + amcol_array(in2) - amcol_array(in6);
    alss = alss + als;
    alss = alss + sbsc*emco;
    alsss = alsss + alss;
    xtemp = (alsss*sbsci - sscal*elu)*fscal;
    yL(i,1) = xtemp;
    
    %% ��������
    in1 = mod(in1,Num)+1;
    in2 = mod(in2,Num)+1;
    in4 = mod(in4,Num)+1;
    in6 = mod(in6,Num)+1;
    in7 = mod(in7,Num)+1;
    in = mod(in,Num)+1;
    %%
    save(i,1) = alu;
    save(i,2) = elu;
    save(i,3) = als;
    save(i,4) = alsss;
end
%%
out = yL;
end
