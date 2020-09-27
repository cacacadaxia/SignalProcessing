
% =========================================================================
%
%                  ���ֹ�������㷨����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.���ȸ�����Ϊ�򵥵Ĺ��򲿷�
%        2.�����˲����Ĳ��֣������жԱȡ���Ҫ�ǹ��򣬽���ܺ�
%        3. ��ʼ�����ǲ�ֵĳ�ʼ������������������ûʲô
% 
% 
%--------------------------------------------------------------------------

% load_txt;
close all;
clear all;
filepath = 'data/0916_1337_x/';
start_pos = 1;
N = 10000;
load_txt;
size(wave_out);
x = 0:0.25:0.25*(N-1);
x = x/1000;
%%
tmp5 = textread([filepath,'tmp2.txt']);
if length(tmp5)>N
    tmp5 = tmp5(start_pos:start_pos+N-1,:);
end
%%
tmp2 = textread([filepath,'tmp_zhongjian_1337.txt']);
if length(tmp2)>N
    tmp2 = tmp2(start_pos:start_pos+N-1,:);
end
% gpxbr,hfcra,lfcrp

%% ���ĶԱ�
chaogao = wave_out(:,6);

% ****************�����趨*********************
delay = 418; %%һֱ��������𣿲�һ����ÿ��������Ŀ��Ӧ�Ĳ�������һ��
%% ������������ɶ��˼������������׼ȷ��
%% ȷ�������������ĺ���
G_par = 3.8259e-14*141500.0;
ht = 3.90398e-05*141500.0*0.268;

tbs = fmctrl_data(:,end);
tbs_s = tbs/1e5;
% ***************step1 ģ�ʹ*************************
rou_l = fmctrl_data(:,9);
rou_r = fmctrl_data(:,11);

% rou_l = rou_l/(129.01);
% rou_r = rou_r/(129.01);
rou_l_dot2(1,1) = rou_l(1);
rou_r_dot2(1,1) = rou_r(1);
rou_l_dot2(2,1) = rou_l(2) - 2*rou_l(1);
rou_r_dot2(2,1) = rou_r(2) - 2*rou_r(1);
for i = 3:length(rou_l)
    rou_l(i) = rou_l(i)*2;
    rou_r(i) = rou_r(i)*2;
    rou_l_dot2(i,1) = rou_l(i) - 2*rou_l(i-1) + rou_l(i-2);
    rou_r_dot2(i,1) = rou_r(i) - 2*rou_r(i-1) + rou_r(i-2);
end

% ******************step2 ���ٶȼƵ��˲� **********************************
ay = fmctrl_data(:,5);

% ---------------------- �����˲��� --------------------
for i = 1:length(ay)
    ay_Gz(i,1) = G(ay(i) , tbs(i));
end

%% Ƶ�׹۲�
% plot_mag(ay,'�˲�ǰ')
% plot_mag(ay_Fz,'�˲���')
% 
% figure;plot(aln(:,2)-ay_Gz);


%% ����
% sita_b = sita_b/3276.8/180*pi;
sita_b = tmp2(:,1);
sita_b_dot2(2,1) = sita_b(2)-sita_b(1) - sita_b(1);
sita_b_dot2(1,1) = sita_b(1);
for i = 3:length(sita_b)
    sita_b_dot2(i,1) = sita_b(i) - 2*sita_b(i-1) + sita_b(i-2);
end

camo = ay_Gz - G_par .* sita_b .* tbs.^2 + ht * sita_b_dot2;
camo = -camo;
camo = quzheng(camo);%%����˵����ȡ������������
% camo = aln(:,1);%%���������ȡ����
%%����˵����camo��Ĳ�׼����Ҫ��ay�е������
amcol = camo + rou_l_dot2;
amcor = camo - rou_r_dot2;
%%
marm = aln(:,5);
figure;plot(marm - ht * sita_b_dot2);
 %%
 figure;plot( amcol - aln(:,3));
 figure;plot(camo - aln(:,1));
 %%
alu = 0;elupp = 0;elup = 0;elu = 0;als = 0;alss = 0;alsss = 0;
sscal = 0.000825;
sbsci = 0.019802;
fscal = 0.1;
sbsc = 101.000;
%% ��Ȼ�������ķ���ȥʵ���ˣ����ǲ�û�л����ͬ�Ľ��
%% ���鲻֪����ɶ�þ�����
%% ����Ļ�����ʲô�����أ�
%% �иĽ��Ŀռ�
Num = 768;
amcol_array = zeros(Num,1);
amcol_arraytmp = zeros(Num,1);
in1 = 533; in2 = 432;in4=382;in6=331;in7=230;
in = 539;

for i = 1:length(amcol)           %%�򵥻��֣��϶��ǲ��Ե�
    %%
%     amcol = aln(:,3);
   amcol_array(in) = amcol(i);
   amcol_array_tmp(i,1) = amcol_array(in1);
   amcol_array_tmp(i,2) = amcol_array(in2);
   amcol_array_tmp(i,3) = amcol_array(in4);
   amcol_array_tmp(i,4) = amcol_array(in6);
   amcol_array_tmp(i,5) = amcol_array(in7);
    %%-----------------------------
%     amcol_array(in1) = tmp5(i,1);
%     amcol_array(in2) = tmp5(i,2);
%     amcol_array(in4) = tmp5(i,3);
%     amcol_array(in6) = tmp5(i,4);
%     amcol_array(in7) = tmp5(i,5);
    %%�ɴ˿ɼ���������������⣬���Բ�������
    %%-----------------------------
    
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
% Ŀǰֻ�л����Ǵ��������

%% ����Ա�
% camco result  amcol[gpami]
% figure;plot(camo(30:end));hold on;
% plot(aln(:,1));
% legend 1 2
% figure;plot(amcol_array_tmp - tmp5(:,1:5));
% figure;plot(camo - aln(:,1));
%% alu elu als alsss
    % tmp4 = textread([filepath,'tmp1.txt']);
    % if length(tmp4)>N
    %     tmp4 = tmp4(start_pos:start_pos+N-1,:);
    % end
    % l = tmp4 - save;
    % figure;plot(l);
    % legend 1 2 3 4
    % figure;plot(l(:,1));
%%
% aln(:,4);
figure;plot(yL,'LineWidth',1);hold on;plot(aln(:,4));legend 1 2
figure;plot((yL - aln(:,4))/103);%%������ȫһ�£�����������ʱ�Ȳ��ÿ�����

%% wave out (����͸ߵ͵Ļ�ͼ)
gaodi_l = wave_out(:,1)/129.01;
guixiang_l = wave_out(:,3)/103.21;

%% ����Ĳ���Ƶ��
aln_l = wave_out(:,3);
aln_r = wave_out(:,4);

function out = F(x,tbs)
%% Fs�ı�׼д��
%% �˲����趨
%% ����x[3]������ͬ��
% x(3)=x_n
% x(2)=x_n-1
% x(1)=x_n-2
%% ����˲���������ʱ����ʮ����
persistent y;
if isempty(y)
    y = zeros(3,1);
end
y(3) = ( y(2)*(2*2^28+2^14*tbs) - y(1)*2^28+tbs^2*x  )/(2^28 + 2^14*tbs + tbs^2);
%% ����temp��
y(1) = y(2);
y(2) = y(3);
out = y(3);
end


function out = R(x_k,tbs)
wd = 0.001;

persistent y x;
if isempty(y)
    y = zeros(2,1);
    x = zeros(2,1);
end
x(2) = x_k;
x_dot = x(2)-x(1);
y(2) = (1 - wd) * (x_dot + y(1));

%% ����temp��
x(1) = x(2);
y(1) = y(2);
out = y(2);

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
tbs(1) = tbs(2);
out = y;

end

function plot_mag(signal_data , tit)
figure;
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



