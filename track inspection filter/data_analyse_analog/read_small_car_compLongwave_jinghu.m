

% =========================================================================
%
%                  �ԱȾ�̬С�����쳵������
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��31��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�����Ķ�׼�����������ݣ����ֶԱȲ��ϣ�����ǰ���ʮ�����
%        2.
%       3. 
%--------------------------------------------------------------------------


close all;
clear all;
% o = textread('data_small_car.txt');
fd = fopen('data_small_car.txt','r');
A = textscan(fd,'%s%f%f%f%f%f%f%f%f%f','Delimiter','\t','CollectOutput',1);
for i = 1:length(A{1})
    a = str2num(A{1}{i}(2:5));
    b = str2num(A{1}{i}(7:end));
    xSmallCar(i) = a + b/1e3;
end
gpproLTrolley = A{2}(:,5);%%��ߵ�
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolley);

gpproLTrolley = out(len+1:end-len);

figure1 = figure('Color',[1 1 1]);
plot(xSmallCar,gpproLTrolley);
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('��� /1km');ylabel('���� /mm')

%% �����Ƴ���0.65mת����Ϊ0.25m������ͬ��
gpproLTrolley25m = [];
for i = 1:1720*4
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolley(p1) + p2/0.65*( gpproLTrolley(p1+1) - gpproLTrolley(p1) );
end
figure;plot(0:0.25:(length(gpproLTrolley25m)-1)*0.25 , gpproLTrolley25m);
hold on;plot(0:0.65:(length(gpproLTrolley)-1)*0.65 , gpproLTrolley,'--k');


%%
% load_txt;

filepath = 'data/1031-1130/';
start_pos = 1;
N = 1e6;
load_txt;
x = 0:0.25:0.25*(N-1);
x = x/1000;

gpvmp = output_wave(end:-1:1,1);%%km
gpvft = output_wave(end:-1:1,2);%%
gplpe3 = output_wave(end:-1:1,43);
gprpe3 = output_wave(end:-1:1,44);
gppl = fmctrl_data(end:-1:1,11);
gppr = fmctrl_data(end:-1:1,12);
TBS = fmctrl_data(end:-1:1,20);
gplae = output_wave(end:-1:1,5);

%% һ�����ŵ�
start_mileage = 1106;  end_mileage = 1155;
index_start = find_index(gpvmp,gpvft,start_mileage);
index_end = find_index(gpvmp,gpvft,end_mileage);

% index_start = 1;
% index_end = 4e4;

x = gpvmp(index_start:index_end) + gpvft(index_start:index_end)/4000;
figure1 = figure('Color',[1 1 1]);
hold on;
tag1 = 0;
plot( x , gplpe3(index_start + tag1:index_end + tag1)/129.01);xlabel('��� /1km');ylabel('���� /mm')
set(gca,'Fontname','Times New Roman','fontsize',14);
legend 1 2;

% figure1 = figure('Color',[1 1 1]);plot(x,gplae(index_start:index_end)/129.01);xlabel('��� /1km');ylabel('���� /mm')
% set(gca,'Fontname','Times New Roman','fontsize',16);

%% ����һ��Ѱ�����ݶ���ķ���
gplpe3 = gplpe3(index_start:index_end)/129.01;
for i = 1:length(gplpe3)-length(gpproLTrolley25m)-1
    res(i) = sum(abs(gplpe3(i:i+length(gpproLTrolley25m)-1)' - gpproLTrolley25m));
end
figure;plot(res);
figure;plot(gplpe3);hold on;plot([zeros(1,77111) , gpproLTrolley25m])
% start_mileage + 23408/4000


%%
l1 = gpproLTrolley25m;
out = conv(gplpe3 , l1(end:-1:1));
[~,index] = max(out);
index = index - length(l1);
figure;plot(out);



%% ����
function index = find_index(gpvmp,gpvft,pos)
par1 = floor(pos);
par2 = pos - floor(pos);
in3 = floor(par2*4000);
in1 = find(gpvmp == par1);
in2 = gpvft(in1(1));
index = in1(1) + in3 - in2;
if index<=0
    index = 1;
end
end

