

% =========================================================================
%
%                  对比静态小车与轨检车的数据(京沪线)
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 10月31日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.尽量的对准两个车的数据，发现对比不上，搜索前后二十公里的
%        2.这里的轨检车的数据有里程跳变，跳变后的里程反而不准
%       3. 京沪线重放的波形是倒着的，所以这一点需要注意
% 
% 
% 
%--------------------------------------------------------------------------



%% step1: 读取小车数据
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
gpproLTrolleyOri = A{2}(:,5);%%左高低
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolleyOri);
gpproLTrolley = out(len+1:end-len);

% 小车的数据
figure1 = figure('Color',[1 1 1]);
plot(xSmallCar,gpproLTrolleyOri);
set(gca,'Fontname','Times New Roman','fontsize',14);
xlabel('里程 /1km');ylabel('长波 /mm')
hold on;plot(xSmallCar,gpproLTrolley);

%% step3: 将手推车的0.65m转换成为0.25m，与轨检同步
gpproLTrolley25m = [];
for i = 1:1720*4
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolley(p1) + p2/0.65*( gpproLTrolley(p1+1) - gpproLTrolley(p1) );
end
% figure;plot(0:0.25:(length(gpproLTrolley25m)-1)*0.25 , gpproLTrolley25m);
% hold on;plot(0:0.65:(length(gpproLTrolley)-1)*0.65 , gpproLTrolley,'--k');

%% step4: 读取轨检车数据
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

%% step2: 截断数据(轨检车)
start_mileage = 1106;  end_mileage = 1155;
index_start = find_index(gpvmp,gpvft,start_mileage);
index_end = find_index(gpvmp,gpvft,end_mileage);

x = gpvmp(index_start:index_end) + gpvft(index_start:index_end)/4000;
index1 = find(x(2:end) - x(1:end-1)>0.001);
km = x(index1+1) - x(index1);
% x(index1+1:end) = x(index1+1:end) - km;
x(1:index1) = x(1:index1) + km;
gplpe3 = gplpe3(index_start:index_end)/129.01;

% figure1 = figure('Color',[1 1 1]);hold on;
% plot( x , gplpe3);xlabel('里程 /1km');ylabel('长波 /mm')
% set(gca,'Fontname','Times New Roman','fontsize',14);

%% step5: 利用相关性找到对齐点
l1 = gpproLTrolley25m;
out = conv(gplpe3 , l1(end:-1:1));
[~,index] = max(out);
index = index - length(l1);
figure;plot(gplpe3);hold on;plot([zeros(1,index) , gpproLTrolley25m])%%相距10km

%%


%% 直接对齐肯定是有偏差的
% figure;plot(x , gplpe3);hold on;plot(xSmallCar,gpproLTrolley);

%% 函数
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

