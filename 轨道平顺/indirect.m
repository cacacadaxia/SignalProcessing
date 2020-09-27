function [w]=indirect()   %函数名
l=1;   %波长下限
v=150/3.6;   %车速
fs=2*v/l;   %采样频率
f=fopen('0412.txt','r');   %打开04年12月轨道高低不平顺实测数据
y=fscanf(f,'%e');
%读取原始数据
nfft=2^15;    %FFT的点数
R=xcorr(y,'unbiased');     %自相关函数
F=fft(R,nfft);     %快速傅立叶变换
p=abs(F);     %求解功率谱密度
index=0:round(nfft/2-1);     %时间坐标取一半 
k=index*fs/nfft;     %对应的频率
figure      %新窗口
w=loglog(k,p(index+1));     %画频率和功率谱密度的双对数坐标图
title('自相关法','fontsize',25)    %标题名,字号为25
xlabel('频率(Hz)','fontsize',20);     %x轴名称
grid on      %加网格
ylabel('功率谱密度(mm2)','fontsize',20);    %y轴名称

