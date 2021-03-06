%%%%%%%%%%%%
clear all;
clc;
wp = 0.2 * pi;
ws = 0.3 *pi;
rp = 3;
rs = 20;

[n,wn]=buttord(wp,ws,rp,rs,'s');
[bn,an] = butter(n,wn,'s');
fprintf('巴特沃斯滤波器的阶数 N = %4d \n',n);

[z,p,k] = buttap(n);
[Bap,Aap] = zp2tf(z,p,k);
[bb,ab] = lp2lp(Bap,Aap,wn);

[Hn,fn] = freqs (bn,an);
[Hb,wb] =freqs(bb,ab);

plot ((wn/pi) ,20*log10(abs(Hn)),'r','linewidth',3);hold on;
plot ((wb/pi) ,20*log10(abs(Hb)),'k')
grid;axis([0 500 -60 5]);hold on;
%plot(fz,20 *log10(abs(Hz)),'k');
legend('模拟滤波器','数字滤波器');
xlabel('频率/Hz');ylabel('幅值/db');
title('不进行预畸的数字滤波器与模拟滤波器响应曲线比较');
set(gcf,'color','w');box on;

%{
fs=30;              %采样频率
N=300;              %N/fs 秒数据
n=0:N-1;            
t=n/fs;             %时间

if 0
fl = 0.4;           %低频
fh = 5;             %高频
s=cos(2*pi*fl*t)+cos(2*pi*fh*t);    %s是0.4Hz和5Hz信号叠加，低通截止频率是1Hz
subplot(121);plot(t,s);
title('输入信号');xlabel('t/s');ylabel('幅度');
sfft=fft(s);
subplot(122);
plot((1:length(sfft)/2)*fs/length(sfft),2*abs(sfft(1:length(sfft)/2))/length(sfft));
title('信号频谱');xlabel('频率/Hz');ylabel('幅度');
%设计低通滤波器，截止频率为1
if 0
Wp=1/fs;Ws=2/fs;                %截止频率为1Hz,阻带截止频率为2Hz
%估算得到Butterworth低通滤波器的最小阶数N和3dB截止频率Wn
[n,Wn]=buttord(Wp,Ws,1,50);     %阻带衰减大于50db,通带纹波小于1db
else
n=4;
Wn=1/(fs/2);
end
%设计Butterworth低通滤波器
[a,b]=butter(n,Wn);
[h,f]=freqz(a,b,'whole',fs);        %求数字低通滤波器的频率响应
f=(0:length(f)-1*fs/length(f));     %进行对应的频率转换
figure;
plot(f(1:length(f)/2),abs(h(1:length(f)/2)));       %绘制幅频响应图
title('巴特沃斯低通滤波器');xlabel('频率/Hz');ylabel('幅度');
grid;
sF=filter(a,b,s);                   %叠加函数s经过低通滤波器以后的新函数
figure;
subplot(121);
plot(t,sF);                         %绘制叠加函数s经过低通后时域图形
title('输出信号');xlabel('t/s');ylabel('幅度');
SF=fft(sF);
subplot(122);
plot((1:length(SF)/2)*fs/length(SF),2*abs(SF(1:length(SF)/2))/length(SF));
title('低通滤波后频谱');xlabel('频率/Hz');ylabel('幅度');
end

%带通滤波 
fl = 2;             %低频
fh = 10;             %高频
s=cos(2*pi*fl*t)+cos(2*pi*fh*t);    % s是2Hz和3Hz信号叠加，带通截止频率是1Hz~4Hz

H=s;
N=3636;              %N/fs 秒数据
n=0:N-1;            
t=n/fs;             %时间

subplot(121);plot(t,s);
title('输入信号');xlabel('t/s');ylabel('幅度');
sfft=fft(s);
subplot(122);
plot((1:length(sfft)/2)*fs/length(sfft),2*abs(sfft(1:length(sfft)/2))/length(sfft));
title('信号频谱');xlabel('频率/Hz');ylabel('幅度');
%设计带通滤波器，截止频率为0.4~5
n=4;
Wn=[0.4/(fs/2) 5/(fs/2)]

%设计Butterworth低通滤波器
[a,b]=butter(n,Wn);
[h,f]=freqz(a,b,'whole',fs);        %求数字低通滤波器的频率响应
f=(0:length(f)-1*fs/length(f));     %进行对应的频率转换
figure;
plot(f(1:length(f)/2),abs(h(1:length(f)/2)));       %绘制幅频响应图
title('巴特沃斯带通滤波器');xlabel('频率/Hz');ylabel('幅度');
grid;
sF=filter(a,b,s);                   %叠加函数s经过低通滤波器以后的新函数
figure; 
subplot(121);
plot(t,sF);                         %绘制叠加函数s经过低通后时域图形
title('输出信号');xlabel('t/s');ylabel('幅度');
SF=fft(sF);
subplot(122);
plot((1:length(SF)/2)*fs/length(SF),2*abs(SF(1:length(SF)/2))/length(SF));
title('带通滤波后频谱');xlabel('频率/Hz');ylabel('幅度');
%}
