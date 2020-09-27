function [q]=china()     %函数名称
t=fopen('BJ.txt','r');    %打开中国三大干线参数文件
p=fscanf(t,'%e');    %读取数据
fs=1000;    %设置采样频率
l=1:1/fs:30;    %波长范围
A=p(1,1);B=p(2,1);C=p(3,1);D=p(4,1);E=p(5,1);F=p(6,1);G=p(7,1);   %给各参数赋值
S=A.*(l.^2+B*l.^3+C*l.^4)./(l.^2.*(G*l.^4+F*l.^3+E*l.^2+D*l+1));               %经单位转化后的轨道谱公式
q=loglog(l,l.^2.*S);     %画波长与功率谱密度的双对数坐标图
title('三大干线左轨高低','fontsize',25)     %题名，并设置字号为25
xlabel('波长/m','fontsize',20)     %x轴名称
grid on    %加网格
ylabel('功率谱密度/[mm2/(rad/m)]','fontsize',20)    %y坐标名称
figure
A=p(8,1);B=p(9,1);C=p(10,1);D=p(11,1);E=p(12,1);F=p(13,1);G=p(14,1);           %给各参数赋值
S=A.*(l.^2+B*l.^3+C*l.^4)./(l.^2.*(G*l.^4+F*l.^3+E*l.^2+D*l+1));
loglog(l,l.^2.*S);
title('三大干线左轨轨向','fontsize',25)
xlabel('波长/m','fontsize',20)
grid on
ylabel('功率谱密度/[mm2/(rad/m)]','fontsize',20)
figure
A=p(15,1);B=p(16,1);C=p(17,1);D=p(18,1);E=p(19,1);F=p(20,1);G=p(21,1);         %给各参数赋值
S=A.*(l.^2+B*l.^3+C*l.^4)./(l.^2.*(G*l.^4+F*l.^3+E*l.^2+D*l+1));
loglog(l,l.^2.*S);
title('三大干线水平','fontsize',25)
xlabel('波长/m','fontsize',20)
grid on
ylabel('功率谱密度/[mm2/(rad/m)]','fontsize',20)
