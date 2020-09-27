function [q1,q2,q3]=usanew()     %������
fs=1000;      %����Ƶ��
l=0.5:1/fs:300;     %������Χ
Av=0.2095;
Aa=0.0762;
Wc=0.8245;
Ws=0.8209;
k=0.25;
Lc=2*pi/Wc;
Ls=2*pi/Ws;        %��������
%�����弶��·
S1=100*k*Av./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %�ߵͲ�ƽ˳
S2=100*k*Aa./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %����ƽ˳
S3=100*2*k*Av./(pi*l.^2*Lc.^2.*(1./l.^2+1/Lc.^2).*(1./l.^2+1/Ls.^2));   %ˮƽ����಻ƽ˳
%����������·
Av2=0.0339;
Ws2=0.4380;
Ls2=2*pi/Ws2;   %��������
S4=100*k*Av2./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %�ߵͲ�ƽ˳
S5=100*k*Av2./(2*pi*(2*pi./(Wc*l)).^2+2*pi);   %����ƽ˳
S6=100*2*k*Av2./(pi*l.^2*Lc.^2.*(1./l.^2+1/Lc.^2).*(1./l.^2+1/Ls2.^2));   %ˮƽ����಻ƽ˳
figure;
loglog(l,l.^2.*S4,'r');   %���Ʋ����빦�����ܶȵ�˫��������ͼ
title('���������ߵͲ�ƽ˳','fontsize',25);   %�������������ֺ�Ϊ25
ylabel('�������ܶ�/[mm2/(rad/m)]','fontsize',20);   %y������
grid on   %������
xlabel('����/m','fontsize',20);   %x������
figure
loglog(l,l.^2.*S5,'r');
title('������������ƽ˳','fontsize',25);
ylabel('�������ܶ�/[mm2/(rad/m)]','fontsize',20);
grid on
xlabel('����/m','fontsize',20);
figure
loglog(l,l.^2.*S6,'r');
title('��������ˮƽ����಻ƽ˳','fontsize',25);
ylabel('�������ܶ�/[mm2/(rad/m)]','fontsize',20);
grid on
xlabel('����/m','fontsize',20);
