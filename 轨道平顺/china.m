function [q]=china()     %��������
t=fopen('BJ.txt','r');    %���й�������߲����ļ�
p=fscanf(t,'%e');    %��ȡ����
fs=1000;    %���ò���Ƶ��
l=1:1/fs:30;    %������Χ
A=p(1,1);B=p(2,1);C=p(3,1);D=p(4,1);E=p(5,1);F=p(6,1);G=p(7,1);   %����������ֵ
S=A.*(l.^2+B*l.^3+C*l.^4)./(l.^2.*(G*l.^4+F*l.^3+E*l.^2+D*l+1));               %����λת����Ĺ���׹�ʽ
q=loglog(l,l.^2.*S);     %�������빦�����ܶȵ�˫��������ͼ
title('����������ߵ�','fontsize',25)     %�������������ֺ�Ϊ25
xlabel('����/m','fontsize',20)     %x������
grid on    %������
ylabel('�������ܶ�/[mm2/(rad/m)]','fontsize',20)    %y��������
figure
A=p(8,1);B=p(9,1);C=p(10,1);D=p(11,1);E=p(12,1);F=p(13,1);G=p(14,1);           %����������ֵ
S=A.*(l.^2+B*l.^3+C*l.^4)./(l.^2.*(G*l.^4+F*l.^3+E*l.^2+D*l+1));
loglog(l,l.^2.*S);
title('�������������','fontsize',25)
xlabel('����/m','fontsize',20)
grid on
ylabel('�������ܶ�/[mm2/(rad/m)]','fontsize',20)
figure
A=p(15,1);B=p(16,1);C=p(17,1);D=p(18,1);E=p(19,1);F=p(20,1);G=p(21,1);         %����������ֵ
S=A.*(l.^2+B*l.^3+C*l.^4)./(l.^2.*(G*l.^4+F*l.^3+E*l.^2+D*l+1));
loglog(l,l.^2.*S);
title('�������ˮƽ','fontsize',25)
xlabel('����/m','fontsize',20)
grid on
ylabel('�������ܶ�/[mm2/(rad/m)]','fontsize',20)
