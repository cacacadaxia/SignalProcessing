function [w]=indirect()   %������
l=1;   %��������
v=150/3.6;   %����
fs=2*v/l;   %����Ƶ��
f=fopen('0412.txt','r');   %��04��12�¹���ߵͲ�ƽ˳ʵ������
y=fscanf(f,'%e');
%��ȡԭʼ����
nfft=2^15;    %FFT�ĵ���
R=xcorr(y,'unbiased');     %����غ���
F=fft(R,nfft);     %���ٸ���Ҷ�任
p=abs(F);     %��⹦�����ܶ�
index=0:round(nfft/2-1);     %ʱ������ȡһ�� 
k=index*fs/nfft;     %��Ӧ��Ƶ��
figure      %�´���
w=loglog(k,p(index+1));     %��Ƶ�ʺ͹������ܶȵ�˫��������ͼ
title('����ط�','fontsize',25)    %������,�ֺ�Ϊ25
xlabel('Ƶ��(Hz)','fontsize',20);     %x������
grid on      %������
ylabel('�������ܶ�(mm2)','fontsize',20);    %y������

