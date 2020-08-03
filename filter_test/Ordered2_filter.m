% ----- ���׿���� + ���� -------------
clear all;
clc;

% ���׵�ͨ�˲���
w2 = (10^5)/(2^14);
v1=  15/3.6;
v2 = 80/3.6;
v3 = 350/3.6;     % ��Щ�ٶ��Ǵ��������ģ�

t1= 0.25/v1;
t2= 0.25/v2;
t3= 0.25/v3;    % ��ͬ���ٶȴ�����ʲô��
w2t1 = w2*t1;
b2 = [(w2t1)^2 0 0];
a2 = [1+w2t1+(w2t1)^2  ,- (2 + w2t1)  ,1];
[h2 f2] = freqz(b2,a2,800000,500);
figure;suptitle ('�������ֿ�����˲����Ͳ����˲���');
semilogx(v1./f2,20*log10(abs(h2)));hold on;
% freqz(b2,a2,800000,500);
w2t2 = w2*t2;
b2 = [(w2t2)^2 0 0];
a2 = [1+w2t2+(w2t2)^2, -(2+w2t2) , 1];
[h2 f2] = freqz(b2,a2,800000,500);
semilogx(v2./f2,20*log10(abs(h2)),'g');hold on;

w2t3 = w2*t3;
b2 = [(w2t3)^2 0 0];
a2 = [1+w2t3+(w2t3)^2 ,-(2+w2t3) 1];
[h2 f2] = freqz(b2,a2,800000,500);
semilogx(v3./f2,20*log10(abs(h2)),'r');hold on;


% ���ײ����˲���

% w2 = (10^5)/(2^14);
% a = w2/2; aa = a^2;
% bb = 3*(w2^2)/4;
% b2 = [1+a*t1+aa*t1*t1 ,-2*(1-aa*t1*t1), 1-a*t1+aa*t1*t1];
% a2 = [(aa+bb)*t1*t1 0 0];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v1./f2,20*log10(abs(h2)));hold on          % Ϊʲô�ٶ�Ҫ�������
% 
% b2 = [1+a*t2+aa*t2*t2 -2*(1-aa*t2*t2) 1-a*t2+aa*t2*t2];
% a2 = [(aa+bb)*t2*t2 0 0];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v2./f2,20*log10(abs(h2)),'g');hold on;
% 
% b2 = [1+a*t3+aa*t3*t3 -2*(1-aa*t3*t3) 1-a*t3+aa*t3*t3];
% a2 = [(aa+bb)*t3*t3 0 0];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v3./f2,20*log10(abs(h2)),'r');hold on;


% % ���׿����+�����˲���
% %a=��_2 l  b=��_2 ��(1-l^2 ) l=1/2 ��_2=10^5/2^14 a=��_2 l  b=��_2 ��(1-l^2 )
% %l=1/2 ��_2=10^5/2^14
% w2 = (10^5)/(2^14);
% a = w2/2; aa = a^2;
% bb = 3*(w2^2)/4;
% 
% b2 = [(w2t1^2)*(1+a*t1+aa*t1*t1) -2*(1-aa*t1*t1)*(w2t1^2) (1-a*t1+aa*t1*t1)*(w2t1^2)];
% a2 = [(1+w2t1+(w2t1)^2)*(aa+bb)*t1*t1 -(2+w2t1)*(aa+bb)*t1*t1 (aa+bb)*t1*t1];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v1./f2,20*log10(abs(h2)));hold on
% 
% b2 =[(w2t2^2)*(1+a*t2+aa*t2*t2) -2*(1-aa*t2*t2)*(w2t2^2) (1-a*t2+aa*t2*t2)*(w2t2^2)];
% a2 = [(1+w2t2+(w2t2)^2)*(aa+bb)*t2*t2 -(2+w2t2)*(aa+bb)*t2*t2 (aa+bb)*t2*t2];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v2./f2,20*log10(abs(h2)),'g');hold on;
% 
% b2 = [(w2t3^2)*(1+a*t3+aa*t3*t3) -2*(1-aa*t3*t3)*(w2t3^2) (1-a*t3+aa*t3*t3)*(w2t3^2)];
% a2 = [(1+w2t3+(w2t3)^2)*(aa+bb)*t3*t3 -(2+w2t3)*(aa+bb)*t3*t3 (aa+bb)*t3*t3];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v3./f2,20*log10(abs(h2)),'r');hold on
% 
% xlabel('������m��');ylabel('��ֵ(dB)');
% % ---------------------------------------------------------