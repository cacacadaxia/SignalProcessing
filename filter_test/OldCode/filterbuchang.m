%filter �����˲����Ŀռ�Ƶ������
% ��������֤��С����
%{
clear all;
w = (10^5)/(2^14);
bs = [0 0 w^2];
as= [1 w w^2];
s=tf('s');
Q = 2^14;
Q = 100000/Q;
G= (Q^2)/(s^2 + Q*s+ Q^2); 
bode(G);
[h f]= freqs(bs,as);
figure();suptitle('����ģ�⿹���')
semilogx(f/(2*pi),20*log10(abs(h)));

%-------------------Word ����� + �����˲���----------------------
%B��z�� = (W1*T)/(1+W1*T - 1/Z)
%  һ�� �����

W1 = (10^5)/(2^17);
T = 0.002;
b = [ W1*T 0];
a = [ 1+W1*T (-1)];
[h1 f1] = freqz(b,a,800000,500);
v = 150/36;v1 = 200/9;v2 = 3500/36;
figure(101);suptitle('һ�����ֿ����')
semilogx(f1,20*log10(abs(h1)));
%semilogx(f1/(v),20*log10(abs(h1)));hold on
%semilogx(f1/(v1),20*log10(abs(h1)),'g');hold on;
%semilogx(f1/(v2),20*log10(abs(h1)),'r');

figure();suptitle('һ�����ֿ����/����');
T = 0.25/v;
b = [ W1*T 0];
a = [ 1+W1*T (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v./f1,20*log10(abs(h1)));hold on
T = 0.01125;
b = [ W1*T 0];
a = [ 1+W1*T (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v1./f1,20*log10(abs(h1)),'g');hold on;
T = 0.25/v2;
b = [ W1*T 0];
a = [ 1+W1*T (-1)];
[h1 f1] = freqz(b,a,800000,500);
semilogx(v2./f1,20*log10(abs(h1)),'r');hold on
% һ�� ����
%bc = [1+ (W1*T/2) -1+(W1*T/2)];
%ac = [W1*T 0];
%bc = [-1+(W1*T/2) 1+ (W1*T/2) ];
%ac = [0 W1*T];
%bc = [-1+(W1*T/2) 1+ (W1*T/2) ];
bc = [ 1+(W1*0.045/2)  (-1)+(W1*0.045/2)  ];
ac = [1 0];
[h1c f1c] = freqz(bc,ac,800000,500);
%figure();suptitle('һ�����ֲ���');
bc = [ 1+(W1*(0.25/v)/2)  (-1)+(W1*(0.25/v)/2)  ];
ac = [W1*(0.25/v) 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v./f1c,20*log10(abs(h1c)));hold on
bc = [ 1+(W1*0.01125/2)  (-1)+(W1*0.01125/2)  ];
ac = [W1*0.01125 0];
[h1c f1c] = freqz(bc,ac,800000,500);
semilogx(v1./f1c,20*log10(abs(h1c)),'g');hold on;
bc = [ 1+(W1*(0.25/v2)/2)  (-1)+(W1*(0.25/v2)/2)  ];
ac = [W1*0.25/v2 0];
[h1c f1c] = freqz(bc,ac,800000,500);hold on;
semilogx(v2./f1c,20*log10(abs(h1c)),'r');

% -----һ�����ֿ���� + ����-------------
%B1 = [1+W1*T/2 W1*T/2-1];
%A1 = [1+W1*T -1];
%B1 = [ 2*W1*T+(W1*T)^2 (W1*T)^2-2*W1*T];
%A1 = [2+2*W1*T (-2)]
T = (0.25/v);%0.045
B1 = [2*W1*T+(W1^2)*T*(0.25/v) (W1^2)*T*(0.25/v)-2*W1*T];
A1 = [2*W1*(0.25/v)+2*(W1^2)*T*(0.25/v) (-2)*W1*(0.25/v)];
[H1 F1] = freqz(B1,A1,800000,500);
%figure();suptitle('һ�����ֿ����+����');
semilogx(v./F1,20*log10(abs(H1)));hold on
T = 0.01125;% 0.0225
B1 = [2*W1*T+(W1^2)*T*0.01125 (W1^2)*T*0.01125-2*W1*T];
A1 = [2*W1*0.01125+2*(W1^2)*T*0.01125 (-2)*W1*0.01125];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v1./F1,20*log10(abs(H1)),'g');hold on;
T = 0.25/v2 ;%0.0075
B1 = [2*W1*T+(W1^2)*T*(0.25/v2) (W1^2)*T*(0.25/v2)-2*W1*T];
A1 = [2*W1*(0.25/v2)+2*(W1^2)*T*(0.25/v2) (-2)*W1*(0.25/v2)];
[H1 F1] = freqz(B1,A1,800000,500);
semilogx(v2./F1,20*log10(abs(H1)),'r');hold on;
%axis([ 0.01 10^3 -5 5 ]);
%}



% -----���׿���� + ����-------------
% ���׿����
% w2 = (10^5)/(2^14);
% t = 9/150;
% w2t = w2*t;
% v = 150/36;v1 = 200/9;v2 = 3500/36;
% %b2 = [0 0 (w2t)^2];
% %a2 = [1 -(2+w2t) 1+w2t+(w2t)^2];
% b2 = [(w2t)^2 0 0];
% a2 = [1+w2t+(w2t)^2 -(2+w2t) 1];
% [h2 f2] = freqz(b2,a2,800000,500);
% figure();suptitle ('�������ֿ�����˲���');
% semilogx(v./f2,20*log10(abs(h2)));hold on;
% t = 0.01125;
% w2t = w2*t;
% b2 = [(w2t)^2 0 0];
% a2 = [1+w2t+(w2t)^2 -(2+w2t) 1];
% [h2 f2] = freqz(b2,a2,800000,500);
% semilogx(v1./f2,20*log10(abs(h2)),'g');hold on;
% t = 9/3500;
% w2t = w2*t;
% b2 = [(w2t)^2 0 0];
% a2 = [1+w2t+(w2t)^2 -(2+w2t) 1];
% [h2 f2] = freqz(b2,a2,800000,500);
% semilogx(v2./f2,20*log10(abs(h2)),'r');hold on;


% ���ײ����˲���
% %a=��_2 l  b=��_2 ��(1-l^2 ) l=1/2 ��_2=10^5/2^14 a=��_2 l  b=��_2 ��(1-l^2 )
% %l=1/2 ��_2=10^5/2^14
% w2 = (10^5)/(2^14);
% t = 9/150;
% w2t = w2*t;
% a = w2/2; aa = a^2;
% bb = 3*(w2^2)/4;
% %b2 = [0 0 (w2t)^2];
% %a2 = [1 -(2+w2t) 1+w2t+(w2t)^2];
% b2 = [1+a*t+aa*t*t -2*(1-aa*t*t) 1-a*t+aa*t*t];
% a2 = [(aa+bb)*t*t 0 0];
% [h2 f2] = freqz(b2,a2,1000000,500);
% %figure();suptitle ('�������ֲ����˲���');
% semilogx(v./f2,20*log10(abs(h2)));hold on
% t = 0.01125;
% b2 = [1+a*t+aa*t*t -2*(1-aa*t*t) 1-a*t+aa*t*t];
% a2 = [(aa+bb)*t*t 0 0];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v1./f2,20*log10(abs(h2)),'g');hold on;
% t = 9/3500;
% b2 = [1+a*t+aa*t*t -2*(1-aa*t*t) 1-a*t+aa*t*t];
% a2 = [(aa+bb)*t*t 0 0];
% [h2 f2] = freqz(b2,a2,1000000,500);
% semilogx(v2./f2,20*log10(abs(h2)),'r');hold on;
% %semilogx(f2/v2,angle(h2),'r');


% % ���׿����+�����˲���
%a=��_2 l  b=��_2 ��(1-l^2 ) l=1/2 ��_2=10^5/2^14 a=��_2 l  b=��_2 ��(1-l^2 )
%l=1/2 ��_2=10^5/2^14
w2 = (10^5)/(2^14);
t = 0.002;
w2t = w2*t;
a = w2/2; aa = a^2;
bb = 3*(w2^2)/4;
%b2 = [0 0 (w2t)^2];
%a2 = [1 -(2+w2t) 1+w2t+(w2t)^2];
t = 9/150;w2t = w2*t;
b2 = [(w2t^2)*(1+a*t+aa*t*t) -2*(1-aa*t*t)*(w2t^2) (1-a*t+aa*t*t)*(w2t^2)];
a2 = [(1+w2t+(w2t)^2)*(aa+bb)*t*t -(2+w2t)*(aa+bb)*t*t (aa+bb)*t*t];
[h2 f2] = freqz(b2,a2,1000000,500);
%figure();suptitle ('�������ֿ���������˲���');
semilogx(v./f2,20*log10(abs(h2)));hold on
%t = 0.0225;
t= 0.01125;w2t = w2*t;
%b2 = [1+a*t+aa*t*t -2*(1-aa*t*t) 1-a*t+aa*t*t];
b2 =[(w2t^2)*(1+a*t+aa*t*t) -2*(1-aa*t*t)*(w2t^2) (1-a*t+aa*t*t)*(w2t^2)];
a2 = [(1+w2t+(w2t)^2)*(aa+bb)*t*t -(2+w2t)*(aa+bb)*t*t (aa+bb)*t*t];
[h2 f2] = freqz(b2,a2,1000000,500);
semilogx(v./f2,20*log10(abs(h2)),'g');hold on;
t = 9/3500;w2t = w2*t;
b2 = [(w2t^2)*(1+a*t+aa*t*t) -2*(1-aa*t*t)*(w2t^2) (1-a*t+aa*t*t)*(w2t^2)];
a2 = [(1+w2t+(w2t)^2)*(aa+bb)*t*t -(2+w2t)*(aa+bb)*t*t (aa+bb)*t*t];
[h2 f2] = freqz(b2,a2,1000000,500);
semilogx(v2./f2,20*log10(abs(h2)),'r');hold on
%semilogx(f2/v2,angle(h2),'r');
%---------------------------------------------------------

%}




















%--------------------pdf�����ƽ˳���---------------------------
% �����˲���
%{
clc;
s = tf('s');
a = (10^5)/(2^14);
s1 = -3.0518;
s2 = -0.8019;

G = (a^2)/(s^2 +a*s +a^2);
bode (G);
%}
%{
Ls = 0.25;
%s1 = -3.0518;
%s2 = -3.0518;
s1 = -0.0819;
s2 = -0.0819;
a = Ls^2 *s1 *s2;

A = [a 0 0];
%A = [0 0 a];
v = 50/9;
b1 = v^2;
b2 = 2* v^2 +Ls*v*(s1+s2);
b3 = v^2 + Ls*v*(s1+s2);
%B = [b1 b2 b3];
B = [b3 b2 b1];
[h f] = freqz(A,B);
semilogx(f,20*log10(abs(h)),'r');hold on

%A = [0 0 a];
v1 = 100/9;
b1 = v1^2;
b2 = 2* v1^2 +Ls*v1*(s1+s2);
b3 = v1^2 + Ls*v1*(s1+s2);
%B = [b1 b2 b3];
B = [b3 b2 b1];
[h1 f1] = freqz(B,A);
semilogx(f1,20*log10(abs(h1)),'g');hold on

%A = [0 0 a];
v2 = 100/3;
b1 = v2^2;
b2 = 2* v2^2 +Ls*v2*(s1+s2);
b3 = v2^2 + Ls*v2*(s1+s2);
%B = [b1 b2 b3];
B = [b3 b2 b1];
[h2 f2] = freqz(B,A);
semilogx(f2,20*log10(abs(h2)));hold on
%}


%�����ģ���˲���
%{
a = (145000)/(2^14);
v = 50/9;
v1 = 100/9;
v2 = 100/3;
A = [1 a/(2*pi*v) (a/(2*pi*v))^2];
B = [0 0 (a/(2*pi*v))^2];
[h,y] = freqs(B,A);
figure();
semilogx(y,20*log10(abs(h)));hold on
A1 = [1 a/(2*pi*v1) (a/(2*pi*v1))^2 ];
B1 = [0 0 (a/(2*pi*v1))^2];
[h1,y1] = freqs(B1,A1);semilogx(y1,20*log10(abs(h1)),'g');hold on;
A2 = [1 a/(2*pi*v2) (a/(2*pi*v2))^2 ];
B2 = [0 0 (a/(2*pi*v2))^2];
[h2,y2] = freqs(B2,A2);semilogx(y2,20*log10(abs(h2)),'r');


a = (235000)/(2^14);
v = 50/9;
v1 = 100/9;
v2 = 100/3;
A = [(a)^2 a 1];
B = [0 0 a^2];
%A = [1 a (a)^2];
%B = [0 0 (a)^2];
[h,y] = freqs(B,A);
figure();
semilogx(y/(2*pi*v),20*log10(abs(h)));hold on
semilogx(y/(2*pi*v1),20*log10(abs(h)),'g');hold on;
semilogx(y/(2*pi*v2),20*log10(abs(h)),'r');
%}


% �󼫵� 
%{
a = (10^5)/(2^14);
num = [0 0 a^2];
den = [1 a a^2];

num = [a^2 0 0];
den = [a^2 a 1];
sys = TF(num,den);
PZMAP(sys);
[p,z] = PZMAP(sys)
%}
%{
num = [1 5 6];
den = [1 2 2 1];
sys = TF(num,den);
PZMAP(sys);
[p,z] = PZMAP(sys)
%}