% ����ƽ�����鷨��SPA�㷨����Ƶ����Ӧ�� ��ͬAֵ����Ӧ
% ƥ������ ������ƽ�����鷨�ı������ź�ȥ���������������������
clear all;
clc;
N = 50; % R(N) D(N-2,N)
D = zeros(48,50);
for i =1:48
    D(i,i)=1;
    D(i,i+1) = -2;
    D(i,i+2) = 1;
end

I = eye(50);

DTD = D'*D;
a = 10;
L_mid = a^2 * DTD + I;

L_mid_1 = inv(L_mid);

L = I - L_mid_1;

%  DFT :����Ҷ�任 

n = [0:1:N-1];
k = [0:1:N-1];
WN =exp(-j*2*pi/N);
nk = n'*k;
WNnk  = WN .^nk;
Xk = L * WNnk;

L_dft = [];
for i =[1:50]
   l=[];
   l = L(i,:) * WNnk;
   L_dft = [L_dft;l];
end

% ����Ƶ��������
L_dft_change = L_dft';
figure();title('L��Ƶ����Ӧͼ');
surf([1:25],[0:0.02:0.49],abs(L_dft_change(1:25,1:25)));
axis([1 25 0 0.5 0 1.5]);
ylabel('Ƶ�ʣ�Hz��');xlabel('Nȡֵ');zlabel('����');
title('L��Ƶ����Ӧͼ');

aa = [1 5 10 50 200 1000 10000];
figure ();
x =1;y=0;z=0.5;
for ai =[1:7]
    L_mid = aa(ai)^2 * DTD + I;
    L_mid_1 = inv(L_mid);
    L = I - L_mid_1;
    Xk = L * WNnk;
    plot([0:0.02:0.49],abs(Xk(25,1:25)),'color',[x y z]);hold on;
    x = x-0.16;
    y = y+0.13;
    z = z + 0.3*(-1)^(ai);
    if ai==4
        x =0;y = 1;z=0;
    elseif ai ==5
        x =0;y =0;z =1;
    elseif ai ==6
        x =1;y=0;z=0;
    end
end
legend('a =1 ','a = 5','a =10','a = 50 ' ,'a = 200' ,' a = 1000', 'a = 10000');
title('��ͬaֵ��Ӧ��Ƶ����Ӧ ');
    