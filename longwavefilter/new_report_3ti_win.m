


clear all;
close all;
N=128^2;
w=linspace(0,pi,N);
lamda = 2*pi./w./4;
lam = [];

[Acc_Trap_42, TriDB_Trap_42, maxV_Trap_42, minV_Trap_42, BW_Trap_42, slope_Trap_42] = tri_paral(42);

fig1=figure;semilogx((lamda), (abs(Acc_Trap_42)), '-b','LineWidth',1);
plot_line_func(abs(Acc_Trap_42),lamda,42);
magnify(fig1)
xlabel('Wavelength /m')
ylabel('Magnitude /dB')



function [ensco1, tridb, maxV, minV, BW, slope] = tri_paral(lam)
%
%分析3并联梯形窗各波长滤波幅频特性,直接绘制法
N=128^2;


% 这个参数是从哪里来的？
w1 = 1.0;
w2 = 0.75;
w3 = -0.75;

alph = 2.56;

M = round(alph*lam);
w = linspace(0,pi,N);

p1 = 1;
p2 = 1.75;
p3 = 2.60;
%2级联矩形窗
%二矩形窗参数__no. 1 Trapezoid
a1 = 0.718;

Mt1=round(p1*M);%143;%65;%256;%100;%64;%256;
Mt2=round(Mt1*a1+0.5);%127;%59;%229;%90;%58;%229;

%二矩形窗参数__no. 2 Trapezoid
a2 = 0.4;

Kt1=round(p2*M);%143;%65;%256;%100;%64;%256;
Kt2=round(Kt1*a2+0.5);%127;%59;%229;%90;%58;%229;

%二矩形窗参数__no. 3 Trapezoid
a3 = 0.475;

St1=round(p3*M);%143;%65;%256;%100;%64;%256;
St2=round(St1*a3+0.5);%127;%59;%229;%90;%58;%229;


%% report from yu(modify)
% w1 = 0.95;
% w2 = 0.6;
% w3 = -0.55;
% Mt1 = 109;
% Mt2 = 79;
% Kt1 = 157;
% Kt2 = 111;
% St1 = 221;
% St2 = 145;
% 
% Trap1 = zeros(1,N);
% Trap2 = zeros(1,N);
% Trap3 = zeros(1,N);
% Trap  = zeros(1,N);


%%
for i=1:N
    Trap1(i)= ((sin(w(i)*Mt1/2)/sin(w(i)/2))/Mt1)*((sin(w(i)*Mt2/2)/sin(w(i)/2))/Mt2);
    Trap2(i)= ((sin(w(i)*Kt1/2)/sin(w(i)/2))/Kt1)*((sin(w(i)*Kt2/2)/sin(w(i)/2))/Kt2);
    Trap3(i)= ((sin(w(i)*St1/2)/sin(w(i)/2))/St1)*((sin(w(i)*St2/2)/sin(w(i)/2))/St2);
end


Trap = w1.*Trap1+w2.*Trap2+w3.*Trap3;
Acc = 1 - Trap;
% Acc1 = 1- Trap1;
% Acc2 = 1- Trap2;
% Acc3 = 1- Trap3;
ensco1 = Acc;
[tridb, maxV, minV, BW, slope] = filter_para(Acc);
end

function plot_line_func(mag,lamda,waveLen)

% mag = abs(Acc_Trap_42);

%% find index
index = find_707(mag);
ind = find(index==1);
horiz = mag(ind(1));
horiz_x = lamda(ind(1));
verti = waveLen;

%%  plot
% figure;
semilogx((lamda), ((mag)),'LineWidth',1);
hold on;
semilogx([1,1000] , [horiz,horiz] ,'--','LineWidth',1)
semilogx( [horiz_x,horiz_x],[0.3,0.9],'--','LineWidth',1)
semilogx( [verti,verti],[0.3,0.9],'--','LineWidth',1)
end

%% function
function index = find_707(mag)
tol = 0.00;
index = zeros(1,length(mag));
while(sum(index)==0)
    index1 = mag <0.7079+tol;
    index2 = mag >0.7079-tol;
    index = index1.*index2;
    tol = tol+0.001;
end
end







