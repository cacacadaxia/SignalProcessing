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

Trap1 = zeros(1,N);
Trap2 = zeros(1,N);
Trap3 = zeros(1,N);
Trap = zeros(1,N);

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

% lamda = 2*pi./w./4;
% cutn = 3;
% ensco = ensco1([cutn:end]);
% lamda1 = lamda([cutn:end]);
% maxV = max(ensco);
% maxL = 0;
% tridb = 0;
% minL = 0;
% tol = 0.02;
% index = abs(ensco-0.707) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));
% % indexM = indexNumeric;
% if ~isnan(indexM) 
%     tridb = (lamda1(indexM));
% end
% 
% dx=diff(lamda1);
% dy=diff(ensco);
% dydx = dy./dx;
% slope = dydx(indexM);
% 
% index = abs(ensco-1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));
% % indexM = indexNumeric;
% if ~isnan(indexM) 
%     maxL = (lamda1(indexM));
% end
% index = abs(ensco-0.1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));
% % indexM = indexNumeric;
% if ~isnan(indexM) 
%     minL = (lamda1(indexM));
% end
% BW = minL - maxL;
% minV = min(ensco);

