




function [ensco1, tridb, maxV, minV, BW, slope] = tri_rec_paral_my(lam, para)
%三角窗+矩形窗并联幅频特性分析
N=128^2;
w=linspace(0,pi,N);

Tri1 = zeros(1, N);
Tri2 = zeros(1,N);
Rec1 = zeros(1,N);
Rec2 = zeros(1,N);

switch (lam)
    case 25
        M1 = para(1,1);
        M2 = para(1,2);
        M3 = para(1,3);
        M4 = para(1,4);
    case 42
        M1 = para(2,1);
        M2 = para(2,2);
        M3 = para(2,3);
        M4 = para(2,4);
    case 70
        M1 = para(3,1);
        M2 = para(3,2);
        M3 = para(3,3);
        M4 = para(3,4);
    case 120
        M1 = para(4,1);
        M2 = para(4,2);
        M3 = para(4,3);
        M4 = para(4,4);
end
        

% m1 = floor(lam*2.4+1.9+0.5);
% M1 = m1 + 1 - mod(m1,2);
% m2 = floor(M1*0.3216 + 0.5);
% M2 = m2 + 1 - mod(m2,2);
% M3 = round(M1*1.1228)*2+1;
% M4 = round(M1*1.5439)*2+1;

% 矩形窗的固有形式

for i=1:N
%     三角窗
    Tri1(i)= ((sin(w(i)*M1/2)/sin(w(i)/2))/M1)*((sin(w(i)*M1/2)/sin(w(i)/2))/M1);
    Tri2(i)= ((sin(w(i)*M2/2)/sin(w(i)/2))/M2)*((sin(w(i)*M2/2)/sin(w(i)/2))/M2);
%     矩形窗
    Rec1(i)= ((sin(w(i)*M3/2)/sin(w(i)/2))/M3);
    Rec2(i)= ((sin(w(i)*M4/2)/sin(w(i)/2))/M4);
end
MulFilter = 1.036*Tri1-0.036*Tri2+0.25*(Rec1-Rec2);
% a = 1.0;
% b = 1.0;
% MulFilter = Tri1+a*0.036*(Tri1-Tri2)+b*0.25*(Rec1-Rec2);
ensco1 = 1 - MulFilter;
Acc = ensco1;
[tridb, maxV, minV, BW, slope] = filter_para(Acc);

% lamda = 2*pi./w./4;
% cutn = 3;
% ensco = ensco1([cutn:end]);
% lamda1 = lamda([cutn:end]);
% [maxV,index] = max(ensco);
% 
% % segensco= abs(ensco([1:index])-1);
% % maxV = max(segensco);
% 
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
% 
% index = abs(ensco-1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));
% % indexM = indexNumeric;
% if ~isnan(indexM) 
%     maxL = (lamda1(indexM));
% end
% 
% index = abs(ensco-0.1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));
% % indexM = indexNumeric;
% if ~isnan(indexM) 
%     minL = (lamda1(indexM));
% end
% BW = minL - maxL;
% minV = min(ensco);











