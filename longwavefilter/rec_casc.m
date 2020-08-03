function [ensco1, tridb, maxV, minV, BW, slope] = rec_casc(lam,C)
%分析各级联矩形窗的频谱特性
N=128^2;
w=linspace(0,pi,N);
if C == 4
    %四矩形窗级联
    RW1 = zeros(1, N);
    RW2 = zeros(1, N);
    RW3 = zeros(1, N);
    RW4 = zeros(1, N);
    CRF = zeros(1,N);
    
    alphCR = 4/1.969;
    pCR1 = 0.895;
    pCR2 = 0.746;
    pCR3 = 0.617;
    MCR = round(lam*alphCR);
    NCR = round(lam*alphCR*pCR1);
    KCR = round(lam*alphCR*pCR2);
    LCR = round(lam*alphCR*pCR3);
    for i=1:N
        RW1(i) = ((sin(w(i)*MCR/2)/sin(w(i)/2))/MCR);
        RW2(i) = ((sin(w(i)*NCR/2)/sin(w(i)/2))/NCR);
        RW3(i) = ((sin(w(i)*KCR/2)/sin(w(i)/2))/KCR);
        RW4(i) = ((sin(w(i)*LCR/2)/sin(w(i)/2))/LCR);
    end
    CRF = RW1.*RW2.*RW3.*RW4;
    Acc = 1 - CRF;    
end
if C==3
     %三矩形窗级联
    RW1 = zeros(1, N);
    RW2 = zeros(1, N);
    RW3 = zeros(1, N);
    
    CRF = zeros(1,N);
    
    alphCR = 4/1.74;
    pCR1 = 0.829;
    pCR2 = 0.646;    
    MCR = round(lam*alphCR);
    NCR = round(lam*alphCR*pCR1);
    KCR = round(lam*alphCR*pCR2);
   
    for i=1:N
        RW1(i) = ((sin(w(i)*MCR/2)/sin(w(i)/2))/MCR);
        RW2(i) = ((sin(w(i)*NCR/2)/sin(w(i)/2))/NCR);
        RW3(i) = ((sin(w(i)*KCR/2)/sin(w(i)/2))/KCR);
        
    end
    CRF = RW1.*RW2.*RW3;
    Acc = 1 - CRF;   
end

if C == 2
     %三矩形窗级联
    RW1 = zeros(1, N);
    RW2 = zeros(1, N);    
    
    CRF = zeros(1,N);
    
    alphCR = 4/1.527;
    pCR1 = 0.718;
        
    MCR = round(lam*alphCR);
    NCR = round(lam*alphCR*pCR1);
   
   
    for i=1:N
        RW1(i) = ((sin(w(i)*MCR/2)/sin(w(i)/2))/MCR);
        RW2(i) = ((sin(w(i)*NCR/2)/sin(w(i)/2))/NCR);
        
        
    end
    CRF = RW1.*RW2;
    Acc = 1 - CRF;     
end

if C == 1
     %三矩形窗级联
    RW1 = zeros(1, N);    
    CRF = zeros(1,N);
    
    alphCR = 4/1.32;
%     pCR1 = 0.718;
        
    MCR = round(lam*alphCR);  
    for i=1:N
        RW1(i) = ((sin(w(i)*MCR/2)/sin(w(i)/2))/MCR);       
    end
    CRF = RW1;
    Acc = 1 - CRF;     
end

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