function [triDB, maxV,minV, BW, slope] = filter_para(F)
%triDB: -3db
%maxV: 通带最大波动
%BW：带宽
%slope：3分贝点斜率
ensco1 = F;
N=128^2;
w=linspace(0,pi,N);
lamda = 2*pi./w./4;
cutn = 3;
ensco = ensco1([cutn:end]);
lamda1 = lamda([cutn:end]);
[maxV1,index] = max(ensco);

maxL = 0;
triDB = 0;
minL = 0;
tol = 0.02;
index = abs(ensco-0.707) < tol;
indexNumeric = find(index); 
indexM = floor(median(indexNumeric));

if ~isnan(indexM) 
    triDB = (lamda1(indexM));
    dx=diff(lamda1);
    dy=diff(ensco);
    dydx = dy./dx;
    slope = dydx(indexM);
else
    slope = NaN;
end

% index = abs(ensco-maxV1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));

index = abs(ensco - 1.0) < tol;
indexM = find(index,1,'first');
if ~isnan(indexM) 
    maxL = (lamda1(indexM));
    
end
tmp1 = indexM;

ensco2 = ensco([indexM:end]);
maxV = max(abs(ensco2-1));

% index = abs(ensco-0.1) < tol;
% indexNumeric = find(index); 
% indexM = floor(median(indexNumeric));

index = abs(ensco - 0.05) < tol;
indexM = find(index,1,'last');
if ~isnan(indexM) 
    minL = (lamda1(indexM));
end
ensco3 = ensco([1:indexM]);
tmp2 = indexM;
minV = max(abs(ensco3-0));


PLOT = 0;
if PLOT
    fprintf('Results of transition zone:\n')
    disp([lamda1(tmp1) , lamda1(tmp2)]);
    disp([(tmp1) , (tmp2)]);
end

BW = minL - maxL;
