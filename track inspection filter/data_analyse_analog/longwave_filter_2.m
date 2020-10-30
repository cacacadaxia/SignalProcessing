
% =========================================================================
%
%                       ���������˲���
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��28��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------

function yL = longwave_filter_2( amcol , varargin )

if nargin == 5
    NTri_1 = varargin{1};
    NTri_2 = varargin{2};
    NRec_1 = varargin{3};
    NRec_2 = varargin{4};

    %% �������N��ֵ
    mainTriM = (NTri_1-1)/2;
    auxTriM = (NTri_2-1)/2;
    mainRecM = (NRec_1-1)/2;
    auxRecM = (NRec_2-1)/2;
else
    %% Ĭ�ϵ�
    mainTriM = 160;
    auxTriM = 40;
    mainRecM = 160;
    auxRecM = 280;
end

%% ��Ȼ�˲�����ʵ����ûʲô���⣬����ʱ����Ҫ����
% FSCAL = 0.1
% -1.036, 0.036, -0.25, 0.25, FSCAL * 2.0
% �����趨
FSCAL = 0.1;
% longWaveFilterBufferSize = 2048;
mainTriWinCoef = -1.036;
auxTriWinCoef = 0.036;
mainRecWinCoef = -0.25;
auxRecWinCoef = 0.25;
scaleCoef =  FSCAL * 2.0;

%%
mainTriFront = 0; auxTriFront = 0; mainTriBehind = 0; auxTriBehind = 0;%%���Ǵ�
mainRecSum = 0; auxRecSum = 0; mainTriSum = 0; auxTriSum = 0;
% MianREC = 0; AuxREC = 0; MainTRI = 0; AuxTRI = 0;
Num = 2048;     %%ѭ�����д洢����
data = zeros(Num,1);
center_pos = 1000;
dot = 0; yL = zeros(2,1);
for i = 1:length(amcol)
    %% ��������
    data_input_index = center_pos + auxRecM;
    data(find_index(data_input_index)) = amcol(i);
    %% ���δ�
    mainRecSum = mainRecSum + data(find_index(center_pos + mainRecM )) - data(find_index(center_pos - mainRecM - 1 ));
    auxRecSum = auxRecSum + data(find_index(center_pos + auxRecM)) -  data(find_index(center_pos - auxRecM - 1 ));
    %% ���Ǵ�
    mainTriFront = mainTriFront + (data(find_index(center_pos + mainTriM-1)) - data(find_index(center_pos-1))) / mainTriM;
    mainTriBehind = mainTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - mainTriM - 1 )) ) / mainTriM;
    
    auxTriFront = auxTriFront + (data(find_index(center_pos + auxTriM-1)) - data(find_index(center_pos-1))) / auxTriM;
    auxTriBehind = auxTriBehind + ( data(find_index(center_pos - 1)) - data(find_index(center_pos - auxTriM - 1 )) ) / auxTriM;
    
    mainTriSum = mainTriSum + mainTriFront - mainTriBehind;
    auxTriSum = auxTriSum + auxTriFront - auxTriBehind;
    %% �������
    yout = data(find_index(center_pos)) + mainTriWinCoef*mainTriSum/mainTriM + auxTriWinCoef*auxTriSum/auxTriM + ...
        mainRecWinCoef*mainRecSum/(2*mainRecM+1) + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
    %         yout = data(center_pos) - 1.036*mainTriSum/mainTriM + 0.036*auxTriSum/auxTriM;
    %     + auxRecWinCoef*auxRecSum/(2*auxRecM+1);
    
    %% 
    
    yout_save(i) = yout;
    dot = dot - yout * 2 * auxTriM / ( auxTriM * 2 + 1 );
    yL(i+1) = yL(i) + dot;
    %% ����
    center_pos = center_pos + 1;
    center_pos = find_index(center_pos);
end

yL = yL * scaleCoef;
for i = 1:length(yL)
    yL(i) = point3filter(yL(i));
end

end


function out = find_index(in)
%% �����˲����Ķ�����
Num = 2048;
if mod(in,Num) == 0
    out = Num;
else
    out = mod(in,Num);
end
end

function out = point3filter(in)
%% �򵥵������˲�
persistent x;
if isempty(x)
    x = zeros(3,1);
end
x(3) = in;
out = sum(x)/3;
%%
x(1) = x(2);
x(2) = x(3);

end

