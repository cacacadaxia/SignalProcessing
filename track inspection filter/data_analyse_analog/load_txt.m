


% =========================================================================
%
%                  ��ȡ����
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 10��27��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.
%        2.
%       3. 
%--------------------------------------------------------------------------
%%
tp{1} = start_pos;tp{2} = N;tp{3} = filepath;

%%
aln = read_txt_1('aln.txt',tp);
fmctrl_data = read_txt_1('fmctrl_data_1337.txt',tp);
output_wave = read_txt_1('output_wave.txt',tp);
LongWaveResultForAln_L = read_txt_1('LongWaveResultForAln_L.txt',tp);
gppro_data = read_txt_1('gppro.txt',tp);

%% function
function out = read_txt_1(name,tp)

start_pos = tp{1} ;N = tp{2};filepath = tp{3};
out = textread([filepath,name]);
if length(out)>(N + start_pos)
    out = out(start_pos:start_pos+N-1,:);
else
    out = out(start_pos:end,:);
end
end