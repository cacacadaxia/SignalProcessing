fmctrl_data = textread('fmctrl_data_1337.txt');

% ¸ßµÍ¡¢¹ìÏò¡¢¹ì¾à¡¢³¬¸ß
wave_out = textread('wave_out_1337.txt');


% sita_b = textread('tmp_zhongjian_1337.txt');

N = 10000;
fmctrl_data = fmctrl_data(1:N,:);
wave_out = wave_out(1:N,:);
% sita_b = sita_b(1:N,:);

aln = textread('aln.txt');
if length(aln)>N
    aln = aln(1:N,:);
end

