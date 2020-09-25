

%%
fmctrl_data = textread([filepath 'fmctrl_data_1337.txt']);

%%
% ¸ßµÍ¡¢¹ìÏò¡¢¹ì¾à¡¢³¬¸ß
wave_out = textread([filepath 'wave_out_1337.txt']);

N = 10000;
fmctrl_data = fmctrl_data(1:N,:);
wave_out = wave_out(1:N,:);
%%
% camco,result,amcol,xtemp,marm 
aln = textread([filepath 'aln.txt']);
aln = aln(1:N,:);


