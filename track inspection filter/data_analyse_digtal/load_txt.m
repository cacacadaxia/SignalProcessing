

%%
fmctrl_data = textread([filepath 'fmctrl_data_1337.txt']);

%%
% ¸ßµÍ¡¢¹ìÏò¡¢¹ì¾à¡¢³¬¸ß
wave_out = textread([filepath 'wave_out_1337.txt']);
%%
% N = 20000;
fmctrl_data = fmctrl_data(start_pos:start_pos+N-1,:);
wave_out = wave_out(start_pos:start_pos+N-1,:);
%%
% camco,result,amcol,xtemp,marm 
aln = textread([filepath 'aln.txt']);
aln = aln(start_pos:start_pos+N-1,:);


