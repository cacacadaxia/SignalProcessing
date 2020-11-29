% =========================================================================
%
%                  测试∫∫(pitch*0.25)对波长幅值的变化到底有多少？
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 9月16日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.目前传递函数与画图得到的结果对不上，需要重新对比
%        2.这边是没毛病
%       3. 
%--------------------------------------------------------------------------


clear all;
close all;
del_x = 0.25;

WaveLen = 10;
WaveMag = 10;
fspace = @(x)(1*(WaveMag*1e-3*sin(2*pi/WaveLen *x) ));
x = 0:del_x:160;
N = length(x);
longwave = fspace(x);


%%
L = 3;
z1z2_2 = (fspace(x) + fspace(x-L))./2;
figure;plot(z1z2_2*1e3);
cau(z1z2_2*1e3)/2/WaveMag

%%
% pitch = (fspace(x) - fspace(x-L))./L;
pitch = (fspace(x+L) - fspace(x))./L;

z = 0;
for i = 1:length(pitch)
    z(i+1) = z(i) + pitch(i)*0.25;
end

figure;plot(z);
cau(z*1e3)/2/WaveMag

%%
function out = cau(in)

out = max(in)-min(in);

end






