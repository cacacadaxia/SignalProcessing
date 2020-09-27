

% =========================================================================
%
%                  生成轨道不平顺谱
%
% =========================================================================
%
%　(C)2019-2022 铁道科学研究院-基础所
%   版本：V1.0
%   日期：2020年 6月15日
%   作者：s.m.
%--------------------------------------------------------------------------
%  功能： 1.
%        2.
%       3. 
%       4. 
%     该函数以美国六级线路轨道高低不平顺功率谱密度为例，采用相对最优的方法逆傅立叶变换法，产生时间域的轨道不平顺激励
%       ，取空间波长为0.5-50m，速度100Km/h（100/3.6=27.78m/s），
%     采样时间最小步长delta为0.0001，上截止频率为V/0.5=27.78/0.5=55.56，下截止频率V/50=27.78/50=0.5566
%     若使用模块From Workspace仿真，需要双击文件USA6_ex.mat或者在matlab命令窗口输入load ex_USA6 
%     该函数被rail_excitation.mdl中的双击模块Double click to Parameterize使用，并且该模块中使用了load ex_USA6
%     参考文献《！！！！！铁路轨道不平顺功率谱分析与数值模拟毕业论文》西南交大的
%--------------------------------------------------------------------------
clc;
v=100/3.6;%列车运行速度
delta=0.0001;%时间间隔
Nr=2^18;%时域和频域采样点数
df=1/(Nr*delta);%频域采样间隔
fu=v/0.5;%上截止频率
fl=v/50;%下截止频率
Nf=fix((fu-fl)/df);%有效频率段内的采样点数
No=round(fl/df);%下截止频率以下点数
t=zeros(1,No);%下截止频率以下点数置为零
m=zeros(1,Nr/2-Nf-No);%上截止频率到Nf置为零
f=linspace(fl,fu,Nf);%将上下截止频率之间等分Nf份
fc=82.45/(7.2*pi);%截断频率对应的时间频率
S=0.25*0.0339*(fc)^2*100/3.6./(2*pi*f.^2.*(f.^2+(fc)^2));%美国六级谱对应的公式
Sx=[t,S,m,0];%补全采样点
k=Nr/2+1;
dfai=rand(1,k)*2*pi;%角度在0~2pi间均匀分布
Y=exp(dfai*i);%独立相位序列
Xk1=Nr*Y.*sqrt(Sx*df);%频谱
Xk2=Xk1(2:k-1);
Xk3=rot90(real(Xk2),2);
Xk4=rot90(imag(Xk2),2);
Xk5=Xk3-i*Xk4;
Xk=[Xk1,Xk5];%补全频谱
x=ifft(Xk,Nr);%傅立叶逆变换求时域样本




%将轨道不平顺的时域数据作成激励做准备，即送入 to Workspace （产生数组ex_USA6，然后在simulink中使用From Workspace，输入ex_USA6即可）
%送入  to File  （产生USA6_ex.mat文件，然后在Simulink中使用From File，输入USA6_ex.mat）
%这两种方法产生的激励源均为离散时间序列，因此仿真的Simulink模型为离散系统，或者混合系统
%仿真设置中，算法应该采用solver 'VariableStepDiscrete' 而不是 solver'ode45'

time=[0:delta:delta*(Nr-1)];                %时间域模拟时间采样后的离散时间序列

ex_USA6=[time',x'];    %产生数组ex_USA6，然后在simulink中使用From Workspace，输入ex_USA6即可作为激励
%傅里叶逆变换后的离散路面不平顺数据需要倒置
%注意需要倒置
%由于要想此程序产生的序列作为轨道不平顺激励源时，From Workspace 
%使用的数据为列的形式，即第一列为事件（如时间序列），
%第二列第三列等为输出的数据序列，因此,time与x均需要倒置


save ex_USA6   %%%%%%%%%%%%%%%%%%%%%%%%%%%将数据 ex_USA6 保存在 USA6_ex.mat文件中，若使用模块From Workspace仿真，需要双击文件USA6_ex.mat或者在matlab命令窗口输入load ex_USA6，将数据读入内存空间 Workspace，使From Workspace有数据可读

USA6=[time;x];  %From File 使用的数据均为行的形式，与From Workspace相反，即第一行为事件（如时间行），
%第二行第三行等为输出的数据序列，因此，time与x不需要倒置，中间用“；”  隔开，表示开启下一行

save  USA6_ex  USA6  %将数据USA6 保存为 USA6_ex.mat文件（相当于to File ）， 产生USA6_ex.mat文件，然后在Simulink中使用From File，输入USA6_ex.mat也可产生激励可
%该指令也可以为 save  D:\USA6_ex  USA6   相当于设置保持路径为D盘下的USA6_ex.mat

l=size(time);
tic=l(1,2)/10000;   %返回该激励的离散采样时间长度


%下面为画图程序

figure
plot(0.0001*(1:64:Nr),x(1:64:Nr));%画时域样本
title('模拟的时间序列','fontsize',16)
xlabel('时间/s','fontsize',16)
ylabel('高低不平顺幅值/cm','fontsize',16)
set(gca,'Fontname','Times New Roman','fontsize',16);
Sxn=(abs(Xk)/Nr).^2/df;
Kn=No+1:No+Nf;
f=Kn*df;


figure;
loglog(f,Sxn(Kn),'b-.');%绘制由时域样本得到的功率谱曲线
hold on
f=linspace(fl,fu,Nf);
fc=82.45/(7.2*pi);%截断频率对应的时间频率
S=0.25*0.0339*(fc)^2*100/3.6./(2*pi*f.^2.*(f.^2+(fc)^2));
loglog(f,S);%绘制解析谱线
grid on
title('不平顺解析值与模拟值的比较','fontsize',16)
xlabel('频率/Hz','fontsize',16)
ylabel('功率谱密度/(cm2/Hz)','fontsize',16)
set(gca,'Fontname','Times New Roman','fontsize',16);
legend('模拟值','解析值')

