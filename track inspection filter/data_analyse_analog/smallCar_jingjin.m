% =========================================================================
%
%                  ����
%
% =========================================================================
%
%��(C)2019-2022 ������ѧ�о�Ժ-������
%   �汾��V1.0
%   ���ڣ�2020�� 9��16��
%   ���ߣ�s.m.
%--------------------------------------------------------------------------
%  ���ܣ� 1.�����õ��Ǿ�������ݣ�������ʱû�в����ļ�
%        2.
%       3. 
%--------------------------------------------------------------------------

clear;
data = load('gaodi_xiaoche.txt');
gpproLTrolleyOri = data(:,2);
xSmallCar = data(:,1);
ll = xSmallCar(2:end)-xSmallCar(1:end-1);
% figure;plot(ll);
b = load('filter_120m.mat');
len = (length(b.Num)-1)/2;
out = conv(b.Num,gpproLTrolleyOri);

gpproLTrolley = out(len+1:end-len);
%% �����Ƴ���0.65mת����Ϊ0.25m������ͬ��
gpproLTrolley25m = [];
for i = 1:448490
% for i = 1:4484
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = gpproLTrolley(p1) + p2/0.65*( gpproLTrolley(p1+1) - gpproLTrolley(p1) );
end
figure;plot(0:0.25:(length(gpproLTrolley25m)-1)*0.25 , gpproLTrolley25m,'r');
hold on;plot(0:0.65:(length(gpproLTrolley)-1)*0.65 , gpproLTrolley,'--k');

%%
 data2 = load('gaodi_guijian.txt');
 gpproGjOri = data2(:,6);
 xGj = data2(:,1);
figure;plot(xGj,gpproGjOri);

%% ������ݶ����ˣ��������
hold on;plot(xSmallCar(81572 : 81572+1538) , gpproLTrolley(81572 : 81572+1538));

%% ��һ�ַ�ʽ����
gpproLTrolley25m = [];
listTp = gpproLTrolley(81572 : 81572+1538);
for i = 1:3998
    p1 = floor(i*0.25/0.65) + 1;
    p2 = mod(i*0.25 , 0.65);
    gpproLTrolley25m(i) = listTp(p1) + p2/0.65*( listTp(p1+1) - listTp(p1) );
end

%%

out = conv(gpproGjOri,gpproLTrolley25m(end:-1:1));
[~,q] = max(out);
q = q - length(gpproLTrolley25m);
figure;plot(gpproGjOri);hold on;plot([zeros(1,q) , gpproLTrolley25m])

%%
function out = find_index(xSmallCar,start_,end_)


end
