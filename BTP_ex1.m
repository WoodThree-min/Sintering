clear;

% load ys.mat
% load ys.mat;YY=ys;
% load nys.mat;
% load YY.mat;
% load YY1.mat;
% load('yy_4_250_3550.mat', 'yy_4_250_3550')
% C =yy_4_250_3550;
% load('temp_328300_333220.mat')
% C=temp_328300_333220;
% load data_cleansed.mat;
% C = (A+B) / 2;
% load data_cleansed5.mat;
% C = (A5+B5) / 2;
% load('A_side_cleaned_data_70450-78400.mat')
% C=A2;
load('data5_A_88850-95000.mat')
C=YY;
% load('data5_A_267000-271950.mat')
% C=YY;
warning off;
nrC = size(C,1);C=[zeros(nrC,1),C];
alfs   =     [55,58,61,64,67,70,72,74]+2;
bets   =     [58,61,64,67,70,72,74,76]+2;
delts  =     bets-alfs;
pos=[1 3 5 4.5+3*(1:22) 73 75 77];
[zb_hat,TB_hat,zb0] = deal(zeros(1,nrC));
m=12;nzb_hat=zeros(1,nrC-m+1);
%%
tic;
for i = 1:nrC
  [~,zb0(i)]=max(C(i,:));
   [zb_hat(i), TB_hat(i), ~,~,~] = BTP_Search((alfs(zb0(i)-20)+bets(zb0(i)-20))/2,C(i,:));
end
toc;
for i = 1:(nrC-m+1)
    nzb_hat(i) = mean(zb_hat(i:i+m-1));
end
 figure(10);plot(nzb_hat,'.');
 hold on;
 plot(pos(zb0),'r-');
 ylim([69,76]);

 