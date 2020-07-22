load('data5_A_267000-271950.mat')
y = YY(1, :);
y = [0 y];

y(21:28) = linspace(450, 300, 8);  % �����ݣ����21�ŷ�����߻����
y(20) = 300;
% y(22:28) = linspace(450, 300, 7);  % �����ݣ����21�ŷ�����߻����

alfs   =     [55,58,61,64,67,70,72,74]+2;
bets   =     [58,61,64,67,70,72,74,76]+2;

[~, zb0] = max(y);
[zb_hat, TB_hat, ~, BET, ~, ~] = BTP_Search((alfs(zb0-20)+bets(zb0-20))/2, y);
plot_BTP(y, zb_hat, TB_hat, BET);

