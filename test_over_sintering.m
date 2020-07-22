load('data5_A_267000-271950.mat')
y = YY(1, :);
y = [0 y];

y(21:28) = linspace(450, 300, 8);  % 造数据，如果21号风箱最高会出错
y(20) = 300;
% y(22:28) = linspace(450, 300, 7);  % 造数据，如果21号风箱最高会出错

alfs   =     [55,58,61,64,67,70,72,74]+2;
bets   =     [58,61,64,67,70,72,74,76]+2;

[~, zb0] = max(y);
[zb_hat, TB_hat, ~, BET, ~, ~] = BTP_Search((alfs(zb0-20)+bets(zb0-20))/2, y);
plot_BTP(y, zb_hat, TB_hat, BET);

