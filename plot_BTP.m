function plot_BTP(y, ZB, TB, BET)
pos = [1.0 3.0 5.0 7.5 10.5 13.5 16.5 19.5 22.5 25.5 28.5 31.5 34.5 37.5 40.5 43.5 46.5 49.5 52.5 55.5 58.5 61.5 64.5 67.5 70.5 73.0 75.0 77.0];

xx = linspace(pos(21), pos(28), 101);
x1 = xx(xx<=ZB);
x2 = xx(xx>ZB);

a1 = BET(1); b1 = BET(2);
a2 = BET(3); b2 = BET(4);
yl = cal_temperature(x1, ZB, TB, a1, b1, 1.5);  % 左侧delta全部当作1.5
yr = cal_temperature(x2, ZB, TB, a2, b2, 1);    % 右侧delta全部当作1
yy = [yl yr];

box_width = [2 2 2 3*ones(1, 22) 2 2 2];
% left_boundary = cumsum([0 0 2 2 3*ones(1, 22) 2 2]);
right_boundary = cumsum(box_width);
boundary = [0 right_boundary];

plot(pos, y, '.-', 'markersize', 20)  % 原始数据
hold on;
plot(ZB, TB, '.r', 'markersize', 30);  % BTP
plot(xx, yy);
plot([boundary; boundary], [zeros(1,29); 100*ones(1, 29)], 'b')  % 风箱边界

xticks(pos)
hold off;
