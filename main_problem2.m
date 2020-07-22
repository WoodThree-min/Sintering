%% 
% data = load('问题数据/异常数据.txt');
% curve = load('问题数据/对应曲线.txt');

data = load('问题数据2/原始数据.txt');
curve = load('问题数据2/计算结果.txt');

%%
box_width = [2 2 2 3*ones(1, 22) 2 2 2];
left_boundary = cumsum([0 2 2 2 3*ones(1, 22) 2 2]);
right_boundary = cumsum(box_width);
boundary = [0 right_boundary];
pos = (left_boundary + right_boundary) / 2;
nPartition = 10;

xx1 = linspace(pos(21), right_boundary(25), 5*nPartition+1-nPartition/2);  % 到21左边界不好看，改为中点
xx2 = linspace(left_boundary(26), right_boundary(28), 3*nPartition+1);
xx = [xx1 xx2(2:end)];

flag_plot = 0;
for i = 250:-1:1
    if min(curve(i, 76)) < 0
        i
        % curve(i, :)
        figure(101); cla; subplot(2, 1, 1); hold on; plot(pos, data(i, :), 'o'); plot(xx, curve(i, :), '.');
        
        figure(102); cla;
        [YY1, YY2, Predict1, Predict2] = call(data(i, :), data(i, :), 0, 0, 0, 1);
        % break
        % drawnow
        pause
    end
end

