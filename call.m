function [YY1, YY2, Predict1, Predict2] = call(Y1_, Y2_, freq1, freq2, speed, flag_modify, stamp, filename)
% y1: 南侧风箱温度
% y2: 北侧风箱温度
% freq1: 南抽风频率
% freq2: 北抽风频率
% speed: 机速
% stamp: 时间戳
%返回值
% YY1: 南侧BTP温度
% YY2: 北侧BTP温度

%example
% load ys; Y1 = ys;
% load nys; Y2 = YY;
% [YY1, YY2, btp1, btp2, errorcode, stamp] = main(Y1, Y2, 43, 43, 1.6, 1568108967602)

if nargin<2
    % load ys; Y1_ = ys; Y1_(21:28) = linspace(300, 450, 8); % Y1_(18) = 580;
    % load nys; Y2_ = YY; Y2_(2) = -3;
    
    Y1_ = [0.0, 107.4, 108.1, 100.7, 85.9, 80.4, 78.5, 76.6, 77.5, 75.6, 75.7, 74.8, 75.5, 136.4, 74.2, 73.1, 73.7, 73.9, 80.0, 104.3, 152.6, 257.1, 301.3, 363.2, 423.6, 450.8, 480.9, 152.6];
    % Y1_ = [0, 95.8, 94, 87.5, 72.9, 67.4, 66.5, 63.9, 64.7, 62.8, 63.7, 63.6, 64.3, 124.3, 66.9, 68.6, 69.3, 82.5, 143.5, 263.8, 372.7, 437.8, 480.3, 473.1, 504.6, 403.6, 361.9, 310];
    Y2_ = [0.0, 35.2, 107.3, 96.9, 85.0, 76.7, 78.0, 76.7, 75.5, 77.0, 75.0, 76.3, 73.4, 73.4, 75.5, 74.3, 76.6, 76.1, 82.6, 96.6, 158.5, 229.2, 265.2, 315.3, 386.9, 517.2, 404.2, 323.0];
    % Y1_(17) = 599;

    freq1 = 43; freq2 = 43;
    speed = 1.6;
    flag_modify = 1;
    % stamp = 1568108967602;
end

if ~exist('flag_modify', 'var') || isempty(flag_modify)
    flag_modify = 0;
end

if ~exist('stamp', 'var') || isempty(stamp)
    s = datenum(datetime());  % matlab
else
    s = convert_stamp(stamp);  % Java timestamp -> matlab
end

if ~exist('filename', 'var') || isempty(filename)
    filename = [datestr(s, 'yyyymmddHHMMss') '.png'];
end

if length(Y1_) == 27
    Y1_ = [0 Y1_];
end

if length(Y2_) == 27
    Y2_ = [0 Y2_];
end

global pos;
global left_boundary;
global right_boundary;
global boundary;
global nPartition;
global H;
global zb1;
global zb2;
global p;
global LENGTH;
global data1;
global data2;

% 数据清洗，把修正的数据用不同颜色表示？
Y1 = data_clean(Y1_);
Y2 = data_clean(Y2_);

if flag_modify
    Y1(2:28) = pre_process(Y1(2:28), pos(2:28)-2);  % 结果依赖pos
    Y2(2:28) = pre_process(Y2(2:28), pos(2:28)-2);

%     Y1(2:28) = pre_process_ccx(Y1(2:28), pos(2:28)-2);  % 结果依赖pos
%     Y2(2:28) = pre_process_ccx(Y2(2:28), pos(2:28)-2);
end

% xx1 = linspace(left_boundary(21), right_boundary(25), 5*nPartition+1);   % 到21左边界
xx1 = linspace(pos(21), right_boundary(25), 5*nPartition+1-nPartition/2);  % 到21左边界不好看，改为中点
xx2 = linspace(left_boundary(26), right_boundary(28), 3*nPartition+1);
xx = [xx1 xx2(2:end)];

% 南侧
[TMAX, IMAX] = max(Y1);
if IMAX<21
    ZB1 = pos(IMAX);
    TB1 = TMAX;
    % YY1 = NaN;  % 这里应该改成其他的东西，直接插值或者样条？
    YY1 = interp1(pos(21:28), Y1(21:28), xx, 'pchip');  % 
else
    ZB0 = pos(IMAX);
    [ZB1, TB1, ~, BET1, ~, ~] = BTP_Search(ZB0, Y1);   % Evaluating Q-Function;
    TB1 = max(TB1, TMAX);  % TB1可能小于TMAX

    % idx = get_box_no(ZB1);
    x1 = xx(xx<=ZB1);
    x2 = xx(xx>ZB1);

    a1 = BET1(1); b1 = BET1(2);
    a2 = BET1(3); b2 = BET1(4);
    yl = cal_temperature(x1, ZB1, TB1, a1, b1, 1.5);  % 左侧delta全部当作1.5
    yr = cal_temperature(x2, ZB1, TB1, a2, b2, 1);    % 右侧delta全部当作1
    YY1 = [yl yr];    
end
%% 如果结果为负
% [TMIN, IMIN] = min(YY1);
% if YY1(end)<0
%     if IMAX == 27
%         % 
%     end
% end


%% 北侧
[TMAX, IMAX] = max(Y2);
if IMAX<21
    ZB2 = pos(IMAX);
    TB2 = TMAX;
    % YY2 = NaN;
    YY2 = interp1(pos(21:28), Y2(21:28), xx, 'pchip');  % 
else
    ZB0 = pos(IMAX);

    [ZB2, TB2, ~, BET2, ~, ~] = BTP_Search(ZB0, Y2);   % Evaluating Q-Function;
    TB2 = max(TB2, TMAX);  % TB2可能小于TMAX

    % idx = get_box_no(ZB1);
    x1 = xx(xx<=ZB2);
    x2 = xx(xx>ZB2);

    a1 = BET2(1); b1 = BET2(2);
    a2 = BET2(3); b2 = BET2(4);
    yl = cal_temperature(x1, ZB2, TB2, a1, b1, 1.5);  % 左侧delta全部当作1.5
    yr = cal_temperature(x2, ZB2, TB2, a2, b2, 1);    % 右侧delta全部当作1
    YY2 = [yl yr];
end

%% 如果结果为负
% [TMIN, IMIN] = min(YY2);


%% 这里处理预测
p = p + 1;
if p>LENGTH  % 超出记忆长度，减半
    HALF = floor(LENGTH/2);
    zb1(1:HALF) = zb1(LENGTH-HALF+1:LENGTH);
    zb2(1:HALF) = zb2(LENGTH-HALF+1:LENGTH);
    p = HALF + 1;
    
    data1(1:HALF, :) = data1(LENGTH-HALF+1:LENGTH, :);
    data2(1:HALF, :) = data2(LENGTH-HALF+1:LENGTH, :);
end

data1(p, :) = Y1_;
data2(p, :) = Y2_;

zb1(p) = ZB1;
zb2(p) = ZB2;

[Z11, Z12, sigma11, sigma12] = predict_btp_y(zb1, p);
[Z21, Z22, sigma21, sigma22] = predict_btp_y(zb2, p);

Predict1 = [TB1, ZB1, Z11, Z12, sigma11, sigma12];
Predict2 = [TB2, ZB2, Z21, Z22, sigma21, sigma22];


%% 打印图片
% 计算ZB所在的风箱号，左边和右边分别画
subplot(2, 1, 1); cla; hold on;
plot(xx, YY1, '.b', 'markersize', 5);      % 拟合的曲线
plot(pos, Y1_, '.r', 'markersize', 20)  % 原始数据
% plot(pos, Y1, '.r', 'markersize', 20)  % 修正后的数据
plot(ZB1, TB1, '.r', 'markersize', 30);  % BTP
plot([boundary; boundary], [zeros(1,29); 100*ones(1, 29)], 'b')  % 风箱边界
% plot_BTP(Y1, ZB1, TB1, BET1);

im1 = (Y1_ ~= Y1);
if ~isempty(im1)
    plot(pos(im1), Y1(im1), 'y.', 'markersize', 20);  % 修正的数据
end
xticks(pos)

subplot(2, 1, 2); cla; hold on;
plot(xx, YY2, '.b', 'markersize', 5);
plot(pos, Y2_, '.r', 'markersize', 20)
% plot(pos, Y2, '.r', 'markersize', 20)
plot(ZB2, TB2, '.r', 'markersize', 30);
plot([boundary; boundary], [zeros(1,29); 100*ones(1, 29)], 'b')
% plot_BTP(Y2, ZB2, TB2, BET2)

im2 = Y2_ ~= Y2;
if ~isempty(im2)
    plot(pos(im2), Y2(im2), 'y.', 'markersize', 20);  % 修正的数据
end
xticks(pos)


% print(H, '-dpng', filename)
% print(filename, '-dpng');

