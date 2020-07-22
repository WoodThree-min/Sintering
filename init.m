function init()

global pos;
global left_boundary;
global right_boundary;
global boundary;
% global alfs;
% global bets;
% global delts;
global nPartition;  % 每个风箱的分段数
global H;

box_width = [2 2 2 3*ones(1, 22) 2 2 2];
left_boundary = cumsum([0 2 2 2 3*ones(1, 22) 2 2]);
right_boundary = cumsum(box_width);
boundary = [0 right_boundary];
pos = (left_boundary + right_boundary) / 2;
nPartition = 10;

% alfs   =     [55,58,61,64,67,70,72,74]';
% bets   =        [58,61,64,67,70,72,74,76]';
% delts  =     0.5 * ( bets - alfs );

global zb1;
global zb2;
global LENGTH;
LENGTH = 12*20;
zb1 = zeros(1, LENGTH);
zb2 = zeros(1, LENGTH);

global data1;
global data2;

nBox = 28;
data1 = zeros(LENGTH, nBox);  % 记录数据
data2 = zeros(LENGTH, nBox);

global p;  % 记录位置
p = 0;

H = figure('visible', 'off'); clf;
% H = figure('visible', 'on'); clf;
set(H, 'position', [0 0 1400 800])

warning('off')
