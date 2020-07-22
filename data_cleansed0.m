function [datad]=data_cleansed0(data,m)
pos = [1, 3, 5.5, 5.5+cumsum(ones(1,21)*3), 71, 73, 75]; % 28个风箱中点。第一个风箱为点火位置，不计数，故设置为0.


% V型一修
for i=1:size(data,1)
   if data(i,24)<mean(data(i,[23,25])) 
        data(i,24)=(data(i,23)+data(i,25))/2;
   elseif data(i,25)<mean(data(i,[24,26]))
        data(i,25)=(data(i,24)+data(i,26))/2;
   end
end 
pre_data=data(size(data,1),:);

%% 二补
% 补24、25 的v型误差

if size(data,1) <= m      
    data_avg=mean(data,1);
else
    data_avg=mean(data(size(data,1)-m:size(data,1)-1,:));
end


if pre_data(24) <= mean(pre_data([23,25]))
    yb=data_avg([22,23,25]);
    z=pos([22,23,25]);
    p=polyfit(z,yb,2);
    a=p(1);b=p(2);c=p(3);
    z_m0=pos(24);
    yb_m0=z_m0^2*a+z_m0*b+c;
    dif=yb_m0-data_avg(24);
    pre_data(24) = pre_data(24)*(1+dif/yb_m0);
elseif pre_data(25) <= mean(pre_data([24,26]))
    yb=data_avg([22,23,24,26]);
    z=pos([22,23,24,26]);
    p=polyfit(z,yb,2);
    a=p(1);b=p(2);c=p(3);
    z_m0=pos(25);
    yb_m0=z_m0^2*a+z_m0*b+c;
    dif=yb_m0-data_avg(25);
    pre_data(25) = pre_data(25)*(1+dif/yb_m0);
end




datad=pre_data;

end

