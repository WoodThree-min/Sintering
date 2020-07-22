function y = data_clean(y0)
LIMIT = 600;
y = y0;
% im = zeros(size(y0));  % index of modified
% s = 0;
for i = 2:18  % 
    if y(i)<0 || y(i) > LIMIT
        y(i) = (y(i-1) + y(i+1)) / 2;
        % s = s + 1;
        % im(s) = i;
    end
end

% im = im(1:s);

