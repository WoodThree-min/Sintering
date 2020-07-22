function y = cal_temperature(z, ZB, TB, a, b, delta)

y = a*(z.^2 + 1/3 * delta^2 - ZB^2) + b*(z-ZB) + TB;

% c1 = TB + a1 * ( - ZB^2) + b1 *(-ZB);   % Wei: Ignore a * (delta ^2 /3);
% y1 = a1 * Z.^2 + b1 * Z + c1;

