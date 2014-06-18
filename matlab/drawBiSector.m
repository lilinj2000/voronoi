function drawBiSector(p1, p2)

xMin = min(p1(1,1), p2(1,1));
xMax = max(p1(1,1), p2(1,1));

x = linspace(xMin, xMax, 100);

% x1 + x2
A = p1(1,1) + p2(1,1); 
% y1 + y2
B = p1(1,2) + p2(1,2);
% x2 - x1
C = p2(1,1) - p1(1,1);
% y2 - y1
D = p2(1,2) - p1(1,2);

K = D ./ C;

y = -1./K.*(x-A./2) + B./2;

% y = (p2(1,1)^2 + p2(1,2)^2 - (p1(1,1)^2 + p1(1,2)^2) - 2.*(p2(1,1)-p1(1,1)).*x)./(2.*(p2(1,2)-p1(1,2)));

plot(x, y, 'r');
end