function drawCircle(center, radius)
 
x = linspace(center.x- radius, center.x + radius, 1000);
y1 = sqrt(radius.^2-(x-center.x).^2) + center.y;
y2 = -sqrt(radius.^2-(x-center.x).^2) + center.y;
plot(x, y1, 'b-', x, y2, 'b-');

end