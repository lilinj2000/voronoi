function drawParabola(p, line, axis_scaling)

x = linspace(axis_scaling.xmin, axis_scaling.xmax, 1000);

% (y-lineY).^2 = (x-pii(1,1)).^2 + (y-pii(1,2)).^2;

y = ((x - p.x).^2 + p.y.^2 - line.^2)./(2.*(p.y-line));

plot(x, y, 'k-');

% plot(x, lineY);

end