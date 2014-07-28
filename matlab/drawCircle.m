function [X, Y] = drawCircle(x, y, radius)

angle = linspace(0, 2.*pi, 1000);

X = x + radius.*cos(angle);
Y = y + radius.*sin(angle);

end