function drawLine(start_p, end_p)

% y - y1 = (y2-y1)/(x2-x1)*(x-x1)

if start_p.x~=end_p.x
    k = (end_p.y-start_p.y)./(end_p.x-start_p.x);
    
    x = linspace(start_p.x, end_p.x, 1000);
    y = start_p.y + k.*(x-start_p.x);
else
    x = start_p.x;
    y = linspace(start_p.y, end_p.y, 1000);
end

plot(x, y, 'r-');

end