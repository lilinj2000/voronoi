function drawParabola(pii, lineY, r)

x = linspace(r(1), r(2), 1000);

% (y-lineY).^2 = (x-pii(1,1)).^2 + (y-pii(1,2)).^2;

y = ((x - pii(1,1)).^2 + pii(1,2).^2 - lineY.^2)./(2.*(pii(1,2)-lineY));

plot(x, y, 'k-');

% plot(x, lineY);

end