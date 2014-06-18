function fortuneAlgorithm(p, r)
% FORTUNEALGORITHM demonstrate the seeping line algorithm for Voronoi
% diagram
% Function FORTUNEALGORITHM demonstrate the seeping line algorithm on
% voronoi diagram which site points is p, and the 2-D x- and y- axes is r
% 
% Define variables:
% p     -- site points
% r     -- the scaling for the x- and y-axes    
% c     -- circle events by y- axes descending
%           y, centerX, centerY, radius, piiX, piiY
% t     -- arc tree, by x- axes ascending
% v     -- voronoi vertex
% e     -- voronoi edge
% Record Of revisions:
% Date Programmer Description of change
% ======== ============== ========================
% 6/16/2014 Linjiang Li Original code
% 6/17/2014 LinJIang Li handleSiteEvent

% sort the site points which in descending order on y axes
p = sortrows(p, -2);

% init circle event
c = [];

% init arc tree
t = [];

% init vertex
v = [];

% go through all the site events
for ii = 1:size(p,1)
    while size(c, 1)> 0
       if  c(1, 1)>p(ii, 2)
           % handle circle event
           [c, t, v] = handleCircleEvent(c(1, :), p, r, t, v);
       end
    end
    
%     handle site event
    [c, t] = handleSiteEvent(p(ii, :), p, r, t, v);
  
end

% axis(r);

% hold off

end


function [C, T] = handleSiteEvent(pii, p, r, t, v)

C = [];

T = t;

% add the site point arc to tree
T(size(T,1)+1, :) = pii;
% sort the arc tree by x- ascending
T = sortrows(T, 1);

% check circle events
for ii = 2 : size(T,1)-1
    
    c = checkCircle(T(ii-1, :), T(ii, :), T(ii+1, :));
    
%     found the circel event
    if size(c,1)~=0
%         add the circle event
        C(size(C, 1)+1, :) = c;
    end
end

% sort the circle events by y- descending
if size(C,1)~=0
    C = sortrows(C, -2);
end

showFigure(p, r, pii(1,2), T, [], v);

waitforbuttonpress;

end

function [C, T, V] = handleCircleEvent(c, p, r, t, v)

C = [];
T = t;
V = v;

showFigure(p, r, c(1,1), t, c, v);
waitforbuttonpress

% remove the site points from arc tree
for ii = 1:size(T,1)
    if T(ii, :)==[c(1,5), c(1,6)]
%         remove the arc from tree
        T(ii, :) = [];
        break;
    end
end

% add new vertex for the voronoi
% the center of the circle, just the vertex
V(size(V, 1)+1, :) = [c(1,2), c(1,3)];

end

function c = checkCircle(pii, pjj, pkk)

% init circle center & radius
c = [];

% x1 + x2
A = pii(1,1) + pjj(1,1);
% x2 + x3
B = pjj(1,1) + pkk(1,1);

% x2 - x1
C = pjj(1,1) - pii(1,1);
% x3 - x2
D = pkk(1,1) - pjj(1,1);

% y1 + y2
E = pii(1,2) + pjj(1,2);
% y2 + y3
F = pjj(1,2) + pkk(1,2);

% y2 - y1
G = pjj(1,2) - pii(1,2);
% y3 - y2
H = pkk(1,2) - pjj(1,2);

if (G==0 & H==0) | (C==0 & D==0) | (G.*D == H.*C)
%   k is 0, or the k is equal
%   no circle
    return ;
end

center(1, 1) = (D.*B./(2.*H) - C.*A./(2.*G) + F./2 - E./2)./(D./H - C./G);
center(1, 2) = ((A-B)./2 + G.*E./(2.*C) - H.*F./(2.*D))./(G./C - H./D);

radius = sqrt((center(1, 1) - pii(1,1)).^2 + (center(1,2) - pii(1,2)).^2) + 0.001;

c = [ center(1,2)-radius, center, radius, pjj ];

end

function showFigure(p, r, y, t, c, v)

% one new figure
figure;

hold on;

% plot the site points
plot(p(:,1), p(:, 2), 'bx');

% plot the sweep line
x = r(1) : 0.01 : r(2);
plot(x, y);

% plot the arc
for ii = 1 : size(t,1)
    if t(ii, 2)~= y
        drawParabola(t(ii, :), y, r);
    end
end

% plot the circle
if size(c,1)~=0
    radius = c(1,4);
    x = linspace(c(1,2)- radius, c(1,2) + radius, 1000);
    y1 = sqrt(radius.^2-(x-c(1,2)).^2) + c(1,3);
    y2 = -1.*sqrt(radius.^2-(x-c(1,2)).^2) + c(1,3);
    plot(x, y1, x, y2);
end

% plot the vertex
if size(v, 1)~=0
plot(v(:, 1), v(:, 2), 'r^');
end


% set the x- and y- axes
axis(r);

hold off;

end