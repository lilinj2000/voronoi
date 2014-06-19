function fortuneAlgorithm(site_points, axis_scaling)
% FORTUNEALGORITHM demonstrate the seeping line algorithm for Voronoi
% diagram
% Function FORTUNEALGORITHM demonstrate the seeping line algorithm on
% voronoi diagram which site points is p, and the 2-D x- and y- axes is r
% 
% Define variables:
% site_points        -- site points from input
% site_points_event  -- site point events by y- axes descending
% axis_scaling       -- the scaling for the x- and y-axes    
% circle_events      -- circle events by y- axes descending 
% arc_list           -- arc tree, by x- axes ascending
% Record Of revisions:
% Date Programmer Description of change
% ======== ============== ========================
% 6/16/2014 Linjiang Li Original code
% 6/17/2014 LinJIang Li handleSiteEvent
% 6/19/2014 LinJiang Li update the var name

% the site point events is site points in descending order on y axes
[~, idx] = sort([site_points.y], 'descend');
site_point_events = site_points(idx);

% init circle event
circle_events = struct([]);

% init arc tree
arc_list = struct([]);

% init vertex
v = struct([]);

% go through all the site events
while length(site_point_events)>0
    while length(circle_events)> 0
       if  circle_events(1).y>site_point_events(1).y
           % handle circle event
           [circle_events, arc_list, v] = handleCircleEvent(circle_events(1), site_points, axis_scaling, arc_list, v);
       else
           break;
       end
    end
    
%     handle site event
    [circle_events, arc_list] = handleSiteEvent(site_point_events(1), site_points, axis_scaling, arc_list, v);
%      remove the site point event
    site_point_events(1) = [];
end

end


function [circle_events, arc_list] = handleSiteEvent(p, site_points, axis_scaling, arc_list_input, v)

arc_list = arc_list_input;

% add the site point arc to tree
if length(arc_list)~=0
    arc_list(length(arc_list)+1) = p;
else
    arc_list = p;
end

% sort the arc tree by x- ascending
[~, idx] = sort( [arc_list.x] );
arc_list = arc_list(idx);

% check circle events
circle_events = checkCircleEvents(arc_list);

% sort the circle events by y- descending
if length(circle_events)~=0
    [~, idx] = sort( [circle_events.y], 'descend');
    circle_events = circle_events(idx);
end

showFigure(site_points, axis_scaling, p.y, arc_list, [], v);

waitforbuttonpress;

end

function [circle_events, arc_list, V] = handleCircleEvent(c, site_points, axis_scaling, arc_list_input, v)

arc_list = arc_list_input;
V = v;

showFigure(site_points, axis_scaling, c.y, arc_list, c, v);
waitforbuttonpress

% remove the site points from arc tree
for ii = 1:length(arc_list)
    if [arc_list(ii).x arc_list(ii).y]==[c.p.x c.p.y]
%         remove the arc from tree
        arc_list(ii) = [];
        break;
    end
end

% add new vertex for the voronoi
% the center of the circle, just the vertex
if length(V)~=0
    V(length(V)+1)  = c.center;
else
    V = c.center;
end

% check circel events
circle_events = checkCircleEvents(arc_list);

end

function circle_events = checkCircleEvents(arc_list)

circle_events = struct([]);

    for ii = 2 : length(arc_list) - 1

        c = checkCircle(arc_list(ii-1), arc_list(ii), arc_list(ii+1));

    %     found the circel event
        if length(c)~=0
    %         add the circle event
            if length(circle_events)~=0
                circle_events(length(circle_events) +1) = c;
            else
                circle_events = c;
            end
        end
    end

end


function circle_event = checkCircle(pii, pjj, pkk)

% set the circle_event is empty
circle_event = struct([]);

% x1 + x2
A = pii.x + pjj.x;
% x2 + x3
B = pjj.x + pkk.x;

% x2 - x1
C = pjj.x - pii.x;
% x3 - x2
D = pkk.x - pjj.x;

% y1 + y2
E = pii.y + pjj.y;
% y2 + y3
F = pjj.y + pkk.y;

% y2 - y1
G = pjj.y - pii.y;
% y3 - y2
H = pkk.y - pjj.y;

if G.*D == H.*C
%   the k is equal
%   no circle
    return ;
end

AA = [ 2.*C, 2.*G; 2.*D, 2.*H];
BB = [ A.*C + E.*G; B.*D + F.*H];

XX = AA\BB;

center.x = XX(1,1);
center.y = XX(2,1);
% center.x = (D.*B./(2.*H) - C.*A./(2.*G) + F./2 - E./2)./(D./H - C./G);
% center.y = ((A-B)./2 + G.*E./(2.*C) - H.*F./(2.*D))./(G./C - H./D);

% add eps 0.001
radius = sqrt((center.x - pii.x).^2 + (center.y - pii.y).^2) + 0.001;

c.y = center.y-radius;
c.center = center;
c.radius = radius;
c.p = pjj;

circle_event = c;

end

function showFigure(site_points, axis_scaling, y, arc_list, circle_event, v)

% one new figure
figure;

hold on;

% plot the site points
plot([site_points.x], [site_points.y], 'bx');

% plot the sweep line
x = linspace(axis_scaling.xmin, axis_scaling.xmax, 1000);
plot(x, y);

% plot the arc
for ii = 1 : length(arc_list)
    if arc_list(ii).y~= y
        drawParabola(arc_list(ii), y, axis_scaling);
    end
end

% plot the circle
if length(circle_event)~=0
    
    radius = circle_event.radius;
    center = circle_event.center;
    
    x = linspace(center.x- radius, center.x + radius, 1000);
    y1 = sqrt(radius.^2-(x-center.x).^2) + center.y;
    y2 = -sqrt(radius.^2-(x-center.x).^2) + center.y;
    plot(x, y1, 'b-', x, y2, 'b-');
end

% plot the vertex
if length(v)~=0
    plot( [v.x], [v.y], 'r^');
end


% set the x- and y- axes
axis([axis_scaling.xmin axis_scaling.xmax axis_scaling.ymin axis_scaling.ymax]);

hold off;

end