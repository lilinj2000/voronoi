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
% seg_list           -- segement list
% Record Of revisions:
% Date Programmer Description of change
% ======== ============== ========================
% 6/16/2014 Linjiang Li Original code
% 6/17/2014 LinJIang Li handleSiteEvent
% 6/19/2014 LinJiang Li update the var name
% 6/26/2014 LinJiang Li fix the site event & circle event process


% site_points { p(x,y) }
% site_points_event { p(x,y) }
% circle_events { y, center(x,y), radius, arc }
% arc_list { p, seg0(start_p, end_p), seg1(start_p, end_p), circle_event }
% seg_list { start_p, end_p, p}


% the site point events is site points in descending order on y axes
[~, idx] = sort([site_points.y], 'descend');
site_point_events = site_points(idx);

% init circle event
circle_events = struct([]);

% init arc tree
arc_list = struct([]);

% init seg list
seg_list = struct([]);

% go through all the site events
while ~isempty(site_point_events)
    while ~isempty(circle_events)
       if  circle_events(1).y>site_point_events(1).y
           % handle circle event
           [arc_list, circle_events, seg_list] = handleCircleEvent(circle_events(1), site_points, axis_scaling, arc_list, circle_events, seg_list);
       else
           break;
       end
    end
    
%     handle site event
    [arc_list, circle_events] = handleSiteEvent(site_point_events(1), site_points, axis_scaling, arc_list, circle_events, seg_list);
%      remove the site point event
    site_point_events(1) = [];
end

% handle the left circle events
while ~isempty(circle_events)
   % handle circle event
   [arc_list, circle_events, seg_list] = handleCircleEvent(circle_events(1), site_points, axis_scaling, arc_list, circle_events, seg_list);    
end

% finish all the arc list
seg_list = finishArcList(arc_list, seg_list, axis_scaling);

% last figure show
showFigure(site_points, axis_scaling, [], arc_list, [], seg_list);

end


function [arc_list, circle_events] = handleSiteEvent(p, site_points, axis_scaling, arc_list_input, circle_events_input, seg_list)

arc_list = arc_list_input;
circle_events = circle_events_input;

% add arc
[arc, arc_list] = addArc(p, arc_list, axis_scaling.ymax);

% [arc_list, circle_events] = updateSegAfterAddArc(p, arc_list, circle_events);

% check & update circle events
[arc_list, circle_events] = updateCircleEvents(arc, arc_list, circle_events);

showFigure(site_points, axis_scaling, p.y, arc_list, [], seg_list);

waitforbuttonpress;

end

function [arc_list, circle_events, seg_list] = handleCircleEvent(c, site_points, axis_scaling, arc_list_input, circle_events_input, seg_list_input)

arc_list = arc_list_input;
circle_events = circle_events_input;
seg_list = seg_list_input;

showFigure(site_points, axis_scaling, c.y, arc_list, c, seg_list);
waitforbuttonpress

% update seg
% seg_list = updateSegBeforeRemoveArc(c.arc, arc_list, seg_list, c.center);

% remove the site points from arc tree
[arc_list circle_events seg_list] = removeArc(c.arc, arc_list, circle_events, seg_list, c.center);

% remove the circle event
% [arc_list, circle_events] = removeCircleEvent(c, arc_list, circle_events);

% showFigure(site_points, axis_scaling, c.y, arc_list, [], seg_list);
% waitforbuttonpress

% add new vertex for the voronoi
% the center of the circle, just the vertex
% if length(V)~=0
%     V(length(V)+1)  = c.center;
% else
%     V = c.center;
% end

% check circel events
% [arc_list, circle_events] = updateCircleEvents(c.arc, arc_list, circle_events);

end

function [arc_list, circle_events] = updateCircleEvents(arc, arc_list_input, circle_events_input)

arc_list = arc_list_input;
circle_events = circle_events_input;

% check arc index
index = findArc(arc, arc_list);
if index>0
%     do update circle events
%  update index-2, index-1, index

    for ii = index-2:index
        if ii>0 && ii+2<=length(arc_list)
    %         remove circle event from arc_list & circle_events
            c = arc_list(ii+1).circle_event;
            
            if ~isempty(c)
                [circle_events] = removeCircleEvent(c, circle_events);
            end
            arc_list(ii+1).circle_event = [];
    
    %         has the circle events
            c = checkCircle(arc_list(ii).p, arc_list(ii+1).p, arc_list(ii+2).p);
            if ~isempty(c)
                c.arc = arc_list(ii+1);
                arc_list(ii+1).circle_event = c;
                % add the full arc info for the circle_event
                c.arc = arc_list(ii+1);
                
                circle_events = addCircleEvent(c, circle_events);
            end
        end
    end
end

end


function circle_event = checkCircle(pii, pjj, pkk)

% set the circle_event is empty
circle_event = struct([]);

% ensure edge pii-pjj on the left pii-pkk
if (pkk.y-pii.y).*(pjj.x-pii.x) - (pjj.y-pii.y).*(pkk.x-pii.x)>0
    return ;
end

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

function node = intersect(p, index, arc_list)

node = struct([]);

if p.y==arc_list(index).p.y
    % it's not intersect
    return ;
end

node0 = struct([]);
node1 = struct([]);

if index-1>0
    node0 = intersection(arc_list(index-1).p, arc_list(index).p, p.y);
end

if index+1<=length(arc_list)
    node1 = intersection(arc_list(index).p, arc_list(index+1).p, p.y);
end

% between node0 and node1, then intersect with this arc
if (isempty(node0) || node0.x<=p.x) && (isempty(node1) || node1.x>=p.x)
    
    node = intersection(arc_list(index).p, p, p.y);
end

end

function  node = intersection(p0, p1, l)

p = p0;
if p0.y==p1.y % paralle with x-
    node.x = (p0.x+p1.x)./2;
elseif p0.y==l % site event,  cross by p0
    node.x = p0.x;
    p = p1;
elseif p1.y==l % site evetn, cross by p1
    node.x = p1.x;
else % circle event
    A = p0.y - l;
    B = p1.y - l;
    
    AA = 1./A - 1./B ;
    BB = -2.*(p0.x ./ A - p1.x ./ B);
    CC = (p0.x.^2+p0.y.^2-l.^2)./A-(p1.x.^2+p1.y.^2-l.^2)./B;
    
    node.x = (-BB + sqrt(BB.^2-4.*AA.*CC))./(2.*AA);    
end

node.y = ((p.x-node.x).^2 + p.y.^2 - l.^2)./(2.*(p.y-l));

end


function [arc, arc_list] = addArc(p, arc_list_input, ymax)

arc_list = arc_list_input;

% new arc
arc.p = p;
arc.circle_event = [];
arc.seg0 = [];
arc.seg1 = [];

if isempty(arc_list)
    arc_list = arc;
    return ;
end

% check insection and insert the arc
for ii=1:length(arc_list)
    node = intersect(p, ii, arc_list);
    
    if ~isempty(node)
        
        % remove this arc circle event, for it's false alarm
%         if ~isempty(arc_list(ii).circle_event)
%             c = arc_list(ii).circle_event ;
%             arc_list(ii).circle_event = [];
%             circle_events = removeCircleEvent(c, circle_events);
%         end
        
        len = length(arc_list);
        arc_list(len+2:-1:ii+2) = arc_list(len:-1:ii);
        
        %update the arc seg
        arc.seg0.start_p = node;
        arc.seg1.start_p = node;
        
        arc_list(ii+1) = arc;
        
        %update the pre arc
        arc_list(ii).seg1.start_p = node;
        
        %update the post arc
        arc_list(ii+2).seg0.start_p = node;
        
        
        return ;
    end
end

% for the specail point, which no intersect
start_p.x = ymax;
start_p.y = (arc_list(ii).p.x+p.x)./2;

arc.seg0.start_p = start_p;

arc_list(ii).seg1.start_p = start_p;
arc_list(ii+1) = arc;

% add the site point arc to tree
% if ~isempty(arc_list)
%     arc_list(length(arc_list)+1) = arc;
% else
%     arc_list = arc;
% end

% sort the arc tree by x- ascending
% fetch all the points from arc list
% pp = [arc_list.p];
% according the pp index, update the arc_list index
% [~, idx] = sort( [pp.x] );
% arc_list = arc_list(idx);

end

function index = findArc(arc, arc_list)

index = 0;

for ii = 1:length(arc_list)
    p = 0;
    if [arc_list(ii).p.x arc_list(ii).p.y] == [arc.p.x arc.p.y]
        p = 1;
    end
    
    seg0 = 0;
    a.seg0 = arc.seg0;
    b.seg0 = arc_list(ii).seg0;
    
    if isempty(a.seg0) && isempty(b.seg0) 
        seg0 = 1;
    end
    
    if ~isempty(a.seg0) && ~isempty(b.seg0) 
        start_p1 = a.seg0.start_p;
        start_p2 = b.seg0.start_p;
        
        if [start_p1.x start_p1.y]==[start_p2.x start_p2.y]
            seg0 = 1;
        end
    end
    
    
    seg1 = 0;
    a.seg1 = arc.seg1;
    b.seg1 = arc_list(ii).seg1;
    if (isempty(a.seg1) && isempty(b.seg1))
        seg1 = 1;
    end
    
    if ~isempty(a.seg1) && ~isempty(b.seg1) 
        start_p1 = a.seg1.start_p;
        start_p2 = b.seg1.start_p;
        
        if [start_p1.x start_p1.y]==[start_p2.x start_p2.y]
            seg1 = 1;
        end
    end
    
    if p && seg0 && seg1
        index = ii;
        break;
    end
end

end

function [arc_list circle_events seg_list] = removeArc(arc, arc_list_input, circle_events_input, seg_list_input, v)

arc_list = arc_list_input;
circle_events = circle_events_input;
seg_list = seg_list_input;

index=findArc(arc, arc_list);

if index==0
    return
end

prev_arc = struct([]);
if index-1>0
%     update the seg1 prev arc
    arc_list(index-1).seg1.start_p = v;
    prev_arc = arc_list(index-1);
end

post_arc = struct([]);
if index+1<=length(arc_list)
%     update the post arc
    arc_list(index+1).seg0.start_p = v;
    post_arc = arc_list(index+1);
end

% update the arc
arc_list(index).seg0.end_p = v;
arc_list(index).seg1.end_p = v;

seg_list = pushSeg(arc_list(index).seg0, seg_list);
seg_list = pushSeg(arc_list(index).seg1, seg_list);



c = arc_list(index).circle_event;
if ~isempty(c)
    circle_events = removeCircleEvent(c, circle_events);
end

arc_list(index) = [];

if ~isempty(prev_arc)
    [arc_list, circle_events] = updateCircleEvents(prev_arc, arc_list, circle_events);
end

if ~isempty(post_arc)
    [arc_list, circle_events] = updateCircleEvents(post_arc, arc_list, circle_events);
end

end


function circle_events = addCircleEvent(c, circle_events_input)

circle_events = circle_events_input;

if ~isempty(circle_events)
    circle_events(length(circle_events)+1) = c;
    
    % sort the circle events by y- descending
    [~, idx] = sort( [circle_events.y], 'descend');
    circle_events = circle_events(idx);
else
    circle_events = c;
end

end

function index = findCircleEvent(circle_event, circle_events)
index = 0;

center = circle_event.center; 
radius = circle_event.radius;

for ii = 1:length(circle_events)
    c = circle_events(ii);
    
    if [center.x center.y radius] == [c.center.x c.center.y c.radius]
        index = ii;
        break;
    end
end

end

function circle_events = removeCircleEvent(circle_event, circle_events_input)

circle_events = circle_events_input;

index = findCircleEvent(circle_event, circle_events);
if index>0
    circle_events(index) = [];
end

end

function seg_list = finishArcList(arc_list, seg_list_input, axis_scaling)

seg_list = seg_list_input;

y = axis_scaling.ymin - (axis_scaling.xmax-axis_scaling.xmin) - (axis_scaling.ymax-axis_scaling.ymin);

while ~isempty(arc_list)
    if ~isempty(arc_list(1).seg1)
        node = intersection(arc_list(1).p, arc_list(2).p, y);
        
        arc_list(1).seg1.end_p = node;
        
        seg_list = pushSeg(arc_list(1).seg1, seg_list);    
    end
    
    % remove the arc
    arc_list(1) = [];
end

end

function seg_list = pushSeg(seg, seg_list_input)

seg_list = seg_list_input;

% store the seg in seg list
if isempty(seg_list)
    seg_list = seg;
else
    seg_list(length(seg_list)+1) = seg;
end

end


function showFigure(site_points, axis_scaling, y, arc_list, circle_event, seg_list)

% one new figure
figure;

hold on;

% the text axis offset by y
text_axis_offset = 0.5;

% plot the site points
for ii = 1:length(site_points)
    plot(site_points(ii).x, site_points(ii).y, 'bx');
    text(site_points(ii).x, site_points(ii).y + text_axis_offset, int2str(ii));
end

% plot the sweep line
x = linspace(axis_scaling.xmin, axis_scaling.xmax, 1000);
if ~isempty(y)
    plot(x, y);
end

% plot the arc
for ii = 1 : length(arc_list)
    if arc_list(ii).p.y~= y
%         plot the intersection point
        if ~isempty(arc_list(ii).seg0)
            plot(arc_list(ii).seg0.start_p.x, arc_list(ii).seg0.start_p.y, 'ro');
        end
        
        if ~isempty(arc_list(ii).seg1)
            plot(arc_list(ii).seg1.start_p.x, arc_list(ii).seg1.start_p.y, 'ro');
        end
        
        drawParabola(arc_list(ii).p, y, axis_scaling);
    end
end

% plot the circle
if ~isempty(circle_event)
    
    radius = circle_event.radius;
    center = circle_event.center;
    
    drawCircle(center, radius);
end


% plot the segments
for ii = 1 : length(seg_list)
%     plot the start & end point
    plot(seg_list(ii).start_p.x, seg_list(ii).start_p.y, 'ro');
    plot(seg_list(ii).end_p.x, seg_list(ii).end_p.y, 'r^');
    
    drawLine(seg_list(ii).start_p, seg_list(ii).end_p);
end

% set the x- and y- axes
axis([axis_scaling.xmin axis_scaling.xmax axis_scaling.ymin axis_scaling.ymax]);

hold off;

end