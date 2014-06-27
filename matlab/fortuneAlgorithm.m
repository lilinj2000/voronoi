function fortuneAlgorithm(site_points, axis_scaling)
% FORTUNEALGORITHM demonstrate the sweeping line algorithm for Voronoi
% diagram
% Function FORTUNEALGORITHM which implement the sweeping line algorithm on
% voronoi diagram and demo it on each step.
% 
% Define variables:
% site_points        -- site points
% axis_scaling       -- the scaling for the x- and y-axes    
% Record Of revisions:
% Date Programmer Description of change
% ======== ============== ========================
% 6/16/2014 Linjiang Li Original code
% 6/17/2014 LinJIang Li handleSiteEvent
% 6/19/2014 LinJiang Li update the var name
% 6/26/2014 LinJiang Li fix the site event & circle event process
% 6/26/2014 LinJiang Li full implement the algorithm
% 6/27/2014 LinJiang Li fix the remove circle event bug

% Internal variables:
% site_points_event  -- site point events by y- axes descending
% circle_events      -- circle events by y- axes descending 
% arc_list           -- arc list, which store all the beach lines
% seg_list           -- segement list

% the internal struct is:
% site_points { p(x,y) }
% site_points_event { p(x,y) }
% arc_list { p, seg0(start_p, end_p), seg1(start_p, end_p), circle_event }
% circle_events { y, center(x,y), radius, arc }
% seg_list { start_p, end_p, p}
% start_p/end_p { x, y }

% main steps:
% 1. gen site events
% 2. handle site events or handle circle events
% 3. handle the left circle events
% 4. finish the left arc list

% sort site points by descending order on y axes
[~, idx] = sort([site_points.y], 'descend');
% gen site point events
site_point_events = site_points(idx);

% init circle event
circle_events = struct([]);

% init arc tree
arc_list = struct([]);

% init seg list
seg_list = struct([]);

% go through all the site events
while ~isempty(site_point_events)
    % go through all the circle events
    while ~isempty(circle_events)
       if  circle_events(1).y>site_point_events(1).y
           % on circle event
           [arc_list, circle_events, seg_list] = handleCircleEvent(circle_events(1), site_points, axis_scaling, arc_list, circle_events, seg_list);
       else
           break;
       end
    end
    
    % on site event
    [arc_list, circle_events] = handleSiteEvent(site_point_events(1), site_points, axis_scaling, arc_list, circle_events, seg_list);
    
    % remove the site point event
    site_point_events(1) = [];
end

% handle the left circle events
while ~isempty(circle_events)
   % handle circle event
   [arc_list, circle_events, seg_list] = handleCircleEvent(circle_events(1), site_points, axis_scaling, arc_list, circle_events, seg_list);    
end

% finish all the left arc list
seg_list = finishArcList(arc_list, seg_list, axis_scaling);

% last figure show
showFigure(site_points, axis_scaling, [], arc_list, [], seg_list);

end

% process site event 
% 1. init the output parameters
% 2. addArc to arc list
% 3. check & update circle events
function [arc_list, circle_events] = handleSiteEvent(p, site_points, axis_scaling, arc_list_input, circle_events_input, seg_list)

% init the output parameters
arc_list = arc_list_input;
circle_events = circle_events_input;

% add arc
[arc_list, circle_events] = addArc(p, arc_list, circle_events, axis_scaling.ymax);

% check & update circle events
% [arc_list, circle_events] = updateCircleEvents(arc, arc_list, circle_events);

showFigure(site_points, axis_scaling, p.y, arc_list, [], seg_list);

waitforbuttonpress;

end

% on circle event
% 1. init output parameters
% 2. remove Arc
function [arc_list, circle_events, seg_list] = handleCircleEvent(c, site_points, axis_scaling, arc_list_input, circle_events_input, seg_list_input)

% init output parameters
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

function [arc_list, circle_events] = updateCircleEvent(index, arc_list_input, circle_events_input)

% init the output parameters
arc_list = arc_list_input;
circle_events = circle_events_input;


c = [];
if index-1>0 && index+1<=length(arc_list)
    c = checkCircle(arc_list(index-1).p, arc_list(index).p, arc_list(index+1).p);
end
        
       
if ~isempty(c)
    c.arc = arc_list(index);
    arc_list(index).circle_event = c;
    
    circle_events = addCircleEvent(c, circle_events);
end
        

% check arc index
% index = findArc(arc, arc_list);
% if index>0
% %     do update circle events
% %  update index-2, index-1, index
% 
%     for ii = index-2:index
%         if ii>0 && ii+2<=length(arc_list)
%     %         remove circle event from arc_list & circle_events
%             c = arc_list(ii+1).circle_event;
%             
%             if ~isempty(c)
%                 [circle_events] = removeCircleEvent(c, circle_events);
%             end
%             arc_list(ii+1).circle_event = [];
%     
%     %         has the circle events
%             c = checkCircle(arc_list(ii).p, arc_list(ii+1).p, arc_list(ii+2).p);
%             if ~isempty(c)
%                 c.arc = arc_list(ii+1);
%                 arc_list(ii+1).circle_event = c;
%                 % add the full arc info for the circle_event
%                 c.arc = arc_list(ii+1);
%                 
%                 circle_events = addCircleEvent(c, circle_events);
%             end
%         end
%     end
% end

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

function circle_event = checkCircle(pii, pjj, pkk)

% set the circle_event is empty
circle_event = struct([]);

% ensure edge pii-pjj on the left pii-pkk
if (pkk.y-pii.y)*(pjj.x-pii.x) - (pjj.y-pii.y)*(pkk.x-pii.x)>0
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

if G*D == H*C
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
% c.p = pjj;

circle_event = c;

end


% on site event, one new arc should be added.
% example: p1 - p3 - p4 - p2 - p1
% on p5 event: ===> p1 - p3 - (p4 - p5 - p4) - p2 - p1
% main steps:
% 1. init outpu parameters
% 2. on the first arc, just add it, and return.
% 3. check intersect, and insert the right position
% 4. remove the false alarm for circle event
% 5. update prev.seg1 & arc.seg0 & arc.seg1 & post.seg0
% 6. update circle events for prev & post
function [arc_list, circle_events] = addArc(p, arc_list_input, circle_events_input, ymax)

% init output parameters
arc_list = arc_list_input;
circle_events = circle_events_input;

% new arc
arc.p = p;
arc.circle_event = [];
arc.seg0 = [];
arc.seg1 = [];

% check whether the fist arc
if isempty(arc_list)
    arc_list = arc;
    return ;
end

% check insection and insert the arc
for ii=1:length(arc_list)
    node = intersect(p, ii, arc_list);
    
    if ~isempty(node)
        
        % remove this arc circle event, for it's false alarm
        % circle events update must before on the seg change.
        if ~isempty(arc_list(ii).circle_event)
            c = arc_list(ii).circle_event ;
            arc_list(ii).circle_event = [];
            circle_events = removeCircleEvent(c, circle_events);
        end

        
        len = length(arc_list);
        arc_list(len+2:-1:ii+2) = arc_list(len:-1:ii);
        
        %update the arc seg
        arc.seg0.start_p = node;
        arc.seg1.start_p = node;
        
        arc_list(ii+1) = arc;
        
        %update the pre arc seg1
        arc_list(ii).seg1.start_p = node;
        
        %update the post arc seg0
        arc_list(ii+2).seg0.start_p = node;
        
        % check the circle event, for it's false alarm
%         c = arc_list(ii).circle_event;
%         if ~isempty(c)
%             circle_events = removeCircleEvent(c, circle_events);
%             
%             arc_list(ii).circle_event = [];
%             arc_list(ii+2).circle_event = [];
%         end
%         
        % update the circle event for prev, and post
        [arc_list, circle_events] = updateCircleEvent(ii, arc_list, circle_events);
        [arc_list, circle_events] = updateCircleEvent(ii+2, arc_list, circle_events);
        
        
        return ;
    end
end

% for the specail point, which no intersect
start_p.x = ymax;
start_p.y = (arc_list(ii).p.x+p.x)./2;

arc.seg0.start_p = start_p;

arc_list(ii).seg1.start_p = start_p;
arc_list(ii+1) = arc;

% check & update the prev arc, circle events.
[arc_list, circle_events] = updateCircleEvent(ii, arc_list, circle_events);

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

function ok = equalArc(arc1, arc2)
ok = 0;

if isempty(arc1) || isempty(arc2)
    return;
end

p_ok = 0;

p1 = arc1.p;
p2 = arc2.p;

if p1.x==p2.x && p1.y==p2.y
    p_ok = 1;
end

seg0_ok =0;

% both of the seg0 is empty
if isempty(arc1.seg0) && isempty(arc2.seg0) 
    seg0_ok = 1;
end

% both of the seg0 is not empty
if ~isempty(arc1.seg0) && ~isempty(arc2.seg0) 
    start_p1 = arc1.seg0.start_p;
    start_p2 = arc2.seg0.start_p;
    % the start point equal
    if start_p1.x==start_p2.x && start_p1.y==start_p2.y
        seg0_ok = 1;
    end
end
    
    
seg1_ok = 0;
if (isempty(arc1.seg1) && isempty(arc2.seg1))
    seg1_ok = 1;
end

if ~isempty(arc1.seg1) && ~isempty(arc2.seg1) 
    start_p1 = arc1.seg1.start_p;
    start_p2 = arc2.seg1.start_p;

    if start_p1.x==start_p2.x && start_p1.y==start_p2.y
        seg1_ok = 1;
    end
end


if p_ok && seg0_ok && seg1_ok
    ok = 1;
end

end

function index = findArc(arc, arc_list)

index = 0;

for ii = 1:length(arc_list)
    
    if equalArc(arc, arc_list(ii))
        index = ii;
        break;
    end
end

end

% on circle event, should do remove arc
% 1. init output parameters.
% 2. find the right arc in list
% 3. remove the circle event of the arc
% 4. update the prev seg1 & post arc seg0
% 5. finish the arc seg & push to seg list
% 6. remove the arc
% 7. check & update circle events

function [arc_list circle_events seg_list] = removeArc(arc, arc_list_input, circle_events_input, seg_list_input, v)

% init output parameters
arc_list = arc_list_input;
circle_events = circle_events_input;
seg_list = seg_list_input;

index=findArc(arc, arc_list);
if index==0
    disp('can''t find the arc on remove arc.');
    return
end

% attentions: the false alarm for the circle events must be removed before
% the seg change
if index-1>0
    % remove the false alarm of the circle events for the prev
    c = arc_list(index-1).circle_event;
    if ~isempty(c)
        circle_events = removeCircleEvent(c, circle_events);
    end
    
    % update the seg1 prev arc
    arc_list(index-1).seg1.start_p = v;
end

if index+1<=length(arc_list)
    % remove the false alarm of the circle events for the prev
    c = arc_list(index+1).circle_event;
    if ~isempty(c)
        circle_events = removeCircleEvent(c, circle_events);
    end
    
    % update the post arc
    arc_list(index+1).seg0.start_p = v;
end

% remove the current circle events from list 
c = arc_list(index).circle_event;
if ~isempty(c)
    circle_events = removeCircleEvent(c, circle_events);
end

% update the arc
arc_list(index).seg0.end_p = v;
arc_list(index).seg1.end_p = v;

seg_list = pushSeg(arc_list(index).seg0, seg_list);
seg_list = pushSeg(arc_list(index).seg1, seg_list);

% remove the arc from list
arc_list(index) = [];

% the index point to the post arc now
% check & update the prev & post arc circle event
[arc_list, circle_events] = updateCircleEvent(index-1, arc_list, circle_events);
[arc_list, circle_events] = updateCircleEvent(index, arc_list, circle_events);

end


function circle_events = addCircleEvent(c, circle_events_input)

circle_events = circle_events_input;

if isempty(c)
    return ;
end

if ~isempty(circle_events)
    circle_events(length(circle_events)+1) = c;
    
    % sort the circle events by y- descending
    [~, idx] = sort( [circle_events.y], 'descend');
    circle_events = circle_events(idx);
else
    circle_events = c;
end

end

function ok = equalCircleEvent(c1, c2)
ok = 0;

% any of the circle event is empty
if isempty(c1) || isempty(c2)
    return ;
end

if equalArc(c1.arc, c2.arc)
    ok = 1;
end

end

function index = findCircleEvent(circle_event, circle_events)
index = 0;

if isempty(circle_event)
    return ;
end

for ii = 1:length(circle_events)
    if equalCircleEvent(circle_event, circle_events(ii))
        index = ii;
        break;
    end
end

% center = circle_event.center; 
% radius = circle_event.radius;
% 
% for ii = 1:length(circle_events)
%     c = circle_events(ii);
%     
%     if center.x==c.center.x && center.y==c.center.y && radius==c.radius
% %      if [center.x center.y radius]==[c.center.x c.center.y c.radius]
%         index = ii;
%         break;
%     end
% end

end

function circle_events = removeCircleEvent(circle_event, circle_events_input)

circle_events = circle_events_input;

index = findCircleEvent(circle_event, circle_events);
if index>0
    circle_events(index) = [];
else
    disp('error on remove circle event');
end

end

function seg_list = finishArcList(arc_list, seg_list_input, axis_scaling)

seg_list = seg_list_input;

y = axis_scaling.ymin - (axis_scaling.xmax-axis_scaling.xmin).*2 - (axis_scaling.ymax-axis_scaling.ymin).*2;

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