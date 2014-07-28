classdef VoronoiFortuneAlgo < handle
    %VORONOIFORTUNEALGO one voronoi diagram by fortune algorithm
    
    
    properties
        %site_points = struct([]);        % all of the site points
        
        site_points;
        axis_scaling = [];
        
        site_point_events = struct([]);  % the site point event queue by y- descend
        circle_events = struct([]);      % the circle event queue by y- descend
        seg_list = struct([]);           % all of the edge list
        arc_list = struct([]);           % all the beach lines 
        
        edge_list = cell(0);            % store the site points => seg list
        
    end
    
    methods
        
        % site_points which is m*2 array
        % axis_opt:
        function VFA = VoronoiFortuneAlgo(site_points, axis_ratio)
            
            [site_points_sort, I, ~] = unique(site_points, 'rows', 'first');
            [~, I] = sort(I);
            site_points = site_points_sort(I, :);
            
            if size(site_points, 1)<2 || size(site_points, 2)~=2
                disp('invalid input for site points, which must be m*2 array');
                return;
            end
            
            % init the array
            VFA.site_points(size(site_points,1)).x = [];
            VFA.site_points(size(site_points,1)).y = [];
            for ii = 1:size(site_points,1)
                VFA.site_points(ii).x = site_points(ii, 1);
                VFA.site_points(ii).y = site_points(ii, 2);
            end

            % calculate [XMIN XMAX YMIN YMAX]
            xmin = min([VFA.site_points.x]);
            xmax = max([VFA.site_points.x]);
            ymin = min([VFA.site_points.y]);
            ymax = max([VFA.site_points.y]);
            
            if ~isnumeric(axis_ratio)
                disp('input invalid, the axis_ratio must be numeric.');
                return;
            end
            
            delta = max(xmax, ymax) - min(xmin, ymin);
            
            VFA.axis_scaling.ymax = max(xmax, ymax) + delta.*axis_ratio + 1;
            VFA.axis_scaling.xmax = VFA.axis_scaling.ymax;
            
            VFA.axis_scaling.ymin = min(xmin, ymin) - delta.*axis_ratio - 1;
            VFA.axis_scaling.xmin = VFA.axis_scaling.ymin;
                        
        end
        
        function do(VFA)
            % sort site points by descending order on y axes
            [~, idx] = sort([VFA.site_points.y], 'descend');
            % gen site point events
            VFA.site_point_events = VFA.site_points(idx);
            
            
            % go through all the site events
            while ~isempty(VFA.site_point_events)
                % go through all the circle events
                while ~isempty(VFA.circle_events)
                    if  VFA.circle_events(1).y>VFA.site_point_events(1).y
                       
                       %VFA.showFigure(VFA.circle_events(1).y, VFA.circle_events(1));
                       %waitforbuttonpress;
                       
                       % on circle event
                       VFA.onCircleEvent(VFA.circle_events(1));
                    else
                       break;
                    end
                end
    
                % on site event
                VFA.onSiteEvent(VFA.site_point_events(1));
                
                %VFA.showFigure(VFA.site_point_events(1).y, []);
                %waitforbuttonpress;
                
                % remove the site point event
                VFA.site_point_events(1) = [];
            end
            
            % handle the left circle events
            while ~isempty(VFA.circle_events)
                
                %VFA.showFigure(VFA.circle_events(1).y, VFA.circle_events(1));
                %waitforbuttonpress;
               
                % handle circle event
               VFA.onCircleEvent(VFA.circle_events(1));    
            end

            % finish all the left arc list
            VFA.finishLeftArcList();
            
            VFA.showFigure([], []);
            %waitforbuttonpress;
            
        end
        
        function onSiteEvent(VFA, p)
            % new arc
            arc.p = p;
            arc.circle_event = [];
            arc.seg0 = [];
            arc.seg1 = [];

            % check whether the fist arc
            if isempty(VFA.arc_list)
                VFA.arc_list = arc;
                return ;
            end

            % check insection and insert the arc
            for ii=1:length(VFA.arc_list)
                node = VFA.intersect(p, ii);

                if ~isempty(node)
                    % remove this arc circle event, for it's false alarm
                    % circle events update must before on the seg change.
                    if ~isempty(VFA.arc_list(ii).circle_event)
                        c = VFA.arc_list(ii).circle_event ;
                        VFA.arc_list(ii).circle_event = [];
                        VFA.removeCircleEvent(c);
                    end


                    len = length(VFA.arc_list);
                    VFA.arc_list(len+2:-1:ii+2) = VFA.arc_list(len:-1:ii);

                    %update the arc seg
                    arc.seg0.start_p = node;
                    arc.seg1.start_p = node;

                    VFA.arc_list(ii+1) = arc;

                    %update the pre arc seg1
                    VFA.arc_list(ii).seg1.start_p = node;

                    %update the post arc seg0
                    VFA.arc_list(ii+2).seg0.start_p = node;

              
                    % update the circle event for prev, and post
                    VFA.updateCircleEvent(ii);
                    VFA.updateCircleEvent(ii+2);

                    return ;
                end
            end
            
            % for the specail point, which no intersect
            start_p.x = ymax;
            start_p.y = (VFA.arc_list(ii).p.x+p.x)./2;

            arc.seg0.start_p = start_p;

            VFA.arc_list(ii).seg1.start_p = start_p;
            VFA.arc_list(ii+1) = arc;

            % check & update the prev arc, circle events.
            VFA.updateCircleEvent(ii);
        end
        
        function onCircleEvent(VFA, c)
            
            index=VFA.findArc(c.arc);
            if index==0
                disp('can''t find the arc on circle event.');
                return
            end

            % attentions: the false alarm for the circle events must be removed before
            % the seg change
            if index-1>0
                % remove the false alarm of the circle events for the prev
                c1 = VFA.arc_list(index-1).circle_event;
                if ~isempty(c1)
                    VFA.removeCircleEvent(c1);
                end

                % update the seg1 prev arc
                VFA.arc_list(index-1).seg1.start_p = c.center;
            end

            if index+1<=length(VFA.arc_list)
                % remove the false alarm of the circle events for the prev
                c2 = VFA.arc_list(index+1).circle_event;
                if ~isempty(c2)
                    VFA.removeCircleEvent(c2);
                end

                % update the post arc
                VFA.arc_list(index+1).seg0.start_p = c.center;
            end

            % remove the current circle events from list 
            VFA.removeCircleEvent(c);

            % update the arc
            VFA.arc_list(index).seg0.end_p = c.center;
            VFA.arc_list(index).seg1.end_p = c.center;

            VFA.pushSeg(VFA.arc_list(index).seg0, VFA.arc_list(index-1).p, VFA.arc_list(index).p);
            VFA.pushSeg(VFA.arc_list(index).seg1, VFA.arc_list(index).p, VFA.arc_list(index+1).p);

            % remove the arc from list
            VFA.arc_list(index) = [];

            % the index point to the post arc now
            % check & update the prev & post arc circle event
            VFA.updateCircleEvent(index-1);
            VFA.updateCircleEvent(index);
        end
        
        function updateCircleEvent(VFA, index)

            c = [];
            if index-1>0 && index+1<=length(VFA.arc_list)
                c = VoronoiFortuneAlgo.checkCircle(VFA.arc_list(index-1).p, ...
                    VFA.arc_list(index).p, VFA.arc_list(index+1).p);
            end


            if ~isempty(c)
                c.arc = VFA.arc_list(index);
                VFA.arc_list(index).circle_event = c;

                VFA.addCircleEvent(c);
            end
        end
        
        function node = intersect(VFA, p, index)

            node = struct([]);

            if p.y==VFA.arc_list(index).p.y
                % it's not intersect
                return ;
            end

            node0 = struct([]);
            node1 = struct([]);

            if index-1>0
                node0 = VoronoiFortuneAlgo.intersection(VFA.arc_list(index-1).p, VFA.arc_list(index).p, p.y);
            end

            if index+1<=length(VFA.arc_list)
                node1 = VoronoiFortuneAlgo.intersection(VFA.arc_list(index).p, VFA.arc_list(index+1).p, p.y);
            end

            % between node0 and node1, then intersect with this arc
            if (isempty(node0) || node0.x<=p.x) && (isempty(node1) || node1.x>=p.x)
                node = VoronoiFortuneAlgo.intersection(VFA.arc_list(index).p, p, p.y);
            end
        end
        
        function index = findArc(VFA, arc)
            index = 0;

            for ii = 1:length(VFA.arc_list)

                if VoronoiFortuneAlgo.equalArc(arc, VFA.arc_list(ii))
                    index = ii;
                    break;
                end
            end
        end
        
        function addCircleEvent(VFA, c)

            if isempty(c)
                return ;
            end

            if ~isempty(VFA.circle_events)
                VFA.circle_events(length(VFA.circle_events)+1) = c;

                % sort the circle events by y- descending
                [~, idx] = sort( [VFA.circle_events.y], 'descend');
                VFA.circle_events = VFA.circle_events(idx);
            else
                VFA.circle_events = c;
            end
        end
        
        function index = findCircleEvent(VFA, c)
            index = 0;

            if isempty(c)
                return ;
            end

            for ii = 1:length(VFA.circle_events)
                if VoronoiFortuneAlgo.equalCircleEvent(c, VFA.circle_events(ii))
                    index = ii;
                    break;
                end
            end
        end
        
        function removeCircleEvent(VFA, c)

            index = VFA.findCircleEvent(c);
            if index>0
                VFA.circle_events(index) = [];
            else
                disp('error on remove circle event');
            end
        end
        
        function finishLeftArcList(VFA)

            y = VFA.axis_scaling.ymin - ...
                (VFA.axis_scaling.xmax-VFA.axis_scaling.xmin).*2 ...
                - (VFA.axis_scaling.ymax-VFA.axis_scaling.ymin).*2;
            
            while ~isempty(VFA.arc_list)
                if ~isempty(VFA.arc_list(1).seg1)
                    node = VoronoiFortuneAlgo.intersection(VFA.arc_list(1).p, VFA.arc_list(2).p, y);

                    VFA.arc_list(1).seg1.end_p = node;

                    VFA.pushSeg(VFA.arc_list(1).seg1, VFA.arc_list(1).p, VFA.arc_list(2).p);    
                end

                % remove the arc
                VFA.arc_list(1) = [];
            end
        end

        function pushSeg(VFA, seg, p1, p2)

            % store the seg in seg list
            if isempty(VFA.seg_list)
                VFA.seg_list = seg;
            else
                VFA.seg_list(length(VFA.seg_list)+1) = seg;
            end
            
            VFA.pushEdge(seg, p1);
            VFA.pushEdge(seg, p2);
                        
        end
        
        function pushEdge(VFA, seg, p)
            
            index = VFA.findEdgeList(p);
            if index~=0
                VFA.edge_list{index, 2} = [VFA.edge_list{index, 2}; seg];
            else
                index = size(VFA.edge_list, 1)+1;
                
                VFA.edge_list{index, 1} = p;
                VFA.edge_list{index, 2} = seg;
            end
        end
        
        function index = findEdgeList(VFA, p)
            index = 0;
            
            if isempty(VFA.edge_list)
                return ;
            end
            
            for ii = 1 : size(VFA.edge_list, 1)
                if VFA.edge_list{ii, 1}.x==p.x && VFA.edge_list{ii, 1}.y==p.y
                    index = ii;
                    return ;
                end
            end
        end
        
        % way: 1 - max value; 2 - mean value
        function radius = calculateRadius(VFA, p, way)
            max_radius = 30000; % 30km
            radius = 0;
            
            xmin = VFA.axis_scaling.xmin;
            xmax = VFA.axis_scaling.xmax;
            ymin = VFA.axis_scaling.ymin;
            ymax = VFA.axis_scaling.ymax;
            
            index = VFA.findEdgeList(p);

            if index~=0
                dist = [];
                for ii = 1 : size(VFA.edge_list{index, 2}, 1)
                    points = [VFA.edge_list{index,2}(ii,1).start_p; ...
                            VFA.edge_list{index,2}(ii,1).end_p ];
                    
                    
                    for jj  = 1: length(points)
                        distance = VFA.distance(p, points(jj));
                        if distance>max_radius
                            distance = max_radius;
                        end

                        dist = [dist; distance];
                    end
                                        
                end

                if way==1 % max value
                    radius = max(dist);
                elseif way==2 % mean value
                    radius = mean(dist);
                end
                
                
            end

        end
        
        function showFigure(VFA, y, c)

            % one new figure
            figure;

            hold on;

            % the text axis offset by y
            text_axis_offset = 0.5;

            % plot the site points
            for ii = 1:length(VFA.site_points)
                plot(VFA.site_points(ii).x, VFA.site_points(ii).y, 'bx');
                text(VFA.site_points(ii).x, VFA.site_points(ii).y + text_axis_offset, int2str(ii));
            end

            % plot the sweep line
            x = linspace(VFA.axis_scaling.xmin, VFA.axis_scaling.xmax, 1000);
            if ~isempty(y)
                plot(x, y);
            end

            % plot the arc
            for ii = 1 : length(VFA.arc_list)
                if VFA.arc_list(ii).p.y ~= y
                    if ~isempty(y) && isempty(c) % just site event
                        % plot the intersection point
                        if ~isempty(VFA.arc_list(ii).seg0)
                            plot(VFA.arc_list(ii).seg0.start_p.x, VFA.arc_list(ii).seg0.start_p.y, 'ro');
                        end


                        if ~isempty(VFA.arc_list(ii).seg1)
                            plot(VFA.arc_list(ii).seg1.start_p.x, VFA.arc_list(ii).seg1.start_p.y, 'ro');
                        end
                    end

                    if ii-1>0
                        node = VoronoiFortuneAlgo.intersection(VFA.arc_list(ii-1).p, VFA.arc_list(ii).p, y);
                        xmin = node.x;
                        
                        VoronoiFortuneAlgo.drawLine(node, VFA.arc_list(ii).seg0.start_p, 'k-.');
                    else
                        xmin = VFA.axis_scaling.xmin;
                    end
                    
                    if ii+1<=length(VFA.arc_list)
                        node = VoronoiFortuneAlgo.intersection(VFA.arc_list(ii).p, VFA.arc_list(ii+1).p, y);
                        xmax = node.x;
                        
                        VoronoiFortuneAlgo.drawLine(VFA.arc_list(ii).seg1.start_p, node, 'k-.');
                    else
                        xmax = VFA.axis_scaling.xmax;
                    end
                    
                    VoronoiFortuneAlgo.drawParabola(VFA.arc_list(ii).p, y, xmin, xmax);
                end
            end
            
            % plot the circle arc
            if ~isempty(c) % this is circle event
                arc = c.arc;
                
                VoronoiFortuneAlgo.drawParabola(arc.p, y, ...
                    VFA.axis_scaling.xmin, VFA.axis_scaling.xmax);
            end

            % plot the circle
            if ~isempty(c)

                radius = c.radius;
                center = c.center;

                VoronoiFortuneAlgo.drawCircle(center, radius);
            end


            % plot the segments
            for ii = 1 : length(VFA.seg_list)
                % plot the start & end point
                % plot(VFA.seg_list(ii).start_p.x, VFA.seg_list(ii).start_p.y, 'ro');
                plot(VFA.seg_list(ii).end_p.x, VFA.seg_list(ii).end_p.y, 'r^');

                VoronoiFortuneAlgo.drawLine(VFA.seg_list(ii).start_p, VFA.seg_list(ii).end_p, 'r-');
            end

            % set the x- and y- axes
            axis([VFA.axis_scaling.xmin VFA.axis_scaling.xmax VFA.axis_scaling.ymin VFA.axis_scaling.ymax]);

            hold off;
        end
        
    end
    
    methods( Static )
        
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
        
        function c = checkCircle(pii, pjj, pkk)
            % set the circle_event is empty
            c = struct([]);

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

            c(1).y = center.y-radius;
            c(1).center = center;
            c(1).radius = radius;
            % c.p = pjj;

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
        
        function ok = equalCircleEvent(c1, c2)
            ok = 0;

            % any of the circle event is empty
            if isempty(c1) || isempty(c2)
                return ;
            end

            if VoronoiFortuneAlgo.equalArc(c1.arc, c2.arc)
                ok = 1;
            end
        end
        
        function drawParabola(p, line, xmin, xmax)

            x = linspace(xmin, xmax, 1000);

            % (y-lineY).^2 = (x-pii(1,1)).^2 + (y-pii(1,2)).^2;

            y = ((x - p.x).^2 + p.y.^2 - line.^2)./(2.*(p.y-line));

            plot(x, y, 'k-');

            % plot(x, lineY);

        end
        
        function drawCircle(center, radius)
 
            angle = linspace(0, 2.*pi, 1000);
            
            plot(center.x + radius.*cos(angle), ...
                center.y + radius.*sin(angle), 'b-');
            
            % x = linspace(center.x- radius, center.x + radius, 1000);
            % y1 = sqrt(radius.^2-(x-center.x).^2) + center.y;
            % y2 = -sqrt(radius.^2-(x-center.x).^2) + center.y;
            % plot(x, y1, 'b-', x, y2, 'b-');

        end
        
        function drawLine(start_p, end_p, s)

            % y - y1 = (y2-y1)/(x2-x1)*(x-x1)

            if start_p.x~=end_p.x
                k = (end_p.y-start_p.y)./(end_p.x-start_p.x);

                x = linspace(start_p.x, end_p.x, 1000);
                y = start_p.y + k.*(x-start_p.x);
            else
                x = start_p.x;
                y = linspace(start_p.y, end_p.y, 1000);
            end

            plot(x, y, s);

        end
        
        function dist = distance(p1, p2)
            dist = sqrt((p1.x-p2.x).^2 + (p1.y-p2.y).^2);
        end
    end
end

