function halfPlaneAlgorithm(site_points, axis_scaling)
% HALFPLANEALGORITHM demonstrate the half-plane algorithm for Voronoi
% diagram
% Function HALFPLANEALGORITHM demonstrate the half-plane algorithm on
% voronoi diagram which site points is p, and the 2-D x- and y- axes is r
% 
% Define variables:
% site_points           -- site points
% axis_scaling          -- the scaling for the x- and y-axes
% demo_point_index      -- demonstrate site point index
% Record Of revisions:
% Date Programmer Description of change
% ======== ============== ========================
% 6/16/2014 Linjiang Li Original code
% 6/18/2014 Linjiang Li update the var name


% input which points should be visualize on the figure
while true
    demo_point_index = input(['Please input which site point should be demonstrate [1 - ' num2str(length(site_points)) ']:\n']) ;
    if demo_point_index<1 | demo_point_index>length(site_points)
        disp('the input invalid !!!');
    else
        break;
    end
end
    
% open a new figure
figure
hold on

demo_point = site_points(demo_point_index);

% the text axis offset by y
text_axis_offset = 0.5;

% go through all the site points
for ii = 1:length(site_points)
        
    % check if the point equa the demonstrate point    
    if ii~=demo_point_index
    
        %    draw the site points
        plot(site_points(ii).x, site_points(ii).y, 'bx');
        text(site_points(ii).x, site_points(ii).y + text_axis_offset, int2str(ii));
    
        x = linspace(axis_scaling.xmin, axis_scaling.xmax, 1000);
        
        p = site_points(ii);
        % x1 + x2
        A = demo_point.x + p.x; 
        % y1 + y2
        B = demo_point.y + p.y;
        % x2 - x1
        C = p.x - demo_point.x;
        % y2 - y1
        D = p.y - demo_point.y;

%         on the same x- axes
        if C==0 
            y = (demo_point.y + p.y)./2;
        elseif D==0
%         on the same y- axes
            y = linspace(axis_scaling.ymin, axis_scaling.ymax, 1000);
            x = (demo_point.x + p.x)./2;
        else
%         not on the same x- or y- axes
            K = D ./ C;
            y = -1./K.*(x-A./2) + B./2;
        end
        
        plot(x, y, 'k-');
    else
    %    draw the site points
        plot(site_points(ii).x, site_points(ii).y, 'bo');
        text(site_points(ii).x, site_points(ii).y + text_axis_offset, int2str(ii));    
    end
end

axis([axis_scaling.xmin axis_scaling.xmax axis_scaling.ymin axis_scaling.ymax]);

hold off

end