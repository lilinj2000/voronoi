% Script file:VoronoiShow.m
%
% Purpose:
% To visualize the algorithm of voronoi, including half-plane and fortune
%
% Record of revisions:
% Date Programmer Description of change
% ==== ========= ================
% 6/16/2014 Linjiang Li Original code
% 6/18/2014 Linjiang Li update the var name
%
% Define variables:
% site_points_input     -- input the site pionts
% site_points           -- site points for struct
% axis_scaling          -- the scaling for the x- and y-axes    
% selected_algorithm    -- half-plane or fortune algorithm

clear;
clc;

% just for debug info
debug = 1;

if debug
    selected_algorithm = 2;
    % two points
    site_points_input = [ 7 6; 6 7 ];
    
    % three points
%     site_points_input = [ 7 7; 5 6; 6 2 ];
    
    % more ponits
%     site_points_input = [ 7 7; 5 6; 6 2; 8 3];
    
    % more and more points
%     site_points_input = [ 7 7; 4 12; 5 10; 9 14; 10 9; 3 9; 6 20; 15 15; 12 13];
end

while true
    
    if ~debug
        selected_algorithm = input('Please select the demonstrate algorithm: \n 1 - half-plane algorithm \n 2 - fortune algorithm \n');
    end
    
    if selected_algorithm~=1 && selected_algorithm~=2
        disp('invalid input !!!');
    else
        break;
    end
end

% input site points
while true
    if ~debug
        site_points_input = input('Please input site points [x1, y2; x2, y2; ...]:\n');
    end
    
    % unique the site points
    site_opints_input = unique(site_points_input, 'rows');
    
    % check the input data
    if size(site_points_input, 1)<2 || size(site_points_input,2)~=2
        % site pionts less than 2
        % or the format is wrong (not the x, y axes)
        disp('the input invalid !!!');
    else
        break;
    end
end

% change the site points input for struct
site_points(size(site_points_input,1)).x = [];
site_points(size(site_points_input,1)).y = [];

for ii = 1:size(site_points_input,1)
    site_points(ii).x = site_points_input(ii, 1);
    site_points(ii).y = site_points_input(ii, 2);
end

% calculate [XMIN XMAX YMIN YMAX]
xmin = min([site_points.x]);
xmax = max([site_points.x]);
ymin = min([site_points.y]);
ymax = max([site_points.y]);

% set the expand ratio 0.25
ratio = 0.25;
% calculat the x- and y- axes delta
delta_x = xmax - xmin;
delta_y = ymax - ymin;
% updat the R by the ratio
% add 1 for prevent equal between min and max value
axis_scaling.xmin = xmin - delta_x.*ratio - 1; 
axis_scaling.xmax = xmax + delta_x.*ratio + 1;
axis_scaling.ymin = ymin - delta_y.*ratio - 1;
axis_scaling.ymax = ymax + delta_y.*ratio + 1;

if selected_algorithm==1
    % draw bisectors for the demonstrate site point
    halfPlaneAlgorithm(site_points, axis_scaling);
else
    fortuneAlgorithm(site_points, axis_scaling);
end


