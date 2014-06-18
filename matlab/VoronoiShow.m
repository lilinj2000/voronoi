% Script file:VoronoiShow.m
%
% Purpose:
% To visualize the algorithm of voronoi, including half-plane and fortune
%
% Record of revisions:
% Date Programmer Description of change
% ==== ========= ================
% 6/16/2014 Linjiang Li Original code
%
% Define variables:
% p     -- site points
% r     -- the scaling for the x- and y-axes    
% selectedAlgorithm -- half-plane or Fortune algorithm


while true
%     selectedAlgorithm = input('Please select the demonstrate algorithm: \n 1 - half-plane algorithm \n 2 - fortune algorithm \n');
    selectedAlgorithm = 2;
    if selectedAlgorithm~=1 & selectedAlgorithm~=2
        disp('invalid input !!!');
    else
        break;
    end
end

% input site points
while true
%     p = input('Please input site points [x1, y2; x2, y2; ...]:\n');

    p = [ 7 7; 5 6; 6 2; 8 3];
    
    % check the input data
    if size(p, 1)<2 | size(p,2)~=2
        % site pionts less than 2
        % or the format is wrong (not the x, y axes)
        disp('the input invalid !!!');
    else
        break;
    end
end

% init the r [XMIN XMAX YMIN YMAX]
r(1) = min(p(:, 1));
r(2) = max(p(:, 1));
r(3) = min(p(:, 2));
r(4) = max(p(:, 2));

% set the expand ratio 0.25
RATIO = 0.25;
% calculat the x- and y- axes delta
deltaX = r(2) - r(1);
deltaY = r(4) - r(3);
% updat the R by the ratio
% add 1 for prevent equal between min and max value
r(1) = r(1) - deltaX.*RATIO - 1; 
r(2) = r(2) + deltaX.*RATIO + 1;
r(3) = r(3) - deltaY.*RATIO - 1;
r(4) = r(4) + deltaY.*RATIO + 1;



if selectedAlgorithm==1
    % draw bisectors for the demonstrate site point
    halfPlaneAlgorithm(p, r);
else
    fortuneAlgorithm(p, r);
end


