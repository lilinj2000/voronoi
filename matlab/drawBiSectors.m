function halfPlaneAlgorithm(p, r)
% HALFPLANEALGORITHM demonstrate the half-plane algorithm for Voronoi
% diagram
% Function HALFPLANEALGORITHM demonstrate the half-plane algorithm on
% voronoi diagram which site points is p, and the 2-D x- and y- axes is r
% 
% dd    -- demonstrate site point

% input which points should be visualize on the figure
while true
    dd = input(['Please input which site point should be demonstrate [1 - ' num2str(size(p,1)) ']:\n']) ;
    if dd<1 | dd>size(p,1)
        disp('the input invalid !!!');
    else
        break;
    end
end
    
% open a new figure
figure
hold on

% go through all the site points
for ii = 1:size(p,1)
%    draw the site points
    plot(p(ii, 1), p(ii, 2), 'bx');
    
% check if the point equa the demonstrate point    
    if ii~=dd
        x = linspace(r(1), r(2), 1000);
        % x1 + x2
        A = p(dd,1) + p(ii,1); 
        % y1 + y2
        B = p(dd,2) + p(ii,2);
        % x2 - x1
        C = p(ii,1) - p(dd,1);
        % y2 - y1
        D = p(ii,2) - p(dd,2);

%         on the same x- axes
        if C==0 
            y = (p(dd,2)+p(ii,2))./2;
        elseif D==0
%         on the same y- axes
            y = linspace(r(3), r(4), 1000);
            x = (p(dd,1)+p(ii,1))./2;
        else
%         not on the same x- or y- axes
            K = D ./ C;
            y = -1./K.*(x-A./2) + B./2;
        end
        
        plot(x, y, 'k-');
    end
end

axis(r);

hold off

end