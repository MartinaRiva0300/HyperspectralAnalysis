function [x_elbow, x_elbow_idx] = elbow_pt(x,y, idx)
% Find the maximum point
[~, maxIndex] = max(y);
maxPoint = [x(maxIndex), y(maxIndex)];
endPoint = [x(idx), y(idx)];

% Define the endpoint
if x(idx)>x(maxIndex)
% Calculate the slope and intercept of the line
    slope = (endPoint(2) - maxPoint(2)) / (endPoint(1) - maxPoint(1));
    intercept = maxPoint(2) - slope * maxPoint(1);
elseif x(idx)<x(maxIndex)
    slope = -(endPoint(2) - maxPoint(2)) / -(endPoint(1) - maxPoint(1));
    intercept = maxPoint(2) - slope * maxPoint(1);
end

% Calculate perpendicular distances
distances = zeros(size(x));
indices = sort([idx, maxIndex]);
for i = indices(1):indices(2)
    point = [x(i), y(i)];
    
    % Perpendicular slope
    perpSlope = -1 / slope;
    
    % Line equation of perpendicular line through the point
    perpIntercept = point(2) - perpSlope * point(1);
    
    % Find intersection of perpendicular line with the original line
    intersectX = (perpIntercept - intercept) / (slope - perpSlope);
    intersectY = slope * intersectX + intercept;
    
    % Calculate the distance between the point and the intersection
    distances(i) = sqrt((point(1) - intersectX)^2 + (point(2) - intersectY)^2);
end


[~,x_elbow_idx] = max(distances);
x_elbow = x(x_elbow_idx);

end
% 
% % Plot the curve and the line
% figure;
% plot(x, y, 'b-', 'DisplayName', 'Curve');
% hold on;
% plot([maxPoint(1), endPoint(1)], [maxPoint(2), endPoint(2)], 'r--', 'DisplayName', 'Line through Max and Endpoint');
% plot(maxPoint(1), maxPoint(2), 'go', 'DisplayName', 'Maximum Point');
% plot(endPoint(1), endPoint(2), 'mo', 'DisplayName', 'Endpoint');
% legend show;
% xlabel('x');
% ylabel('y');
% title('Curve and Line with Perpendicular Distances');
% grid on;
% plot(x(idx_max), y(idx_max), "ro");
% hold off;
% 
% % Plot the distances
% figure;
% plot(x, distances, 'k-', 'DisplayName', 'Perpendicular Distances');
% xlabel('x');
% ylabel('Distance');
% title('Perpendicular Distances to the Line');
% legend show;
% grid on;
% 
% end