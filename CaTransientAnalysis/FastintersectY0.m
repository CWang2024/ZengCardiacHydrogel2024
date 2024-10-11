function [xout, yout] = FastintersectY0(y0, x, y)

    % INTERSECT_WITH_Y0 Finds intersections of a curve with a horizontal line y = y0
    %   [XOUT, YOUT] = INTERSECT_WITH_Y0(Y0, X, Y) computes the points of 
    %   intersection between the curve defined by X and Y, and the horizontal line y = Y0.
    %
    %   X and Y are vectors where each pair (X(i), Y(i)) defines a point on the curve.
    %   XOUT and YOUT are the x and y coordinates of the intersection points.

    % Initialize output vectors
    xout = [];
    yout = [];

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    %startInd=find(y>y0, 1, "first")-1;
    %EndInd=find(y<y0, 1, "last")+1;
    diff=y-y0;
    % Loop through each segment of the curve
    for i = 1:length(y)-1
        % Check if the curve crosses y = y0 between points i and i+1
        if (diff(i)) * (diff(i+1)) < 0
            % Linear interpolation to find the exact intersection
            x1 = x(i);
            x2 = x(i + 1);
            y1 = y(i);
            y2 = y(i + 1);

            % Compute the x-coordinate of the intersection
            x_intersect = x1 + (y0 - y1) * (x2 - x1) / (y2 - y1);
            y_intersect = y0;

            % Append to output
            xout(end+1, 1) = x_intersect;
            yout(end+1, 1) = y_intersect;
        elseif y(i) == y0
            % If a point is exactly on the line y = y0
            xout(end+1, 1) = x(i);
            yout(end+1, 1) = y(i);
        end
    end
    
    % Check the last point if it's exactly on the line y = y0
    if y(end) == y0
        xout(end+1, 1) = x(end);
        yout(end+1, 1) = y(end);
    end

    xout=median(xout);
    yout=y0;
end
