function  varargout = DrawMoranScatterplot(value, spatialLag,I_global, axes_handle)
% This function draws the Moran Scatterplot in the axes handle specified.
% It returns the handle of the points in the scatterplot and that of the
% line of I_global.

%reshape input value and spatial lag in such a way that they are vectors
value = reshape(value, [],1);
spatialLag = reshape(spatialLag, [], 1);

maxValue = max(max(value));
minValue = min(min(value));
maxAbsValue = abs( max([maxValue minValue], [], 'ComparisonMethod','abs'));
maxSpatialLag = max(max(spatialLag));
minSpatialLag = min(min(spatialLag));
maxAbsSpatialLag = abs( max([maxSpatialLag minSpatialLag], [], 'ComparisonMethod','abs'));
maxAll = max([maxAbsValue maxAbsSpatialLag]);

moranColors = [1 0.41 0.16; 1 1 0; 0.07 0.62 1; 0 1 0 ];

verticesQuadrant1 = 1.1* [0 0; maxAll 0; maxAll maxAll; 0 maxAll];

%maxThreshold = max(max(abs(I_local)));

if isempty(axes_handle)
    figure
    axes_handle = axes;
end

title(axes_handle, 'Moran''s Scatterplot')
xlabel(axes_handle, "Pixel Value")
ylabel(axes_handle, "Spatial lag")
grid(axes_handle, 'on')
hold(axes_handle, 'on')

% xlim(1.1*[minA maxA]);
% ylim(1.1*[minB maxB]);
xlim(axes_handle, 1.05*[-maxAll maxAll]);
ylim(axes_handle, 1.05*[-maxAll maxAll]);
xticks(axes_handle, -ceil(maxAll):ceil(maxAll));
yticks(axes_handle, -ceil(maxAll):ceil(maxAll));
% pbaspect([1 1 1])

% draw the colored backgrounds
patch(axes_handle, verticesQuadrant1(:,1), verticesQuadrant1(:,2), moranColors(1,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.2);
patch(axes_handle, -verticesQuadrant1(:,1), verticesQuadrant1(:,2), moranColors(2,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.2);
patch(axes_handle, -verticesQuadrant1(:,1), -verticesQuadrant1(:,2), moranColors(3,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.2);
patch(axes_handle, verticesQuadrant1(:,1), -verticesQuadrant1(:,2), moranColors(4,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.2);

% draw the points 

%scatterplotHandle = plot(axes_handle, value, spatialLag, '.', 'color', 'k', 'MarkerSize', 0.5);
scatterplotHandle = scatter(axes_handle, value, spatialLag, 3, 'k', 'filled');
% draw the linear regression line given by the global Moran's I   
MoranLineHandle = plot(axes_handle, axes_handle.XLim, I_global*axes_handle.XLim, '-r');


if nargout > 0
    varargout{1} = scatterplotHandle;
    if nargout == 2
        varargout{2} = MoranLineHandle;
    end
end

end