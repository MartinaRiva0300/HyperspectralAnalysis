function blackNan(M, varargin)
    % This function creates a figure or plots in target axes a matrix in
    % which NaN values are coloured according to the appropriate argument.
    % The colormap is set by default to parula, but can also be changed.

    % arg1: the matrix to be drawn
    % arg2: colormap (matrix n*3 with RGB coordinates)
    % arg3: nan color (vector 1*3 with RGB coordinates)
    % arg4: target axes
    
    % Default values
    defaultCmap = parula;
    defaultColorNan = [0, 0, 0];
    
    % Initialize parameters with default values
    cmap = defaultCmap;
    colorNan = defaultColorNan;

    % Check the number of input arguments
    if nargin >= 2
        if ~isempty(varargin{1})
            cmap = varargin{1};
        end
    end
    
    if nargin >= 3
        if ~isempty(varargin{2})
            colorNan = varargin{2};
        end
    end
    
    % Normalize and convert matrix to indexed image
    indM = uint8(256 * mat2gray(M));
    rgbM = ind2rgb(indM, cmap);
    
    % Separate RGB channels
    rgbM1 = rgbM(:,:,1);
    rgbM2 = rgbM(:,:,2);
    rgbM3 = rgbM(:,:,3);
    
    % Replace NaN values with specified color
    nanMask = isnan(M);
    rgbM1(nanMask) = colorNan(1);
    rgbM2(nanMask) = colorNan(2);
    rgbM3(nanMask) = colorNan(3);
    
    % Combine channels back into an RGB image
    M_out = cat(3, rgbM1, rgbM2, rgbM3);
    
    % Plot the image
    if nargin >= 4
        if ~isempty(varargin{3})
            targetAxes = varargin{3};
            imagesc(targetAxes, M_out)
        end
    else % create a new figure and plot it there
        figure, imagesc(M_out);
    end
end

%% OLD VERSION
% function blackNan(M)
%     indM = uint8(256*mat2gray(M));
%     rgbM = ind2rgb(indM, parula);
%     rgbM1 = rgbM(:,:,1);
%     rgbM2 = rgbM(:,:,2);
%     rgbM3 = rgbM(:,:,3);
% 
%     rgbM1(isnan(M)) = 0;
%     rgbM2(isnan(M)) = 0;
%     rgbM3(isnan(M)) = 0;
% 
%     M_out = cat(3, rgbM1, rgbM2, rgbM3);
%     figure, imagesc(M_out)
% end
