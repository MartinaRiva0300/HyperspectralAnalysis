function neighbors = getNeighbors(index, sz, range)
% index is the index of the element we want to find the neighbors of
% sz is the size of the original matrix as a vector
% range is the range within which we are looking for the neighbors


    % Convert the linear index to row and column subscripts
    rows = sz(1);
    columns = sz(2);
    [row, col] = ind2sub([rows, columns], index);
    
    % Initialize an empty array to store neighbor indices
    neighbors = [];
    
    % Loop through each possible neighbor within the specified range
    for dRow = -range:range
        for dCol = -range:range
            % Skip the center element itself
            if dRow == 0 && dCol == 0
                continue;
            end
            
            newRow = row + dRow;
            newCol = col + dCol;

            % Check if the new position is within matrix bounds
            if all([newRow >= 1, newRow <= rows, newCol >= 1, newCol <= columns])
                % Convert the subscript to a linear index and add to the list
                neighborIndex = sub2ind([rows, columns], newRow, newCol);
                neighbors = [neighbors; neighborIndex];
            end

        end
    end
end


%%
% function neighbors = getNeighbors(index, rows, columns)
%     % Convert the linear index to row and column subscripts
%     [row, col] = ind2sub([rows, columns], index);
% 
%     % Define the relative positions of the 8 neighbors
%     dRow = [-1, -1, -1, 0, 0, 1, 1, 1];
%     dCol = [-1, 0, 1, -1, 1, -1, 0, 1];
% 
%     % Initialize an empty array to store neighbor indices
%     neighbors = [];
% 
%     % Loop through each neighbor position
%     for k = 1:length(dRow)
%         newRow = row + dRow(k);
%         newCol = col + dCol(k);
% 
%         % Check if the new position is within matrix bounds
%         if newRow >= 1 && newRow <= rows && newCol >= 1 && newCol <= columns
%             % Convert the subscript to a linear index and add to the list
%             neighborIndex = sub2ind([rows, columns], newRow, newCol);
%             neighbors = [neighbors; neighborIndex];
%         end
%     end
% end
