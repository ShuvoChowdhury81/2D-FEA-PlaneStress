function PlotMesh(Nodes, Elements, U_global, scale_factor)
% =========================================================================
%  Function: PlotMesh
%  Purpose:  Plots the finite element mesh. Can show deformed shape.
%
%  Inputs:
%    Nodes        - N x 3 matrix [ID, x, y]
%    Elements     - M x 5 matrix [ID, n1, n2, n3, n4]
%    U_global     - (Optional) Displacement vector. Pass [] if not needed.
%    scale_factor - (Optional) Scale factor for deformation. Default = 0.
% =========================================================================

    % Handle optional arguments
    if nargin < 3
        U_global = []; 
        scale_factor = 0;
    elseif nargin < 4
        scale_factor = 1.0;
    end

    % 1. Prepare Vertices (Coordinates)
    % Extract X and Y columns (Cols 2 and 3)
    coords = Nodes(:, 2:3);
    
    % If displacements are provided, add them to coordinates
    if ~isempty(U_global) && scale_factor ~= 0
        % Reshape U into (N_nodes x 2) to match coords
        % U_global is [u1; v1; u2; v2...]
        
        disp_x = U_global(1:2:end); % Odd indices
        disp_y = U_global(2:2:end); % Even indices
        
        coords(:, 1) = coords(:, 1) + disp_x * scale_factor;
        coords(:, 2) = coords(:, 2) + disp_y * scale_factor;
        
        title_str = sprintf('Deformed Mesh (Scale: %.1f)', scale_factor);
        color_val = 'b'; % Blue for deformed
    else
        title_str = 'Undeformed Mesh';
        color_val = 'w'; % White for undeformed
    end

    % 2. Prepare Faces (Connectivity)
    % Extract Node IDs (Cols 2 to 5)
    faces = Elements(:, 2:5);

    % 3. Create Plot
    %figure;
    patch('Faces', faces, 'Vertices', coords, ...
          'FaceColor', color_val,'FaceAlpha', 0.7, ...
          'EdgeColor', 'k', ...     % Black edges
          'LineWidth', 1.0);
      
    axis equal;  % Crucial: Keeps x and y scales the same (circles look like circles)
    grid on;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    %title(title_str);
    
    % Optional: Number the nodes (good for debugging small meshes)
    if size(Nodes, 1) < 50
        for i = 1:size(Nodes, 1)
            text(coords(i,1), coords(i,2), num2str(Nodes(i,1)), ...
                 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end

end