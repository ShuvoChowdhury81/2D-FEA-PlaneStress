% ===============================================
% ME 5310 -  Finite Element Project
% 2D Plane Stress Elasticity Solver
% File: Main_Project.m (Driver Script)
% ==========================================
clc; clear; close all;
%% 1. Input Data (Set up for test case A: Single Element Validation)
%-----------------------------------------------------------------
% Material Properties
E  = 30e6;
nu = 0.3;
t = 1; %Thickness

% Mesh Definition
% Nodes List: [Node ID, x-coord, y-coord]
Nodes = [
1,	0,	0	;
2,	1.0,	0.0	;
3,	1.0,		1.0;
4,	0.,	1.0	

    ];

% Element Connectivity: [Elem ID, Node 1, Node 2, Node 3, Node 4]
% Note: Nodes must be listed in counter-clockwise order.
Elements = [
    1,	1,	2,	3,	4
];

% Problem Size Parameters
num_nodes = size(Nodes, 1);
num_elems = size(Elements, 1);
num_dofs = num_nodes * 2; % 2 DOFs per node (u, v)

%% 2. INITIALIZATION
% -------------------------------------------------------------------------
% Initialize global stiffness matrix (Sparse is better for large problems)
K_global = sparse(num_dofs, num_dofs);
F_global = zeros(num_dofs, 1);
U_global = zeros(num_dofs, 1);

%% 3. ASSEMBLY PROCESS
% -------------------------------------------------------------------------
disp('Assembling Global Stiffness Matrix...');

for e = 1: num_elems
    % A. extract node IDs for this element 
    node_ids = Elements(e, 2:5);

    % B. Extract coordinates for this node
    ele_coords = Nodes(node_ids, 2:3);

    % C. Compute element stiffness matrix
    k_e = Quad4Stiffness (ele_coords, E, nu, t);

    % D. Mapping to the global matrix by a scatter vector

    scatter = zeros (1,8);

    for n = 1:4
        global_node= node_ids(n);
        scatter(2*n-1) = 2*global_node - 1; % x-dof (Odd index)
        scatter(2*n)   = 2*global_node;     % y-dof (Even index)
    end

    % E. Add to global matrix
    K_global(scatter, scatter) = K_global(scatter, scatter) + k_e;

end


%% 4. APPLY BOUNDARY CONDITIONS (BCs) & LOADS
% -------------------------------------------------------------------------

fixed_nodes_u = [1,4]; % Nodes fixed in X
fixed_nodes_v = [1,4]; % Nodes fixed in Y

% Convert Node IDs to DOF Indices
fixed_dofs = [ (2*fixed_nodes_u - 1), (2*fixed_nodes_v) ]; % [1, 2, 4, 7]
fixed_dofs = unique(fixed_dofs); % Remove duplicates % [1, 2, 4, 7]

% Free DOFs are all other DOFs
all_dofs = 1:num_dofs;  % [1, 2, 3, 4, 5, 6,7,8]
free_dofs = setdiff(all_dofs, fixed_dofs);   % [3, 5, 6,8]

% B. Apply External Forces (Point Loads)
% Test Case Setup: Apply Tension in Y-direction at top nodes (3 & 4)
% Load P = 1000 lbs split between top nodes.
P_load = 1971;  %Arbitrary
load_nodes = [2,3];

for i = 1:length(load_nodes)
    node = load_nodes(i);
    dof_x = 2*node-1; % y-dof
    F_global(dof_x) = F_global(dof_x) + P_load;
end

%% 5. SOLVER
% -------------------------------------------------------------------------
disp('Solving system equations...');

% Partition the system (Extract only free DOFs)
K_ff = K_global(free_dofs, free_dofs);
F_f  = F_global(free_dofs);

% Solve for unknown displacements [cite: 65]
U_f = K_ff \ F_f;

% Reconstruct the full displacement vector
U_global(free_dofs) = U_f;
U_global(fixed_dofs) = 0; % Enforce zero at fixed supports

disp('Nodal Displacements:');
disp(U_global);

% %% 6. POST-PROCESSING
% % -------------------------------------------------------------------------
% disp('Calculating Element Stresses...');
% 
% for e = 1:num_elems
%     node_ids = Elements(e, 2:5);
%     ele_coords = Nodes(node_ids, 2:3);
% 
%     % Extract displacements for this element
%     u_ele = zeros(8,1);
%     for n = 1:4
%         global_node = node_ids(n);
%         u_ele(2*n-1) = U_global(2*global_node - 1);
%         u_ele(2*n)   = U_global(2*global_node);
%     end
% 
%     % Calculate Stress
%     % calls the function CalcStress.m (To be written)
%     [sigma_vec, vm_stress] = CalcStress(ele_coords, u_ele, E, nu);
% 
%     % Print Results
%     if e==150
%         fprintf('Element %d Centroid Stress:\n', e);
%         fprintf('  Sig_xx: %.2f  Sig_yy: %.2f  Tau_xy: %.2f  VonMises: %.2f\n', ...
%                 sigma_vec(1), sigma_vec(2), sigma_vec(3), vm_stress);
%     end
% end


%% 6. POST-PROCESSING (Updated for Gauss Points)
% ---------------------------------------------
disp('Calculating Element Stresses at Integration Points...');

for e = 1:num_elems
    node_ids = Elements(e, 2:5);
    ele_coords = Nodes(node_ids, 2:3);

    % Extract displacements
    u_ele = zeros(8,1);
    for n = 1:4
        global_node = node_ids(n);
        u_ele(2*n-1) = U_global(2*global_node - 1);
        u_ele(2*n)   = U_global(2*global_node);
    end

    % Calculate Stress at 4 Gauss Points
    [stress_mat, vm_vec] = CalcStress_Gauss(ele_coords, u_ele, E, nu);

    if e==150
        %Print Results
        fprintf('\nElement %d Results:\n', e);
        fprintf('  GP | Sig_xx   | Sig_yy   | Tau_xy   | VonMises\n');
        fprintf('  ---|----------|----------|----------|---------\n');
        for gp = 1:4
            fprintf('  %d  | %8.2f | %8.2f | %8.2f | %8.2f\n', ...
                gp, stress_mat(gp,1), stress_mat(gp,2), stress_mat(gp,3), vm_vec(gp));
        end
    end
end


disp('Analysis Complete.');
disp('Plotting input mesh...');
PlotMesh(Nodes, Elements);

disp('Plotting deformed shape...');
% scale = 1000; % Exaggerate deformation to make it visible
% PlotMesh(Nodes, Elements, U_global, 1);
disp('Plotting Superimposed shape')

figure;
hold on;
% 1. Plot undeformed mesh
PlotMesh(Nodes, Elements);
scale = 10000; % Exaggerate deformation to make it visible
PlotMesh(Nodes, Elements, U_global, scale);

%title(['Superimposed Deformation (Scale: ', num2str(scale), ')']);
axis equal;
grid on;
hold off;    % Release the hold