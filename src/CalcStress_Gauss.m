function [stress_matrix, vm_vec] = CalcStress_Gauss(coords, u_ele, E, nu)
% =========================================================================
%  Function: CalcStress_Gauss
%  Purpose:  Calculates stress at ALL 4 Integration Points (Gauss Points)
%
%  Inputs:
%    coords - 4x2 matrix of node coordinates
%    u_ele  - 8x1 vector of element displacements
%    E, nu  - Material properties
%
%  Outputs:
%    stress_matrix - 4x3 matrix containing stress at each Gauss point.
%                    Rows: Gauss Point 1 to 4
%                    Cols: [Sig_xx, Sig_yy, Tau_xy]
%    vm_vec        - 4x1 vector of Von Mises stress at each Gauss point
% =========================================================================

    % 1. CONSTITUTIVE MATRIX (D)
    C = E / (1 - nu^2);
    D = C * [ 1   nu  0 ;
              nu  1   0 ;
              0   0   (1-nu)/2 ];

    % 2. GAUSS POINTS (Same as Stiffness Matrix)
    pt = 1 / sqrt(3); 
    GP_xi  = [-pt,  pt,  pt, -pt];
    GP_eta = [-pt, -pt,  pt,  pt];
    
    % Initialize Outputs
    stress_matrix = zeros(4, 3);
    vm_vec = zeros(4, 1);

    % 3. LOOP OVER GAUSS POINTS
    for i = 1:4
        xi  = GP_xi(i);
        eta = GP_eta(i);
        
        % A. Derivatives at this Gauss Point
        dN_dxi  = 0.25 * [ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ];
        dN_deta = 0.25 * [ -(1-xi),  -(1+xi),   (1+xi),   (1-xi)  ];
        
        % B. Jacobian
        J = [dN_dxi; dN_deta] * coords;
        invJ = inv(J);
        
        % C. Global Derivatives
        dN_dxy = invJ * [dN_dxi; dN_deta];
        dN_dx = dN_dxy(1, :);
        dN_dy = dN_dxy(2, :);
        
        % D. B-Matrix Construction
        B = zeros(3, 8);
        for n = 1:4
            col_u = 2*n - 1;
            col_v = 2*n;
            B(1, col_u) = dN_dx(n);
            B(1, col_v) = 0;
            B(2, col_u) = 0;
            B(2, col_v) = dN_dy(n);
            B(3, col_u) = dN_dy(n);
            B(3, col_v) = dN_dx(n);
        end
        
        % E. Calculate Stress
        epsilon = B * u_ele;
        sigma = D * epsilon;
        
        % Store Sig_xx, Sig_yy, Tau_xy
        stress_matrix(i, :) = sigma'; 
        
        % Calculate Von Mises for this point
        sig_x = sigma(1);
        sig_y = sigma(2);
        tau_xy = sigma(3);
        vm_vec(i) = sqrt(sig_x^2 + sig_y^2 - sig_x*sig_y + 3*tau_xy^2);
    end
end