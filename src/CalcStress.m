function [sigma, vm_stress] = CalcStress(coords, u_ele, E, nu)
% =========================================================================
%  Function: CalcStress
%  Purpose:  Calculates stress at the element CENTROID (xi=0, eta=0)
%            for a Q4 element in Plane Stress.
%
%  Inputs:
%    coords - 4x2 matrix of node coordinates [x1 y1; ... ; x4 y4]
%    u_ele  - 8x1 vector of element displacements [u1; v1; ...; u4; v4]
%    E      - Young's Modulus
%    nu     - Poisson's Ratio
%
%  Outputs:
%    sigma     - 3x1 vector [Sig_xx; Sig_yy; Tau_xy]
%    vm_stress - Scalar Von Mises Equivalent Stress
% =========================================================================

    % 1. CONSTITUTIVE MATRIX (D) - Plane Stress
    % Same D matrix as in Stiffness calculation
    C = E / (1 - nu^2);
    D = C * [ 1   nu  0 ;
              nu  1   0 ;
              0   0   (1-nu)/2 ];

    % 2. EVALUATE AT CENTROID
    % We calculate stress at the center of the element.
    xi = 0;
    eta = 0;

    % 3. DERIVATIVES OF SHAPE FUNCTIONS (at Centroid)
    % For Q4, at (0,0), the derivatives are simple constants:
    % dN_dxi  = [-0.25,  0.25,  0.25, -0.25]
    % dN_deta = [-0.25, -0.25,  0.25,  0.25]
    
    dN_dxi  = 0.25 * [ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ];
    dN_deta = 0.25 * [ -(1-xi),  -(1+xi),   (1+xi),   (1-xi)  ];

    % 4. JACOBIAN AND MAPPING
    J = [dN_dxi; dN_deta] * coords;
    
    detJ = det(J);
    invJ = inv(J);
    
    % Transform derivatives to global x,y
    dN_dxy = invJ * [dN_dxi; dN_deta];
    dN_dx = dN_dxy(1, :);
    dN_dy = dN_dxy(2, :);

    % 5. ASSEMBLE B-MATRIX (Strain-Displacement)
    % Size: 3x8
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

    % 6. CALCULATE STRAIN & STRESS
    % Strain = B * u
    epsilon = B * u_ele;
    
    % Stress = D * Strain
    sigma = D * epsilon;
    
    % Extract components for clarity
    sig_x = sigma(1);
    sig_y = sigma(2);
    tau_xy = sigma(3);

    % 7. CALCULATE VON MISES STRESS
    % Formula for 2D Plane Stress
    vm_stress = sqrt(sig_x^2 + sig_y^2 - sig_x*sig_y + 3*tau_xy^2);

end