% =========================================================================
%  Function: Quad4Stiffness
%  Purpose:  Calculates the 8x8 Element Stiffness Matrix for a 
%            4-Node Isoparametric Quadrilateral in Plane Stress.
%
%  Inputs:
%    coords - 4x2 matrix of node coordinates [x1 y1; x2 y2; x3 y3; x4 y4]
%    E      - Young's Modulus
%    nu     - Poisson's Ratio
%    t      - Thickness
%
%  Outputs:
%    ke     - 8x8 Element Stiffness Matrix
% =========================================================================

function [ke] = Quad4Stiffness(coords, E, nu, t)

    % 1. CONSTITUTIVE MATRIX (D) - Plane Stress Setup
    C = E/(1-nu^2);
    D = C* [1 nu 0;
            nu 1 0;
            0 0 (1-nu)/2];
    pt = 1 / sqrt(3);           % Location of points +/- 0.57735
    w  = 1.0;                   % Weight is 1.0 for all points in 2x2 rule

    % Gauss points (xi, eta)
    GP_xi  = [-pt,  pt,  pt, -pt];
    GP_eta = [-pt, -pt,  pt,  pt];

    % Initialize Element Stiffness Matrix
    ke = zeros(8, 8);
    
    % 3. Integration points

    for i = 1:4
        xi = GP_xi(i);
        eta = GP_eta(i);

        % N1 = 0.25*(1-xi)*(1-eta)
        % N2 = 0.25*(1+xi)*(1-eta)
        % N3 = 0.25*(1+xi)*(1+eta)
        % N4 = 0.25*(1-xi)*(1+eta)

        % dN_dxi (1x4 vector)
        dN_dxi  = 0.25 * [ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ];
        
        % dN_deta (1x4 vector)
        dN_deta = 0.25 * [ -(1-xi),  -(1+xi),   (1+xi),   (1-xi)  ];

        % B. JACOBIAN MATRIX (J)
        % J = [ dx/dxi   dy/dxi ]
        %     [ dx/deta  dy/deta]
        % computed as (dN_dnatural * coords)
        
        J = [dN_dxi; dN_deta] * coords;

        detJ = det(J); % Determinant of Jacobian (Area scaling factor)
        
        % Check for distorted elements (Area must be positive)
        if detJ <= 0
            error('Element Jacobian is non-positive. Check node numbering (Counter-Clockwise).');
        end

        % C. DERIVATIVES wrt PHYSICAL COORDS (x, y)
        % [dN_dx] = inv(J) * [dN_dxi ]
        % [dN_dy]            [dN_deta]
        
        %invJ = inv(J);
        dN_dxy = J\ [dN_dxi; dN_deta]; % invJ * [dN_dxi; dN_deta]; % 2x4 Matrix
        
        dN_dx = dN_dxy(1, :); % Top row is dN/dx
        dN_dy = dN_dxy(2, :); % Bottom row is dN/dy

        % D. ASSEMBLE B-MATRIX (Strain-Displacement Matrix) [cite: 62]
        % B is 3x8. Each node "j" fills 2 columns.
        % Format:
        % [ N1,x   0     N2,x   0     ... ]
        % [ 0      N1,y  0      N2,y  ... ]
        % [ N1,y   N1,x  N2,y   N2,x  ..

        B = zeros(3, 8);
        for n = 1:4
            col_u = 2*n - 1; % Column for u-displacement
            col_v = 2*n;     % Column for v-displacement
            
            B(1, col_u) = dN_dx(n);
            B(1, col_v) = 0;
            
            B(2, col_u) = 0;
            B(2, col_v) = dN_dy(n);
            
            B(3, col_u) = dN_dy(n);
            B(3, col_v) = dN_dx(n);
        end

        % E. STIFFNESS ACCUMULATION
        % ke = Sum ( B' * D * B * detJ * weight * thickness )
        ke = ke + (B' * D * B) * detJ * w * w * t;
    end
end

    
