function [ D ] = mat_stress( E,nu )
    % Tensión Plana
    D = E/(1-nu^2) * [  1, nu,     0; 
                        nu, 1,     0; 
                        0,  0, (1-nu)/2  ];
    
    % Deformación Plana
%     D = E/((1+nu)*(1-2*nu))*[   1-nu, nu,   0;
%                                 nu,   1-nu, 0;
%                                 0,    0,  (1-2*nu)/2 ]; 

end

