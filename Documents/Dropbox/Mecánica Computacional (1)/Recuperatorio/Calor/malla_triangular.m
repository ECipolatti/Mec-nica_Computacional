function [icone,xnod,B] = malla_triangular(Nx,XA,XB,Ny,YA,YB,bc,esq)
	coord_x = linspace(XA,XB,Nx);
    coord_y = linspace(YA,YB,Ny);
    
    [X,Y] = meshgrid(coord_x,coord_y);
    xnod = [reshape(X,Nx*Ny,1),reshape(Y,Nx*Ny,1)];

    icone = delaunay(X,Y);
    
    B = []; %#ok<*AGROW>
    
    % Ajuste de los bordes
    for i = 1 : Nx*Ny
        if xnod(i,1) == XA
            if xnod(i,2) ~= YA && xnod(i,2) ~= YB
                B = [B; i,bc(1,:)];
            elseif xnod(i,2) == YA
                B = [B; i,esq(1,:)];
            else
                B = [B; i,esq(2,:)];
            end
        end
        
        if xnod(i,1) == XB
            if xnod(i,2) ~= YA && xnod(i,2) ~= YB
                B = [B; i,bc(2,:)];
            elseif xnod(i,2) == YA
                B = [B; i,esq(3,:)];
            else
                B = [B; i,esq(4,:)];
            end
        end
        
        if xnod(i,2) == YA
            if xnod(i,1) ~= XA && xnod(i,1) ~= XB
                B = [B; i,bc(3,:)];
            end
        end
        
        if xnod(i,2) == YB
            if xnod(i,1) ~= XA && xnod(i,1) ~= XB
                B = [B; i,bc(4,:)];
            end
        end
    end 
end
