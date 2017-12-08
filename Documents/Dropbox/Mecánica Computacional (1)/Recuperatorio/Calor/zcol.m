function [x] = zcol(r,N)
    % Devuelve un vector columna de N elementos con un '1' en la fila 'r' y '0' en el resto
    if r > N
        r = N;
    end
    
    if r > 1
        x = [zeros(r-1,1); 1; zeros(N-r,1)];
    else
        x = eye(N,1);
    end
end

