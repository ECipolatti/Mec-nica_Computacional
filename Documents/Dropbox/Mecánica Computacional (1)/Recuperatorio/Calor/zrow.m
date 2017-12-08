function [x] = zrow(c,N)
    % Devuelve un vector fila de N elementos con un '1' en la columna 'c' y '0' en el resto
    if c > N
        c = N;
    end
    
    if c > 1
        x = [zeros(1,c-1) 1 zeros(1,N-c)];
    else
        x = eye(1,N);
    end
end

