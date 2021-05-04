function z = normz (z, dim)

if dim ~= 0
    if dim == 2
        zb = sum(z,1);
        for i = 1:size(z,2)
            z(:,i) = z(:,i)/zb(i);
        end
    elseif dim == 1
        zb = sum(z,2);
        for i = 1:size(z,1)
            z(i,:) = z(i,:)/zb(i);
        end
    end
    
    for i = 1:size(z,1)
        for j = 1:size(z,2)
            if isnan(z(i,j))
                z(i,j) = 0;
            end
        end
    end
end