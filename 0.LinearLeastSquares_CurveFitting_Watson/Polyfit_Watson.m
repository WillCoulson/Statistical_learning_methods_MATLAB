function [p] = Polyfit_Watson(fit_x,fit_y,n)

    x = fit_x';
    y = fit_y';

    V = ones(length(x),n+1);
    
    % Construct the Vandermonde matrix V = [x.^n ... x.^2 x ones(size(x))]
    V(:,n+1) = ones(length(x),1,class(x));
    for j = n:-1:1
        V(:,j) = x.*V(:,j+1);
    end

    % 线性最小二乘法拟合实现的核心 
    [Q,R] = qr(V,0);
    % 线性最小二乘法拟合实现的核心 
 
    p = R\(Q'*y);               % Same as p = V\y

    p = p.'; % Polynomial coefficients are row vectors by convention.

end
