function K = rbf2(X, X2, sigma)
    %%% function to generate Gaussian Kernel Matrix with 2 X's
    % https://code.google.com/p/altoolbox/source/browse/kernelmatrix.m
    % Gustavo Camps-Valls
    % 2006(c)
    % Jordi (jordi@uv.es), 2007
    % 2007-11: if/then -> switch, and fixed RBF kernel
    
    % return a n1 x n1 OR n2 x n1 matrix
    n1sq = sum(X.^2, 2);
    n1 = size(X, 1);

    if isempty(sigma) 
        D = ones(n1,1) * n1sq' + n1sq * ones(1, n1) -2 * X*X';
        sigma = ones(1, n1) * (D).^0.5 * ones(n1, 1)/n1^2;
    end
    
    if isempty(X2);
        D = ones(n1,1) * n1sq' + n1sq * ones(1, n1) -2 * X*X';
        %sigma = ones(1, n1) * (D).^0.5 * ones(n1, 1)/n1^2;
    else
        n2sq = sum(X2.^2, 2);
        n2 = size(X2, 1);
        D = ones(n2,1) * n1sq' + n2sq * ones(1, n1) -2 * X2*X';
        %sigma = ones(1, n2) * (D).^0.5 * ones(n1, 1)/(n1 * n2); 
    end
    
    K = exp(-D/(2*sigma^2));
end

