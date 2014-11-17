function [y, performance] = test(model, Xtest, varargin)
    if nargout > 1
        [y, performance] = gurls_test(model, Xtest, varargin{:});
    else
       y = gurls_test(model, Xtest, varargin{:});
    end
end
