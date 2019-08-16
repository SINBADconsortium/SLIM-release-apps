function y = ifelse_func(true_func,false_func,condition,varargin)
    if condition
        y = true_func(varargin{:});
    else
        y = false_func(varargin{:});
    end
end