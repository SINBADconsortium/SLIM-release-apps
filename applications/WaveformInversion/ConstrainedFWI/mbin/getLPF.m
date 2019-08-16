function LPF = getLPF(model, params, k)
    M = model.n(1);
    N = model.n(2);
    scale   = k * model.d(1);

    lbound = scale * params.smoothpars(1);
    ubound = scale * params.smoothpars(2);
    aspect = params.smoothpars(3)*(model.d(1)/model.d(2));
    LPF = lpFilter(M, N, lbound, ubound, aspect);

end
