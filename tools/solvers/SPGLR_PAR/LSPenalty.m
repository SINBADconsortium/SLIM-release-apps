function [f,g] = LSPenalty(X)
    f = 0.5 * norm(X,'fro')^2;
    g = X;
end    