%--------------------------------------------------------------------------
% This script implements the example in Appendix A, showing that the minimizer in the 
% hypothesis testing characterization is not taken at any positive operator.
%--------------------------------------------------------------------------
% author: Kun Fang (fangfred11@gmail.com)
% based on the paper (https://arxiv.org/abs/1706.06221)
%--------------------------------------------------------------------------

% input state initialization
da = 2;
db = 2;
theta = 0.5;
phi1 = [cos(theta) 0 0 sin(theta)]';
phi1 = phi1*phi1';
phi2 = [0 0 1 0]';
phi2 = phi2*phi2';
rho = 3/4*phi1 + 1/4*phi2;
% error tolerance
e = 1-sqrt(3)/2;

% optimization one
cvx_begin sdp quiet
    variable X(da*db, da*db) hermitian
    variable G(da*db, da*db) hermitian
    variable t
    r1 = trace(X) + t*(1-e);
    maximize r1
    subject to
        G - X - t*rho >= 0;
        X <= 0;
        t >= 0;
        SchattenNorm(PartialTranspose(G,2,[da db]),1) <= 1;
        G == G';
cvx_end
opt1 = -log2(r1);

% optimization two
cvx_begin sdp quiet
    variable X(da*db, da*db) hermitian
    variable G(da*db, da*db) hermitian
    variable t
    r2 = trace(X) + t*(1-e);
    maximize r2
    subject to
        G - X - t*rho >= 0;
        X <= 0;
        t >= 0;
        SchattenNorm(PartialTranspose(G,2,[da db]),1) <= 1;
        G >= 0;
cvx_end
opt2 = -log2(r2);