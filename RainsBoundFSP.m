%--------------------------------------------------------------------------
%   Calculate the Rains bound based on CVXQUAD package
%--------------------------------------------------------------------------
% This function has one required arguments:
%   rho: a bipartite quantum state
%--------------------------------------------------------------------------
% [r tau] = RainsBoundFazi(rho) produces the Rains bound of the quantum
% state rho. It outputs the Rains bound r and the minimizer tau.
%--------------------------------------------------------------------------
% requires: 
%    CVX (http://cvxr.com/cvx/) 
%    QETLAB (http://www.qetlab.com/Main_Page)
%    CVXQUAD (https://github.com/hfawzi/cvxquad)
%--------------------------------------------------------------------------
% References: 
% [1] Fawzi, H., Saunderson, J. and Parrilo, P.A., 2016. 
%     Semidefinite approximations of the matrix logarithm. 
%     Foundations of Computational Mathematics, pp.1-38.
% [2] Fawzi, H. and Fawzi, O., 2018. 
%     Efficient optimization of the quantum relative entropy. 
%     Journal of Physics A: Mathematical and Theoretical, 51(15), p.154003.
%--------------------------------------------------------------------------
% author: Kun Fang (fangfred11@gmail.com)

function [r tau] = RainsBoundFSP(rho)
    d = size(rho);
    d = d(1);
    cvx_begin sdp quiet
    variable tau(d,d) hermitian
    minimize (quantum_rel_entr(rho,tau)/log(2));
    subject to
        tau >= 0; 
        SchattenNorm(PartialTranspose(tau,2),1) <= 1;
    cvx_end
    r = cvx_optval;
end