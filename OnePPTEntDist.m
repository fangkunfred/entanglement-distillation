%--------------------------------------------------------------------------
%   Calculate the one-shot distillable entanglement under PPT operations
%--------------------------------------------------------------------------
% This function has four required arguments:
%   RHO: a bipartite quantum state
%   DA: a real integer number
%   DB: a real integer number
%   E: a real number between zero and one
%--------------------------------------------------------------------------
% R = OnePPTEntDist(RHO,DA,DB,E) produces the one-shot distillable entanglement under PPT operations
% for the quantum state RHO with local dimension DA, DB, input arguement E
% indicates the error tolerance
%--------------------------------------------------------------------------
% requires: CVX (http://cvxr.com/cvx/) and QETLAB package (http://www.qetlab.com/Main_Page)
% author: Kun Fang (fangfred11@gmail.com)
% based on the paper (https://arxiv.org/abs/1706.06221)

function r = OnePPTEntDist(rho,da,db,e)

    cvx_begin sdp quiet
    cvx_precision best 
    variable M(da*db, da*db) hermitian
    variable eta    
    minimize eta
    subject to
        0 <= M <= eye(da*db);
        trace(rho*M) >= 1-e;
        -eta*eye(da*db) <= PartialTranspose(M,2,[da db]) <= eta*eye(da*db);    
    cvx_end
    r = log2(floor(1/eta));
    
end
