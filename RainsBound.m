%--------------------------------------------------------------------------
%              Calculate the Rains bound within given tolerence
%--------------------------------------------------------------------------
% This algorithm is completely based on the one by Zinchenko,Friedland,Gour [1]
% and Girard,Zinchenko,Friedland,Gour [2] with slight adjustments.
% In [1,2], they intend to calculate relative entropy of entanglment w.r.t. PPT states. 
% We adapt their script to calculating the Rains bound. 
%--------------------------------------------------------------------------
%  requires: CVX (http://cvxr.com/cvx/) and QETLAB package (http://www.qetlab.com/Main_Page)
%--------------------------------------------------------------------------
%   Standard usage: [Xopt,rains]=RainsBound(m,n,A)
%   Variables:
%        m,n       - dimensions of the subsystems
%        A         - density matrix whose Rains bound we are trying to compute
%        Xopt      - optimal PPT' matrix that minimizes the relative entropy
%        rains     - output upper bound of Rains bound
%
%   Optional inputs with defaults:
%   RainsBound(m,n,A,eps,maxIter,lineSearchEps)
%        eps           - precision such that |uBound-lBound|<eps 
%                        (default: eps = 1e-5)
%        maxIter       - max number of iterations 
%                        (default: maxIter = 200)
%        lineSearchEps - precision of intermediate line search 
%                        (default: lineSearchEps = 1e-10;)
%
% We define a function traceAlogmA to compute trace(A*logm(A))to accept
% rank-defficnent matrices.
%--------------------------------------------------------------------------
%   References:
%       [1] Zinchenko, Y., Friedland, S., & Gour, G. (2010).
%           Numerical estimation of the relative entropy of entanglement. 
%           Physical Review A, 82(5), 052336.
%       [2] Girard, M. W., Zinchenko, Y., Friedland, S., & Gour, G. (2015).
%           Erratum: Numerical estimation of the relative entropy of entanglement
%           [Phys. Rev. A 82, 052336 (2010)]. Physical Review A, 91(2), 029901. 
%--------------------------------------------------------------------------
% author: Kun Fang (fangfred11@gmail.com)

function [Xopt,uBound,lBound,outerCount] = RainsBound(m,n,A,eps,maxIter,lineSearchEps)

%--------------------------------------------------------------------------
%                       Check the input arguments
%--------------------------------------------------------------------------
if nargin<3
    error('Not enough arguments; input m,n,A')
end
% check if A is indeed positive semi-definite
if min(eig(A))<-10^(-10) || max(max(abs(A-A')))>1e-12 || max(size(A)~=[m*n,m*n])
    error('A must be positive semi-definite (mn x mn)-matrix');
end
% check if A is trace-1 (within some numerical tolerance level)
if abs(1-trace(A))>1e-12
    error('A must be trace-1; |1-trace(A)| exceeds the allowed 1e-12-tolerance');
end 

%--------------------------------------------------------------------------
%                 Set the optional parameters if not specified
%--------------------------------------------------------------------------
% If not specified, set the default precision.
if nargin<4
    eps=1e-5;
end
% if not specified, set the default maximum number of outer iterations.
if nargin<5
    maxIter=200;
end
% If not specified, set the default line-search precision.
if nargin<6
    lineSearchEps=1e-10;
end

%--------------------------------------------------------------------------
%                    Check if input state belongs to PPT' 
%--------------------------------------------------------------------------
if SchattenNorm(PartialTranspose(A,2,[m,n]),1) <= 1
    lBound=0;
    uBound=0;
    Xopt=A;
    disp('A belongs to Rains set, thus Xopt=A and RainsBound=0');
    return;
end

%--------------------------------------------------------------------------
%                        Initialize search
%--------------------------------------------------------------------------

% set blank output
lBound=-Inf;
uBound=Inf;
Xopt=[];
% set outer iteration counter
outerCount=0;
% set initial numerical status
status=0;
% initialze N
N=0;

% Initialize list of X. 
% This is the list of points in the interior of the PPT' states. 
% We construct the tangent planes at each of these points 
% to create a polytope approximation of the epigraph.
% X{N+1}=eye(m*n)/(m*n);
X{N+1}=Tensor(RandomDensityMatrix(m,1),RandomDensityMatrix(n,1));

% Generate list of E.
% Used for the approximate epigraph of the objective. These generate the
% 'Gateaux derivatives' of the tr(A*logm(X)) function at each X.
% (Here we provide a generic script that would work for any N, 
% although for now we always start with N=0)

for i=0:N
    % build E^{(i)} from X^{(i)}
    [U,L]=eig(X{i+1});
    for j=1:m*n
        for k=1:m*n
            if L(j,j)==L(k,k) D(j,k)=1/L(j,j);
            else D(j,k)=(log(L(j,j))-log(L(k,k)))/(L(j,j)-L(k,k));
            end
        end
    end
    E{i+1}=U*(D.*(U'*A*U))*U';
end

% symmetrize E
for i=1:N+1
    E{i}=(E{i}+E{i}')/2;
end

% Define list of b's. 
%   These are the values b{i}=-trace(A*logm(X))+trace(E*X) for each X and E.
for i=1:N+1
    b{i}=-trace(A*logm(X{i}))+trace(E{i}*X{i});
end
% Make into a vector so it can be used in cvx
bvect=zeros([N+1,1]);
for i=1:N+1
    bvect(i)=b{i};
end

% Start iterating until we reach the prescribed precision eps
% or unitl we exceed maximum number of iterations, i.e., outerCount>maxIter

% Set bestN index to point to the best upper bound out of X{i}, i=0,...,N
bestN=N;
% Re-initialize the bounds (we know something already)
lBound=-real(traceAlogmA(A));
uBound=-trace(A*logm(eye(m*n)/(m*n)));
% Set probSolved_flag to indicate wether the problem is solved yet
probSolved_flag=0;

%--------------------------------------------------------------------------
%                       Start optimization program
%--------------------------------------------------------------------------
while ~probSolved_flag && outerCount<maxIter %&& ~status 
    % Formulate and solve the approximation SDP problem. 
    % Variables are for SDP are:
    %     (Y,t)   - where Y is n*m by n*m PPT' matrix
    %               and (Y,t) is in approximation to epigraph of -trace(A*log(Y))
    %               so that t>=lBound.
    %     s       - s>=0
    % Main constraint is
    %     -trace(A*logm(X{i}))+trace(E{i}*X{i})-trace(E{i}*Y) <= t
    
    % Use CVX to solve SDP problem
    cvx_begin sdp quiet
%     cvx_precision high;
    variable t 
    variable s(N+1)
    variable Y(m*n,m*n) hermitian
    variable Yp(m*n,m*n) hermitian
    variable Ym(m*n,m*n) hermitian
    expression V(N+1)
    for i=1:N+1
        V(i)=trace(E{i}*Y);
    end
    minimize t
    subject to
        t>=lBound;
        s>=0;
        PartialTranspose(Y,1,[m,n])== Yp - Ym; % Rains set condition
        trace(Yp) + trace(Ym) <= 1; % Rains set condition
        Y>= 0; Yp>=0; Ym>= 0;
        s+bvect-V==t*ones([N+1,1]);        
    cvx_end
    
    % After SDP, use line search to find better optimal.
    Ystart=X{bestN+1};
    Yend=Y;
    % set the search ray direction
    dY=Yend-Ystart;
    % set the mid point and the objective derivative
    Ynext=(Ystart+Yend)/2;
    [U,L]=eig(Ynext);
    for j=1:m*n
        for k=1:m*n
          if L(j,j)==L(k,k) D(j,k)=1/L(j,j);
          else D(j,k)=(log(L(j,j))-log(L(k,k)))/(L(j,j)-L(k,k));
          end
        end
    end
    Enext=U*(D.*(U'*A*U))*U';
    df=-trace(Enext*dY);
    % iterate (with 'cheap' norm)
    while (norm(Yend(:)-Ystart(:))>lineSearchEps)
       if df<0 Ystart=Ynext;
       else Yend=Ynext;
       end
       % recompute the mid point and the objective derivative
       Ynext=(Ystart+Yend)/2;
       [U,L]=eig(Ynext);
       for j=1:m*n
            for k=1:m*n
                 if L(j,j)==L(k,k) D(j,k)=1/L(j,j);
                 else D(j,k)=(log(L(j,j))-log(L(k,k)))/(L(j,j)-L(k,k));
                 end
            end
       end
       Enext=U*(D.*(U'*A*U))*U';
       df=-trace(Enext*dY);
    end
    % Use Xnext for next point in list of P
    Y=Ynext;
    N=N+1;
    X{N+1}=Y;
    lBound=max(lBound,t);
    if uBound>-trace(A*logm(Y))
         uBound=-trace(A*logm(Y));
         bestN=N;
    end
    % output iteration result
    s=sprintf('[%d] lower bound: %e, upper bound: %e, gap: %e, relGap: %d%%',outerCount,lBound,uBound,uBound-lBound,round(100*(uBound-lBound)/uBound));
    disp(s);

    [U,L]=eig(X{N+1});
    for j=1:m*n
        for k=1:m*n
           if L(j,j)==L(k,k) D(j,k)=1/L(j,j);
           else D(j,k)=(log(L(j,j))-log(L(k,k)))/(L(j,j)-L(k,k));
           end
        end
    end
    E{N+1}=U*(D.*(U'*A*U))*U';

    % Symmetrize E
    E{N+1}=(E{N+1}+E{N+1}')/2;

    % Define vector of b's
    b{N+1}=real(-trace(A*logm(X{N+1}))+trace(E{N+1}*X{N+1}));
    bvect=zeros([N+1,1]);
    for i=1:N+1
       bvect(i)=b{i};
    end
    
    % Verify if we found a solution 
    if (uBound-lBound)<eps 
        probSolved_flag=1;
    end

    % Increment outer iteration counter
    outerCount=outerCount+1;
end

%--------------------------------------------------------------------------
%                        Output final result
%--------------------------------------------------------------------------
% Set the output
lBound=lBound+traceAlogmA(A);
lBound = lBound/log(2);
uBound=real(uBound+traceAlogmA(A));
uBound = uBound/log(2);
Xopt=X{bestN+1};

%--------------------------------------------------------------------------
% define a function traceAlogmA to compute trace(A*logm(A))to accept
% rank-defficnent matrices.
%--------------------------------------------------------------------------
function [val] = traceAlogmA(A)
% overloads MATLAB's trace(A*logm(A)) to accept
% rank-defficient positive semi-definite A
% by computing the limiting value
D=eig(A);
idx=find(D);
val=sum(D(idx).*log(D(idx)));

