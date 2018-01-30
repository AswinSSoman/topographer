function [ mrna, prob ] = fspstat( param, max_level )
% generate the expression distribution from given parameters. 
% 
% INPUT:
%   max_level: the maximum possible expression level of mRNA   
%   param.N: the number of states. For the simplest model, N = 2. 
%   param.mu: a vector of N; trancription rates for each state. 
%   param.A: a N-by-N transition matrix, corresponding to number fo states. 
%   param.delta: the degradation rate of mRNA transcripts. Set to 1 foe time
%       scale-invariant results. 
%

% originally, A(i, j) represents the rate from state-i to state-j for i!=j 
param.A = param.A'; % transpose it to be ready for right multiplication 
param.A = param.A - diag(sum(param.A)); % now represents drift in probability ?

D = param.delta * eye(param.N, param.N); 
T = diag(param.mu);

% set up matrix Q 
K = max_level + 1;
for i = 1:K
    for j = 1:K
        P{i, j} = zeros(param.N, param.N);
    end
end
P{1, 1} = param.A - T;
P{1, 2} = D;
for iter = 2:K-1
    P{iter, iter-1} = T;
    P{iter, iter  } = param.A - T - (iter-1) * D; 
    P{iter, iter+1} = iter * D;
end
P{K, K-1} = T;
P{K, K} = param.A - T - (K-1) * D;
Q = cell2mat(P);

% solve for steady state 
p0 = [Q; ones(1, K * param.N)] \ [zeros(K * param.N, 1); 1];
mrna = 0:K-1;
prob = sum( reshape(p0, param.N, K) );

end