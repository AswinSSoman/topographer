function rho = local_density(dist, R, STRICT_RHO, LOCAL_WIDHTH)

if ~exist('STRICT_RHO', 'var'), STRICT_RHO = false; end
if ~exist('LOCAL_WIDHTH', 'var'), LOCAL_WIDHTH = false; end

assert(issymmetric(dist));
N = size(dist, 1);
rho = zeros(1,N);

% set dist function
if STRICT_RHO, % use the strict definition of 'rho'
    dFunc = @(d, sigma) double( d <= sigma );
else % use Gaussian kernel. smoother.
    dFunc = @(d, sig) exp(-d^2 / (2 * sig^2) ) / sig;
end


if ~LOCAL_WIDHTH % regular choice !
    
    for i = 1:N-1
        for j = i+1:N
            dDen = dFunc(dist(i, j), R);
            rho(i) = rho(i) + dDen;
            rho(j) = rho(j) + dDen;
        end
    end
    
else % dynamic width, don't use R 
    
    for i = 1:N
        % find localized width
        if isinteger(LOCAL_WIDHTH)
            q = double(LOCAL_WIDHTH) / N;
            sigma = quantile(dist(i, :), q);
        else
            k = LOCAL_WIDHTH;
            ordereddist = sort(dist(i, :));
            sigma = ordereddist(k);
        end
        % get density with local sigma
        for j = 1:N
            dDen = dFunc(dist(i, j), sigma);
            rho(i) = rho(i) + dDen;
        end
    end
    
end

end