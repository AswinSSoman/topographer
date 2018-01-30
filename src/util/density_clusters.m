function [ assign, peaks ] = density_clusters( dist, density, k, HARD )
%DENSITY_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT:
%   HARD: if true, the result cluster number would be precicely k. 

if ~exist('HARD', 'var'), HARD = false; end
BY_DISTANCE = false;

% get density peaks
[ peaks, ~, gamma, pointer ] = density_peaks( dist, density, k );

if HARD && (numel(peaks) > k)
    % restrict the cluster number to k
    sorted_peaks = sort_by(peaks, gamma(peaks), 'descend');
    peaks = sorted_peaks(1:k);
end

% fix the pointer for peaks; they should point to themself !
pointer(peaks) = peaks;

% assgien each point to a closest peak
if BY_DISTANCE
    
    [~, assign] = min( dist(:, peaks) , [] , 2 );
    
else
    
    [~, rho_order] = sort(density, 'descend');

    % update the link until all pointer points to a peak
    for c = rho_order(:)'      
        %{ 
        By doing it in density order, we ensure that this operation would be 
        efficient 
        %}
        while ~ismember(pointer(c), peaks)
            pointer(c) = pointer( pointer(c) );
        end 
    end
    
    % make up the result 
    assign = zeros(numel(density), 1);
    for i = 1:numel(peaks)
        p = peaks(i);
        assign(pointer == p) = i;
    end
    
end

end

