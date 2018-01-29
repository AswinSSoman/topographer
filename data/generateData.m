
clear all
close all
clc

decay          = 0.25;
synth          = 250;
K_d_inhibition = 200;
K_d_activation = 400;
n_inhibition   = 2;
n_activation   = 20;

init = [300 300 0 0 0 0]';

sigma   = 100e-0;
dt      = 1;
T       = 240; 
a(:,1)  = init;
nRealiz = 400;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

MM = 100000;
for i = 1:MM
    rp = randperm(2); rp = rp(1);
    if rp==1
        np = randn/8 + 0.22;
    else
        np = randn/8 + 0.77;
    end
    Np(i) = np;
end
k = 1;
for i = 1:length(Np)
    if (Np(i)<1)&&(Np(i)>0)
        SampleDistribution(k) = Np(i);
        k = k+1;
    end
end
figure, hist(SampleDistribution,100), title('Sampling strategy')
MM = length(SampleDistribution);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

for j = 1:nRealiz
    
    clear a
    a(:,1) = init;
    
    for i = 1:T-1

        a(1,i+1) = a(1,i) + dt*(synth*(1/(1+ (a(2,i)/K_d_inhibition)^n_inhibition)) - decay*a(1,i)) + sqrt(sigma)*randn;
        a(2,i+1) = a(2,i) + dt*(synth*(1/(1+ (a(1,i)/K_d_inhibition)^n_inhibition)) - decay*a(2,i)) + sqrt(sigma)*randn;
        
        a(3,i+1) = a(3,i) + dt*(synth*(1/(1+(K_d_activation/a(1,i))^n_activation))*(1/(1+ (a(4,i)/K_d_inhibition)^n_inhibition)) - decay*a(3,i)) + sqrt(sigma)*randn;
        a(4,i+1) = a(4,i) + dt*(synth*(1/(1+(K_d_activation/a(1,i))^n_activation))*(1/(1+ (a(3,i)/K_d_inhibition)^n_inhibition)) - decay*a(4,i)) + sqrt(sigma)*randn;
    
        a(5,i+1) = a(5,i) + dt*(synth*(1/(1+(K_d_activation/a(2,i))^n_activation))*(1/(1+ (a(6,i)/K_d_inhibition)^n_inhibition)) - decay*a(5,i)) + sqrt(sigma)*randn;
        a(6,i+1) = a(6,i) + dt*(synth*(1/(1+(K_d_activation/a(2,i))^n_activation))*(1/(1+ (a(5,i)/K_d_inhibition)^n_inhibition)) - decay*a(6,i)) + sqrt(sigma)*randn;
        
        if sum(a([3 4 5 6],i+1)>950)
            ind = ceil(i*SampleDistribution(round(rand*MM)));
            DataDiff(1:6,j) = a(1:6,ind);
            break
        end
    
    end

end

figure, plot(a'), title('Single realisation from toggle switch network')

% save diffusionInput DataDiff
