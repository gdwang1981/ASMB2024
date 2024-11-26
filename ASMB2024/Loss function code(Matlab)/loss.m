clc
clear
t0=50;
t1=5000;
data=xlsread('beta and eta.xlsx');
beta_matrix = data(:,2:13);
eta_matrix = data(:,14:25);

for group_index = 1:15
    eta = eta_matrix(group_index, :); 
    beta = beta_matrix(group_index, :); 
    
    Et = ones([length(eta), 1]); 
    for i = 1:length(eta)
        miu(i) = log(eta(i));
        sigma(i) = 1./beta(i);
        w_0(i) = (t0./exp(miu(i))).^(1./sigma(i));
        w_1(i) = (t1./exp(miu(i))).^(1./sigma(i));
        Em(i) = 1 - exp(-w_0(i));
        Ec(i) = t1./(t1 - t0).*(exp(-w_0(i)) - exp(-w_1(i))) - 1/(t1 - t0).*exp(miu(i)).*gamma(1 + sigma(i)).*(chi2cdf(2*w_1(i), 2*(1 + sigma(i))) - chi2cdf(2*w_0(i), 2*(1 + sigma(i))));
        Et(i) = Em(i) + Ec(i);
    end
    
    logloss(group_index, :, 1) = log(Et * 100); 
end