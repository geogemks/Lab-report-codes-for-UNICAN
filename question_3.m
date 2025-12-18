%Re-writing the first code, as in the file question_2.m
mat1_uncleaned = readtable('SCPUnion2.1_AllSNe.txt', 'ReadVariableNames', false);

for col = 2:width(mat1_uncleaned)
    mat1_uncleaned{:, col} = str2double(string(mat1_uncleaned{:, col}));
end
if width(mat1_uncleaned) > 4
    mat1_uncleaned(:, end) = []; 
end

rowsToRemove = any(mat1_uncleaned{:, 2:end} == -99.9, 2);
data = mat1_uncleaned(~rowsToRemove, :);

z_obs = data{:, 2};       % Redshift
mag_obs = data{:, 3};     % Observed Magnitude
mag_err = data{:, 4};     % Magnitude Error


mu_obs = mag_obs + 19.3; 

c = 299792.458;   % km/s
H0 = 67.4;        % km/s/Mpc

Om_range = linspace(0, 1, 100); 
chi2_values = zeros(size(Om_range));

for i = 1:length(Om_range)
    Om_curr = Om_range(i);
    OL_curr = 1 - Om_curr; 
    
    mu_theory = zeros(size(z_obs));
    
    for j = 1:length(z_obs)
        integrand = @(z) 1./sqrt(Om_curr*(1+z).^3 + OL_curr);
        dL_Mpc = (1+z_obs(j)) * (c/H0) * integral(integrand, 0, z_obs(j));
        mu_theory(j) = 5 * log10(dL_Mpc) + 25;
    end
    
    chi2_values(i) = sum( ((mu_obs - mu_theory) ./ mag_err).^2 );
end

[min_chi2, idx_best] = min(chi2_values);
best_Om = Om_range(idx_best);
best_OL = 1 - best_Om;


fprintf('------------------------------------------------\n');
fprintf('BEST FIT RESULTS (Flat Universe, H0=67.4)\n');
fprintf('------------------------------------------------\n');
fprintf('Best Matter Density (Omega_m): %.4f\n', best_Om);
fprintf('Best Dark Energy (Omega_Lambda): %.4f\n', best_OL);
fprintf('Minimum Chi-Squared: %.2f\n', min_chi2);
fprintf('------------------------------------------------\n');


figure('Color','w');
plot(Om_range, chi2_values, 'k-', 'LineWidth', 2);
hold on;
plot(best_Om, min_chi2, 'ro', 'MarkerFaceColor', 'r');
xlabel('\Omega_m (Matter Density)');
ylabel('\chi^2 (Goodness of Fit)');
title('Chi-Squared Minimization for Flat Universe');
grid on;
