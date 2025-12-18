% Best Fit Cosmological Parameters (Flat Universe)
% Fits Om (Matter Density) assuming Flatness (Om + OL = 1)
% Uses H0 = 67.4 (Planck 2020)

clear; clc; close all;

% 1. Load Data
% Ensure 'SCPUnion2.1_AllSNe.txt' is in your folder
mat1_uncleaned = readtable('SCPUnion2.1_AllSNe.txt', 'ReadVariableNames', false);

% Convert columns to double
for col = 2:width(mat1_uncleaned)
    mat1_uncleaned{:, col} = str2double(string(mat1_uncleaned{:, col}));
end
if width(mat1_uncleaned) > 4
    mat1_uncleaned(:, end) = []; 
end

% Remove invalid rows
rowsToRemove = any(mat1_uncleaned{:, 2:end} == -99.9, 2);
data = mat1_uncleaned(~rowsToRemove, :);

% Extract variables
z_obs = data{:, 2};       % Redshift
mag_obs = data{:, 3};     % Observed Magnitude
mag_err = data{:, 4};     % Magnitude Error

% Absolute Magnitude Correction (Approximate)
% Distance Modulus mu = m - M. 
% We fit mu_theoretical to (mag_obs + 19.3).
mu_obs = mag_obs + 19.3; 


% 2. Constants
c = 299792.458;   % km/s
H0 = 67.4;        % km/s/Mpc

% 3. Fitting Loop
% We will test Omega_m values from 0 to 1
Om_range = linspace(0, 1, 100); 
chi2_values = zeros(size(Om_range));

for i = 1:length(Om_range)
    Om_curr = Om_range(i);
    OL_curr = 1 - Om_curr; % Enforce Flatness (K=0)
    
    % Calculate Theory mu for this Omega_m
    mu_theory = zeros(size(z_obs));
    
    for j = 1:length(z_obs)
        % Integral for luminosity distance in Flat Universe
        integrand = @(z) 1./sqrt(Om_curr*(1+z).^3 + OL_curr);
        dL_Mpc = (1+z_obs(j)) * (c/H0) * integral(integrand, 0, z_obs(j));
        mu_theory(j) = 5 * log10(dL_Mpc) + 25;
    end
    
    % Calculate Chi-Squared for this model
    % Chi2 = sum( (Observed - Theory)^2 / Error^2 )
    chi2_values(i) = sum( ((mu_obs - mu_theory) ./ mag_err).^2 );
end

% 4. Find Best Parameter
[min_chi2, idx_best] = min(chi2_values);
best_Om = Om_range(idx_best);
best_OL = 1 - best_Om;

% 5. Output Results
fprintf('------------------------------------------------\n');
fprintf('BEST FIT RESULTS (Flat Universe, H0=67.4)\n');
fprintf('------------------------------------------------\n');
fprintf('Best Matter Density (Omega_m): %.4f\n', best_Om);
fprintf('Best Dark Energy (Omega_Lambda): %.4f\n', best_OL);
fprintf('Minimum Chi-Squared: %.2f\n', min_chi2);
fprintf('------------------------------------------------\n');

% 6. Plot the Chi-Squared Curve
figure('Color','w');
plot(Om_range, chi2_values, 'k-', 'LineWidth', 2);
hold on;
plot(best_Om, min_chi2, 'ro', 'MarkerFaceColor', 'r'); % Mark the minimum
xlabel('\Omega_m (Matter Density)');
ylabel('\chi^2 (Goodness of Fit)');
title('Chi-Squared Minimization for Flat Universe');
grid on;
