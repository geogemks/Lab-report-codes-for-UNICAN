%Re-writing the question_1 code to finalise it here

mat1_uncleaned = readtable('SCPUnion2.1_AllSNe.txt', 'ReadVariableNames', false);


for col = 2:width(mat1_uncleaned)
    mat1_uncleaned{:, col} = str2double(string(mat1_uncleaned{:, col}));
end


if width(mat1_uncleaned) > 4
    mat1_uncleaned(:, end) = [];
end

% Filtering out rows with -99.9 values
rowsToRemove = any(mat1_uncleaned{:, 2:end} == -99.9, 2);
mat1 = mat1_uncleaned(~rowsToRemove, :);

mat1.Properties.VariableNames = ["Supernova", "Redshift", "Magnitude", "Magnitude_error"];


y = mat1.Magnitude + 19.3; 
x = mat1.Redshift;
y_err = mat1.Magnitude_error;


dark_orange = [0.8500 0.3250 0.0980]; 

% --- THEORETICAL CALCULATIONS ---
c = 299792.458;         % Speed of light in km/s
H0 = 67.4;              % Hubble constant (Planck 2020)


z_theory = linspace(min(x), max(x), 200);

% 1. Lambda-CDM (Planck 2020)
Om_A = 0.315;
OL_A = 0.685;
Ok_A = 1 - Om_A - OL_A;

dL_A = zeros(size(z_theory));
for i = 1:length(z_theory)
    integrand = @(zz) 1./sqrt(Om_A*(1+zz).^3 + Ok_A*(1+zz).^2 + OL_A);
    comoving_dist = (c/H0) * integral(integrand, 0, z_theory(i));
    dL_A(i) = (1+z_theory(i)) * comoving_dist;
end
mu_A = 5 * log10(dL_A) + 25;

% 2. Einstein-de Sitter (Om=1, OL=0)
dL_B = (2*c/H0) * (1+z_theory) .* (1 - 1./sqrt(1+z_theory));
mu_B = 5 * log10(dL_B) + 25;


% --- PLOTTING ---
figure('Color', 'w');
hold on;

% Plotting of Supernova Data
h1 = errorbar(x, y, y_err, 'Color', dark_orange, 'LineStyle', 'none', ...
    'Marker', '.', 'CapSize', 0, 'DisplayName', 'Union2.1 SNe Data');

% Plotting of Lambda-CDM Curve
h2 = plot(z_theory, mu_A, 'k-', 'LineWidth', 2, ...
    'DisplayName', '\LambdaCDM (Planck 2020)');

% Plotting of Einstein-de Sitter Curve
h3 = plot(z_theory, mu_B, 'b--', 'LineWidth', 2, ...
    'DisplayName', 'Einstein-de Sitter (\Omega_m=1)');


xlabel('Redshift (z)');
ylabel('Distance Modulus (\mu)');
title('Supernova Distance-Redshift Relation');
legend([h1, h2, h3], 'Location', 'southeast');
grid on;
hold off;
