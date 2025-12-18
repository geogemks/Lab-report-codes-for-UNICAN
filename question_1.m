mat1_uncleaned = readtable('SCPUnion2.1_AllSNe.txt', 'ReadVariableNames', false);


for col = 2:width(mat1_uncleaned)
    mat1_uncleaned{:, col} = str2double(string(mat1_uncleaned{:, col}));
end


mat1_uncleaned(:, end) = [];
rowsToRemove = any(mat1_uncleaned{:, 2:end} == -99.9, 2);

mat1 = mat1_uncleaned(~rowsToRemove, :);
mat1.Properties.VariableNames = ["Supernova", "Redshift", "Magnitude", "Magnitude_error"];


y = mat1.Magnitude + 19.3;
x = mat1.Redshift;
y_err = mat1.Magnitude_error;

%Question 1 Plotting%
errorbar(x, y, y_err,'Color', dark_orange, 'LineStyle', 'none')
xlabel('Redshift')
ylabel('Distance Modulus')
title('Supernova Distance-Redshift Relation')