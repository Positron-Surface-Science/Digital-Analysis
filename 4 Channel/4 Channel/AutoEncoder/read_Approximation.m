function [x t Nv] = read_Approximation(training_file, N, M)
%[x t Nv] = read_approx_file(training_file, N, M)
% The following code reads a text file and stores all the paterns in 
% an Nv by (N+M) matrix
fid = fopen(training_file, 'r');
if(fid == -1)
    fprintf('Could not open file %s\n', training_file)
end
training_file_values = fscanf(fid, '%f');
fclose(fid);

NCOLS = N+M;
Nv = numel(training_file_values)/NCOLS;
fprintf('# of patterns in %s = %d\n', training_file, Nv);
training_file_values = reshape(training_file_values, [NCOLS Nv])';

% Store the inputs in variable x and the currect class ID in ic
x = training_file_values(:, 1:N);
t = training_file_values(:, N+1:end);
