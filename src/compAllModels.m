% Computes all models successively
% NB: Each model can be run separately
% runMeFirst.sh has to be executed to create load_paths.m

% First we load the variables of interest, namely the motif of interest, an
% initial PWM for this motif and a set of N-mers to be searched (here, all
% N-mers from a ChIP study)
fprintf('Loading main variables...\n')
load_params;

% PWM model
if (~exist(filePWM,'file'))
    fprintf('Computing Position Weight Matrix (PWM) model...\n')
    PWM_model;
else
    fprintf('Skipping PWM model, result file already created...\n')
end

% MIX
if (~exist(fileMIX,'file'))
    fprintf('Computing mixture model (MIX)...\n')
    MIX_model;
else
    fprintf('Skipping mixture model, result file already created...\n')
end

% PIM
if (~exist(filePIM,'file'))
    fprintf('Computing Pairwise Interaction Model (PIM)...\n')
    nearest_neighbors = false;
    PIM_model;
else
    fprintf('Skipping PIM, result file already created...\n')
end

% NNM
if (~exist(fileNNM,'file'))
    fprintf('Computing Nearest-Neighbors Model (NNM)...\n')
    nearest_neighbors = true;
    NNM_model;
else
    fprintf('Skipping NNM, result file already created...\n')
end