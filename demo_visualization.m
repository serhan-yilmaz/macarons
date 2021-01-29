addpath(genpath([cd(), '/src/matlab/']));
rng(1, 'twister'); % For reproducibility
%% Prepare Data
phenotypeIndex = 1;     % Index of the phenotype - 4W
load('data/data.mat');

% Valid indices for the selected phenotype
sampleIndices = ~isnan(Y(:, phenotypeIndex));
% X : genotype  - n x m  matrix - n : Samples, m : Features (SNPs)
X = double(X(sampleIndices, :)); 
% Y : phenotype - n x 1 column vector
Y = Y(sampleIndices, phenotypeIndex);
% S: positions and chromosomes of the SNPs - n x 2 matrix
% S(:, 1) = chromosome indices
% S(:, 2) = Position on the chromosome
S  = double(snp);
% C : scores    - m x 1 column vector
C = corr(X, Y).^2; % Use squared correlation

%% Visualize the baseline results without taking dependencies into account
k    = 100;             % Number of Features to be selected
D    = 0;               % Intra-chromosomal distance in base pairs
                        % to limit the search space of Macarons
          
[I] = macarons(C, X, S, k, D);

figure(1);
clf();
set(gcf, 'Position', [0 0 680 600]);
movegui('center');
set(gcf, 'Color', [1 1 1]);
visualizevariants(C, X, S, I);
h = suptitle('Baseline (D=0)');
h.Position(2) = h.Position(2) - 0.035;

%% Visualize the results of Macarons for a moderate D parameter
k    = 100;             % Number of Features to be selected
D    = 2e4;             % Intra-chromosomal distance in base pairs
                        % to limit the search space of Macarons
          
% The first output of Macarons (I) is an m x 1 logical column vector.
% I contains the indicators for a subset of features
% i.e. feature i is selected iff S[i] = true
[I] = macarons(C, X, S, k, D);

figure(2);
clf();
set(gcf, 'Position', [0 0 680 600]);
movegui('center');
set(gcf, 'Color', [1 1 1]);
visualizevariants(C, X, S, I);
h = suptitle('Macarons (D=20 kbp)');
h.Position(2) = h.Position(2) - 0.035;

%% Visualize the results of Macarons for a high D parameter
k    = 100;             % Number of Features to be selected
D    = 1e6;             % Intra-chromosomal distance in base pairs
                        % to limit the search space of Macarons
          
[I] = macarons(C, X, S, k, D);

figure(3);
clf();
set(gcf, 'Position', [0 0 680 600]);
movegui('center');
set(gcf, 'Color', [1 1 1]);
visualizevariants(C, X, S, I);
h = suptitle('Macarons (D=1e6)');
h.Position(2) = h.Position(2) - 0.035;




