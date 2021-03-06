addpath(genpath([cd(), '/src/matlab/']));
%% Prepare Data
phenotypeIndex = 1;     % Index of the phenotype - 4W
nPCs           = 1;     % Number of Principal Components used
                        % to correct population stratification

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
C = computeSKAT(X, Y, 'k', nPCs);

%% Run Macarons with intra-chromosomal distance (D) parameter
k    = 100;             % Number of Features to be selected
D    = 2e4;             % Intra-chromosomal distance in base pairs
                        % to limit the search space of Macarons
          
% The first output of Macarons (I) is an m x 1 logical column vector.
% I contains the indicators for a subset of features
% i.e. feature i is selected iff S[i] = true
[I] = macarons(C, X, S, k, D);


%% Run Macarons using a dependency network
load('data/network.mat');

% W : dependency network - n x n logical sparse matrix
W = logical(W);
% k : Number of Features to be selected
k = 100;             
          
                        
% Instead of using an intra-chromosomal distance, use a dependency network
% to limit the search space for macarons:
[I2] = macarons_network(C, X, W, k);


