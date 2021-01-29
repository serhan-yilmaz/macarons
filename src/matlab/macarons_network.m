function [ indicators, info ] = macarons_network( C, X, W, k, varargin)
    p = inputParser;
    p.CaseSensitive = false;
    p.KeepUnmatched = true;
    validIntegerScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','positive','integer'});
    addRequired(p, 'C', @isnumeric);
	addRequired(p, 'X', @(x) isnumeric(x) || islogical(x));
	addRequired(p, 'W', @(x) isnumeric(x) || islogical(x));
	addRequired(p, 'k', validIntegerScalar);
    addParameter(p, 'UseAdjustedR2', false, @islogical);
    addParameter(p, 'InitialActiveRegionSize', 1000, @isnumeric);
    addParameter(p, 'GrowthFactor', 2, @isnumeric);
    addParameter(p, 'Verbose', false, @islogical);
    parse(p, C, X, W, k, varargin{:});
    param = p.Results;
	param.NumberOfFeatures = k;
    
    C = reshape(C, [], 1);
    [~, sortIndices] = sort(-C, 'ascend');    
    param.SortIndices = sortIndices;

    [indicators, info] = run_macarons(X, C, W, ...
        param.NumberOfFeatures, param);
end

function [funcCalcR2, C_subset, W_subset, indicators, Cact] = initializeActiveSet(X, C, W, nActive, param)
    nAll = size(X, 2);
    validSNPs = param.SortIndices(1:min(nActive, nAll));
    X_subset = double(X(:, validSNPs));
    C_subset = C(validSNPs);
    W_subset = W(validSNPs, validSNPs);   
    if(param.UseAdjustedR2)
        n = size(X, 1);
        funcCalcR2 = @(ind1, ind2) 1 - ...
            (1 - corr(X_subset(:, ind1), X_subset(:, ind2)).^2) .* (n-1) ./ (n-2);
    else
        funcCalcR2 = @(ind1, ind2) corr(X_subset(:, ind1), X_subset(:, ind2)).^2;
    end
    indicators = zeros(nActive, 1, 'logical');
    Cact = -Inf;
    if(nActive < nAll)
        Cact = C(param.SortIndices(nActive+1));
    end
end

function [I, info] = run_macarons(X, C, W, k, param)
    timeStart = tic();

    nAll = size(X, 2);
    nActive = param.InitialActiveRegionSize;
    [funcCalcR2, C_subset, W_subset, indicators, Cact] = ...
        initializeActiveSet(X, C, W, nActive, param);
    
    Cx = C_subset;
    IG = ones(nActive, 1);
    
    nSelected = 0;
    selectedIndices = zeros(1, k);
    timeTotalGrow = 0;
    while(true)
        [c, iSNP] = max(Cx);
        
        if((c < Cact) && (nActive < nAll)) % Grow the set
            timeGrow = tic();
            nActivePrevious = nActive;
            indicatorsPrevious = indicators;
            IG_previous = IG;
%             Cx_previous = Cx;
            nActive = min(nAll, round(nActive * param.GrowthFactor));
            if(param.Verbose)
                fprintf('Growing, nSelected: %d, Size: %d, cAct: %.2e\n', nSelected, nActivePrevious, Cact);
            end
            [funcCalcR2, C_subset, W_subset, indicators, Cact] = ...
                initializeActiveSet(X, C, W, nActive, param);
            
            indicators(1:nActivePrevious) = indicatorsPrevious;
            IG = ones(nActive, 1);
            IG(1:nActivePrevious) = IG_previous;
            for iSelected = 1:nSelected
                iSNP = selectedIndices(iSelected);
                indices = W_subset(:, iSNP) & (~indicators);
                indices(1:nActivePrevious) = false;
                R2 = funcCalcR2(iSNP, indices)';
                IG(indices) = IG(indices) .* (1 - R2);
            end
            Cx = C_subset .* IG;
            Cx(indicators) = -Inf;
            if(param.Verbose)
                fprintf('[Completed] Growing, New Size: %d, Time Passed: %.1fs\n', nActive, toc(timeGrow));
            end
            timeTotalGrow = timeTotalGrow + toc(timeGrow);
            continue;
        end
        
        Cx(iSNP) = -Inf;
        indicators(iSNP) = true;
        selectedIndices(nSelected+1) = iSNP;
        nSelected = nSelected + 1;
        if(nSelected >= min(k, nAll))
           break; 
        end
        indices = W_subset(:, iSNP) & (~indicators);
        
        R2 = funcCalcR2(iSNP, indices)';
        IG(indices) = IG(indices) .* (1 - R2);
        Cx(indices) = C_subset(indices) .* IG(indices);
    end
    
    I = false(nAll, 1);
    I(param.SortIndices(1:nActive)) = indicators;
    
    info = struct();
    info.TimeTotal = toc(timeStart);
    info.TimeGrow = timeTotalGrow;
    info.TimeRun = info.TimeTotal - info.TimeGrow;
    
end













