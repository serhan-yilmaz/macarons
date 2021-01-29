function [ output_args ] = visualizevariants( C, X, S, I)
    nSNP = nnz(I);
    Xs = X(:, I);
    R = corr(Xs).^2;
    Wsnp2snp = R >= 0.5;

    Sx = S(I, :);
    Cx = C(I); %/ max(C);
    [chromosomes, ~, chr_indices] = unique(Sx(:, 1));
    nChr = length(chromosomes);
    Wsnp2chr = sparse((1:nSNP)', chr_indices, true, nSNP, nChr);

    chromosome_names = cell(nChr, 1);
    for iChr = 1:nChr
        chr = chromosomes(iChr);
        chromosome_names{iChr} = sprintf('Chr %d', chr);
    end

    meanVal = mean(Cx);

    hnet = heteronetwork();
    hnet = addNodeClass(hnet, nSNP, 'SNP');
    hnet = addNodeClass(hnet, nChr, 'Chromosome');
    hnet = addEdges(hnet, R>=0.7, 'SNP', 'SNP', 'HighlyRedundantEdges');
    hnet = addEdges(hnet, R>=0.6, 'SNP', 'SNP', 'ModeratelyRedundantEdges');
    hnet = addEdges(hnet, Wsnp2snp, 'SNP', 'SNP');
    hnet = addEdges(hnet, Wsnp2chr, 'SNP', 'Chromosome');
    hnet = hnet.setNodeProperty('NodeScaling', ones(nChr, 1)*meanVal*2.5, 'Chromosome');
    hnet = hnet.setNodeProperty('NodeScaling', Cx, 'SNP');
    hnet = hnet.setNodeProperty('NodeLabel', chromosome_names, 'Chromosome');

    net = networkvisualizer(hnet.W);

    edges = net.Edges;
    hEdges = ismember(edges, hnet.getEdgeClass('HighlyRedundantEdges'));
    mEdges = ismember(edges, hnet.getEdgeClass('ModeratelyRedundantEdges'));

    colorBlue = [0 0.447 0.741];
    colorOrange = [0.929 0.694 0.125];

    net = net.addNodeClass(hnet.getNodeClasses(), 'nodeType');
    net = net.scaleNodeSizes(hnet.getNodeProperty('NodeScaling'), 0.6);
    net = net.setNodeLabels(hnet.getNodeProperty('NodeLabel'));
    net = setNodeCurvature(net, 0.2, {'Chromosome'}, 'nodeType');
    net = setNodeColors(net, colorBlue, 'SNP', 'nodeType');
    net = setNodeColors(net, colorOrange, 'Chromosome', 'nodeType');
    net = setNodeFontSize(net, 11);
    net = setNodeSizes(net, 'auto');
    net = setEdgeLineWidth(net, 0.75);
    net = addEdgeClass(net, hEdges, 'HighlyRedundantEdges');
    net = addEdgeClass(net, mEdges, 'ModeratelyRedundantEdges');
    net = createEdgeClass(net, 'Chr2SNP', 'Chromosome', 'SNP', 'nodeType');
    net = setEdgeColors(net, 0.75*[1 1 1], true, 'Chr2SNP');
    net = setEdgeLineWidth(net, 0.5, true, 'Chr2SNP');
    net = setEdgeColors(net, 0.15*[1 1 1], true, 'HighlyRedundantEdges');
    net = setEdgeLineWidth(net, 1.75, true, 'HighlyRedundantEdges');
    net = setEdgeColors(net, 0.4*[1 1 1], true, 'HighlyRedundantEdges');
    net = setEdgeLineWidth(net, 1.25, true, 'HighlyRedundantEdges');
    
    net = dolayout(net);
    plot(net);
end

