classdef heteronetwork
    
    properties (SetAccess = private)
        W = [];
        NodeClasses = [];
        NodeClassLabels = {};
        NodeLabels = {};
        NodeProperties = {};
        NodePropertyLabels = {};
        EdgeClasses = {};
        EdgeClassLabels = {};
    end
    
    
    methods
        function [obj] = heteronetwork()
            
        end
        
        function [obj] = addNodeClass(obj, n, classLabel)
            if(nargin < 3)
                p = length(obj.NodeClassLabels) + 1;
                classLabel = sprintf('Class%d', p);
            end
            nNode = length(obj.W);
            nClass = length(obj.NodeClassLabels);
            if(nClass == 0)
                obj.W = sparse(n, n);
            else
                obj.W = [obj.W           sparse(nNode, n); ...
                        sparse(n, nNode) sparse(n, n)];
            end
            if(any(strcmpi(obj.NodeClassLabels, classLabel)))
               error('All class labels have to be unique.'); 
            end
            obj.NodeClassLabels = [obj.NodeClassLabels; {classLabel}];
            classIndex = length(obj.NodeClassLabels);
            obj.NodeClasses = [obj.NodeClasses; ones(n, 1).*classIndex];
            obj.NodeLabels = [obj.NodeLabels; cell(n, 1)];
        end
        
        function [classLabels] = getNodeClasses(obj)
            classLabels = obj.NodeClassLabels(obj.NodeClasses); 
        end
        
        function [indices] = isNodeClass(obj, classLabel)
            indices = obj.mapNodes(classLabel); 
        end
        
        function [Wsub] = subnetwork(obj, classLabelA, classLabelB)
            indicesA = obj.mapNodes(classLabelA);
            indicesB = obj.mapNodes(classLabelB);
            Wsub = obj.W(indicesA, indicesB);
        end
        
        function [obj] = createNodeProperty(obj, propertyLabel, propertyType, defaultValues)
            if(nargin < 3)
               propertyType = 'cell'; 
            end
            nNode = length(obj.W);
            if(nargin < 4)
                switch(propertyType)
                    case 'logical'
                        defaultValues = false(nNode, 1);
                    case 'numeric'
                        defaultValues = zeros(nNode, 1);
                    case 'cell'
                        defaultValues = cell(nNode, 1);
                    otherwise
                        error('Invalid node property type.');
                end 
            end
            if(isscalar(defaultValues))
                defaultValues = repmat(defaultValues, nNode, 1);
            end
            
            obj.NodeProperties = [obj.NodeProperties, {defaultValues}];
            obj.NodePropertyLabels = [obj.NodePropertyLabels, {propertyLabel}];
        end
        
        function [obj] = setNodeProperty(obj, propertyLabel, values, classLabel)
            [nodePropertyIndex] = mapNodeProperty(obj, propertyLabel, false);
            if(nodePropertyIndex == 0)
               propertyType = class(values);
               if(isnumeric(values))
                   propertyType = 'numeric';
               end
                obj = createNodeProperty(obj, propertyLabel, propertyType);
                nodePropertyIndex = length(obj.NodeProperties);
            end
            if(nargin < 4)
                nNode = length(obj.W);
                if(nNode ~= length(values))
                    error('Number of node labels must be equal to the number of nodes');
                end
                obj.NodeProperties{nodePropertyIndex} = values;
            else
                indices = obj.mapNodes(classLabel);
                if(nnz(indices) ~= length(values))
                    error('Number of node labels must be equal to the number of nodes in the specified class.');
                end
                obj.NodeProperties{nodePropertyIndex}(indices) = values;
            end
        end
        
        function [values] = getNodeProperty(obj, propertyLabel, classLabel)
            [nodePropertyIndex] = mapNodeProperty(obj, propertyLabel);
            if(nargin < 4)
                values = obj.NodeProperties{nodePropertyIndex};
            else
                indices = obj.mapNodes(classLabel);
                values = obj.NodeProperties{nodePropertyIndex}(indices);
            end
        end
        
        function [indices] = getEdgeClass(obj, classLabel, subscriptFlag)
            if(nargin < 3)
               subscriptFlag = false; 
            end
            [edgeClassIndex] = mapEdgeClass(obj, classLabel);
            indices = obj.EdgeClasses{edgeClassIndex};
            if(~subscriptFlag)
                indices = sub2ind(size(obj.W), indices(:, 1), indices(:, 2));
            end
            
        end
        
        function [obj] = addNodeLabels(obj, labels, classLabel)
            if(isstring(labels))
               labels = cellstr(labels); 
            end
            if(~iscell(labels))
                error('Labels must be a string or a cell array of characters.');
            end
            if(nargin < 3)
                if(length(obj.NodeLabels) ~= length(labels))
                    error('Number of node labels must be equal to the number of nodes');
                end
                obj.NodeLabels = labels;
            else
                indices = obj.mapNodes(classLabel);
                if(nnz(indices) ~= length(labels))
                    error('Number of node labels must be equal to the number of nodes in the specified class.');
                end
                obj.NodeLabels(indices) = labels;
            end
        end
        
        function [obj] = addEdges(obj, Wa2b, classLabelA, classLabelB, edgeClassLabel)
            if(nargin < 4)
                classLabelB = classLabelA;
            end
            if(nargin < 5)
               edgeClassLabel = sprintf('%s2%s', classLabelA, classLabelB); 
            end
            indicesA = obj.mapNodes(classLabelA);
            indicesB = obj.mapNodes(classLabelB);
            
            edges1 = find(obj.W(indicesA, indicesB));
%             edges2 = find(Wa2b);
%             if(~isempty(intersect(edges1, edges2)))
%                 error('Overlapping edges are not supported.');
%             end
            if(islogical(Wa2b))
                obj.W(indicesA, indicesB) = obj.W(indicesA, indicesB) | Wa2b;
                obj.W(indicesB, indicesA) = obj.W(indicesB, indicesA) | Wa2b';
            else
%                 edges1 = find(obj.W(indicesA, indicesB));
                edges2 = find(Wa2b);
                if(~isempty(intersect(edges1, edges2)))
                    error('Overlapping edges are not supported for weighted networks.');
                end
                if(strcmpi(classLabelA, classLabelB))
                    Ws = 0.5*(Wa2b + Wa2b');
                    obj.W(indicesA, indicesB) = Ws;
                else
                    obj.W(indicesA, indicesB) = obj.W(indicesA, indicesB) + Wa2b;
                    obj.W(indicesB, indicesA) = obj.W(indicesB, indicesA) + Wa2b';
                end
%                 obj.W(indicesA, indicesB) = Wa2b;
%                 obj.W(indicesB, indicesA) = Wa2b';
            end
            edges2 = find(obj.W(indicesA, indicesB));
            addedEdges = setdiff(edges2, edges1);
            [i1, i2] = ind2sub(size(Wa2b), addedEdges);
            indfA = find(indicesA);
            indfB = find(indicesB);
%             indices = sub2ind(size(obj.W), indfA(i1), indfB(i2));
            obj = addEdgeClass(obj, indfA(i1), indfB(i2), edgeClassLabel);
        end
        

        
    end
    
    methods (Access = private)
        function [nodeClassIndex] = mapNodeClass(obj, classLabel)
            [isMapped, nodeClassIndex] = ...
                ismember(classLabel, obj.NodeClassLabels);
            if(~isMapped)
               error('Node class with the specified class label is not found.'); 
            end
        end
        
        function [edgeClassIndex] = mapEdgeClass(obj, classLabel)
            [isMapped, edgeClassIndex] = ...
                ismember(classLabel, obj.EdgeClassLabels);
            if(~isMapped)
               error('Edge class with the specified class label is not found.'); 
            end
        end
        
        function [nodePropertyIndex] = mapNodeProperty(obj, propertyLabel, errorFlag)
            if(nargin < 3)
               errorFlag = true;
            end
            [isMapped, nodePropertyIndex] = ...
                ismember(propertyLabel, obj.NodePropertyLabels);
            if(~isMapped && errorFlag)
                error('Node property with the specified label is not found.'); 
            end
        end
        
        
        
        function [indices, nodeClassIndex] = mapNodes(obj, classLabel)
            [nodeClassIndex] = mapNodeClass(obj, classLabel);
            indices = obj.NodeClasses == nodeClassIndex;
        end
        
        function [obj] = addEdgeClass(obj, indices1, indices2, classLabel)
            if(any(strcmpi(obj.EdgeClassLabels, classLabel)))
               error('All edge class labels have to be unique.'); 
            end
            obj.EdgeClassLabels = [obj.EdgeClassLabels; {classLabel}];
%             classIndex = length(obj.EdgeClassLabels);
            indices1 = reshape(indices1, [], 1);
            indices2 = reshape(indices2, [], 1);
            obj.EdgeClasses = [obj.EdgeClasses; {[indices1, indices2]}];
        end
        
        
    end
end



