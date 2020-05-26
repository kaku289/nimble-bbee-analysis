classdef Interval < handle
    properties (Access = public)
        
        leftEnd = nan;
        rightEnd = nan;
        size = nan; % # of points in the interval
        type = {'closed', 'closed'}; % By default closed on both sides
        unfitness = nan; % measure of how fit an interval is e.g., RMSE error for a constant fit within the interval
        
    end
     methods
         function obj = Interval(a,b)
             if nargin ~= 0
                 assert(a>0 && b>a); 
                 obj.leftEnd = a;
                 obj.rightEnd = b;
                 obj.size = b-a+1;
             end
         end
     end
     methods(Static)
        function intervalArray = createIntervalArray(input, unfitness_values)
            % a - N by 2 matrix containing end points of intervals
            % intervalArray - N by 1 interval array 
            N = size(input,1);
            assert(size(input,2) == 2 && sum(input(:,1)>0) == N && sum(input(:,2)>input(:,1)) == N && numel(unfitness_values)==N);
            intervalArray(size(input,1),1) = Interval();
            
            leftendpoints = num2cell(input(:,1));
            [intervalArray.leftEnd] = leftendpoints{:};
            
            rightendpoints = num2cell(input(:,2));
            [intervalArray.rightEnd] = rightendpoints{:};
            
            s = num2cell(input(:,2)-input(:,1) + 1);
            [intervalArray.size] = s{:};
            
            unfitness = num2cell(unfitness_values);
            [intervalArray.unfitness] = unfitness{:};
        end
        
        function neighbors = findOverlappingIntervals(intervalArray, indx)
            % intervalArray - N by 1 interval array
            % indx - Index of an interval whose neighbors (the ones it
            % overlaps with) are to be found
            % neighbors - N by 1 binary array containing true at indices
            % where intervals are neighbors of intervalArray(indx)
            
            a = intervalArray(indx);
            neighbors = ~([intervalArray.rightEnd] <= a.leftEnd | [intervalArray.leftEnd] >= a.rightEnd);
            neighbors(indx) = false;
            
        end
        
        function visited = findRepresentativeIntervals(intervalArray)
            % intervalArray - N by 1 interval array
            % output = N by 1 interval array containing true for intervals
            % that do not directly or indirectly overlap with any other
            % The best among the ones that overlap are selected using the
            % below algorithm.
            
            
            % Algorithm
            % 1. Find the most fit interval (the one with lowest value of rmse_rdot)
            % 2. Prune all intervals that it is connected to
            % 3. Go to the next "most fit" interval among the remaining
            % intervals
            N = numel(intervalArray);
            
            prune = false(N,1);
            visited = false(N,1);
            
            ct = 1;
            while sum(prune | visited) ~= N
                intervalsLeft = intervalArray(~prune & ~visited);
                [~, indx] = min([intervalsLeft.unfitness]);
                indx = find(intervalArray == intervalsLeft(indx));
                visited(indx) = true;
                
                neighbors = Interval.findOverlappingIntervals(intervalArray, indx); % It contains all neighbors irrespective of the fa they are visited/pruned
                prune(neighbors) = true;          % neighbors & ~visited?    
                
                
                if ct > N
                   error('Number of while loop evaluations exceeded number of intervals. Should NOT happen!! Check Algorithm if it happens');
                end
                ct = ct+1;
            end
            
        end
        
        function visited = findRepresentativeIntervalsPerTimeWindow(intervalArray)
            % intervalArray - N by 1 interval array
            % output = N by 1 boolean array containing true for intervals
            % that do not directly or indirectly overlap with any other
            % The best among the ones that overlap are selected using the
            % below algorithm.
            
            
            % Algorithm - FOR EACH TIME WINDOW
            % 1. Find the most fit interval (the one with lowest value of unfitness score)
            % 2. Prune all intervals that it is connected to
            % 3. Go to the next "most fit" interval among the remaining
            % intervals
            N = numel(intervalArray);
            
            visited = false(N,1);
            
            s = [intervalArray.size]; % # of data points in each interval
            unique_s = unique(s);
            for ct=1:length(unique_s)
                current_s = unique_s(ct);
                intervalArray_subset = intervalArray(s==current_s);
                if length(intervalArray_subset)==1
                    visited(s==current_s) = true;
                else
                    N_sub = length(intervalArray_subset);
                    prune_sub = false(N_sub,1);
                    visited_sub = false(N_sub,1);
                    
                    ct1 = 1;
                    while sum(prune_sub | visited_sub) ~= N_sub
                        intervalsLeft = intervalArray_subset(~prune_sub & ~visited_sub);
                        [~, indx] = min([intervalsLeft.unfitness]);
                        indx = find(intervalArray_subset == intervalsLeft(indx));
                        visited_sub(indx) = true;
                        
                        neighbors = Interval.findOverlappingIntervals(intervalArray_subset, indx); % It contains all neighbors irrespective of the fa they are visited/pruned
                        prune_sub(neighbors) = true;          % neighbors & ~visited?
                        
                        if ct1 > N_sub
                            error('Number of while loop evaluations exceeded number of intervals. Should NOT happen!! Check Algorithm if it happens');
                        end
                        ct1 = ct1+1;
                    end
                    
                    indices = arrayfun(@(x) find(intervalArray==x) ,intervalArray_subset(visited_sub));
                    visited(indices) = true;
                end
            end
        end
        
        function output = findOverlappingChain(intervalArray)
            % intervalArray - N by 1 interval array
            % output = N by 1 containing values from 1 onwards
            % Intervals that directly or indirectly overlap get the same
            % number in output
            
            % Write it in a graph form and perform BFS or DFS
            
            
            
        end
    end
end