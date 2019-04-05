classdef RunningStdDev < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        weights = [];
        weights_squared = [];
        mean = [];
        variance = [];
    end
    
    properties(SetAccess = protected)
        was_init = false;
    end
    
    methods
        function obj = RunningStdDev(varargin)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
                obj.setupArrays(varargin{1});
            end
        end
        
        function addData(obj, values, weights)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            if nargin < 3
                weights = ones(size(values));
            end
            
            % By default, give zero weight to NaN values as they should not
            % really contribute to the average.
            weights(isnan(values)) = 0;
            
            if ~obj.was_init
                obj.setupArrays(size(values));
            else
                obj.checkArrays(size(values));
            end
            
            obj.weights = obj.weights + weights;
            obj.weights_squared = obj.weights_squared + weights .^ 2;
            old_mean = obj.mean;
            obj.mean = old_mean + (weights ./ obj.weights) .* (values - old_mean);
            obj.variance = obj.variance + weights .* (values - old_mean) .* (values - obj.mean);
        end
        
        function sigma = getReliabilityStdDev(obj)
            if isequal(obj.weights, obj.weights_squared)
                warning('RunningStdDev:unity_weights', 'It looks like you may not have provided weights for the running std. dev. In that case, the getFrequencyStdDev method may be better');
            end
            sigma = sqrt(obj.variance ./ (obj.weights - obj.weights_squared ./ obj.weights));
        end
        
        function sigma = getFrequencyStdDev(obj)
            sigma = sqrt(obj.variance ./ (obj.weights - 1));
        end
    end
    
    methods(Access = protected)
        function setupArrays(obj, shape)
            obj.weights = zeros(shape);
            obj.weights_squared = zeros(shape);
            obj.mean = zeros(shape);
            obj.variance = zeros(shape);
            obj.was_init = true;
        end
        
        function checkArrays(obj, shape)
            curr_size = size(obj.variance);
            if ~isequal(shape, curr_size)
                error('RunningStdDev:inconsistent_shape', 'The input values have a different size (%s) than the existing values (%s)', mat2str(shape), mat2str(curr_size));
            end
        end
    end
end

