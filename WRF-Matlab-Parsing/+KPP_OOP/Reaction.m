classdef Reaction
    %REACTION This class represents a single reaction in a KPP mechanism
    
    properties(SetAccess = private)
        reactants = struct();
        reactant_names = {};
        products = struct();
        product_names = {};
        species_names = {};
        rate_handle = @(T,C_M) 0;
        is_photolysis = false;
    end
    
    methods
        function obj = Reaction(kpp_mech_line)
            %REACTION Construct a Reaction instance from a single mechanism
            %line in the KPP mechanism.
            [obj.reactants, obj.products, obj.rate_handle, obj.is_photolysis] = obj.ParseLine(kpp_mech_line);
        end
        
        function n = NumReactants(obj)
            n = numel(fieldnames(obj.reactants));
        end
        
        function n = NumProducts(obj)
            n = numel(fieldnames(obj.products));
        end
        
        function bool = IsReactant(obj, reactant)
            %IsReactant(REACTANT) Returns true if REACTANT is a reactant in
            %this reaction. REACTANT must be a char array, cell arrays of
            %chars, or a string array. If given an array, returns a logical
            %array of the same size.
            bool = isfield(obj.reactants, reactant);
        end
        
        function bool = IsProduct(obj, product)
            %IsProduct(PRODUCT) Returns true if PRODUCT is a product in
            %this reaction. PRODUCT must be a char array, cell arrays of
            %chars, or a string array. If given an array, returns a logical
            %array of the same size.
            bool = isfield(obj.products, product);
        end
        
        function bool = IsSpecies(obj, specie)
            %IsSpecies(SPECIE) Returns true if SPECIE is a reactant or
            %product in this reaction. SPECIE must be a char array, cell
            %arrays of chars, or a string array. If given an array, returns
            %a logical array of the same size.
            bool = obj.IsReactant(specie) || obj.IsProduct(specie);
        end
        
        function names = get.reactant_names(obj)
            names = fieldnames(obj.reactants);
        end
        
        function names = get.product_names(obj)
            names = fieldnames(obj.products);
        end
        
        function names = get.species_names(obj)
            names = veccat(obj.reactant_names, obj.product_names);
        end
        
        function c = char(obj)
            % char() Returns a character representation of this reaction
            
            reactant_fmt_spec = repmat({'%.3f %s'},1,obj.NumReactants());
            product_fmt_spec = repmat({'%.3f %s'},1,obj.NumProducts());
            fmt_spec = [strjoin(reactant_fmt_spec, ' + '), ' -> ',...
                        strjoin(product_fmt_spec, ' + '), ' : ',...
                        '%s'];
            reactant_cell = obj.CellStruct(obj.reactants);
            product_cell = obj.CellStruct(obj.products);
            c = sprintf(fmt_spec, reactant_cell{:}, product_cell{:}, func2str(obj.rate_handle));
        end
    end
    
    methods(Static, Access = private)
        function [reactants, products, rate_fxn, is_photo] = ParseLine(kpp_mech_line)
            % Use the existing mechanism line parser
            [prod,prod_coeffs,react,react_coeffs,rate_value,is_photo] = read_wrf_mech_line(kpp_mech_line);
            % Create the product and reactant structures in a single call
            % to struct(), which requires that the species name and
            % coefficient alternate in a cell array.
            reactants = make_struct_from_field_values(react, react_coeffs);
            products = make_struct_from_field_values(prod, prod_coeffs);
            
            if is_photo
                rate_fxn = @tuv_wrapper;
            else
                rate_fxn = rate_value;
            end
            
            function j_value = tuv_wrapper(date_in, hour_in, lon, lat, is_utc, varargin)
                if isempty(varargin)
                    alt = 0.5;
                else
                    alt = varargin{1};
                end
                
                j_value = nan(size(lon));
                
                for i=1:numel(j_value)
                    j_value(i) = call_tuv(rate_value, date_in, hour_in, lon(i), lat(i), is_utc, alt);
                end
            end
        end
        
        
        
        function cell_out = CellStruct(struct_in)
            fns = fieldnames(struct_in);
            cell_out = cell(1,2*numel(fns));
            for i_fn = 1:numel(fns)
                i_cell = (i_fn-1)*2 + 1;
                cell_out{i_cell+1} = fns{i_fn};
                cell_out{i_cell} = struct_in.(fns{i_fn});
            end
        end
    end
end

