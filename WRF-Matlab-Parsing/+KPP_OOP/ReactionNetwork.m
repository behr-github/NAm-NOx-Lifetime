classdef ReactionNetwork < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        reactions;
    end
    
    methods
        function obj = ReactionNetwork(mechanism_file)
            %ReactionNetwork Create an instance of a reaction network that
            %represents a full numerical chemical mechanism.
            %   ReactionNetwork(mechanism_file) Parses the KPP mechanism
            %   given by MECHANISM_FILE, which can either be a path to a
            %   .eqn file or the name of a mechanism stored in the
            %   WRF_Parsing directory.
            
            % Allow the "mechanism file" to be either a path to a file or a
            % mechanism name.
            % At least in R2017b, exist( ___, 'file' ) will return 7 if the
            % "file" is actually a directory and 2 if it is an actual file.
            if exist(mechanism_file, 'file') ~= 2
                mechanism_file = fullfile(obj.MyDir(), '..', upper(mechanism_file), sprintf('%s.eqn', mechanism_file));
                if ~exist(mechanism_file, 'file')
                    error('kpp:unknown_mechanism', '%s is not a file nor a recognized mechanism', mechanism_file);
                end
            end
            obj.reactions = obj.ParseMechFile(mechanism_file);
        end
        
        function [rxns, rxn_inds] = FindReactionsWithReactants(obj, reactants)
            rxn_inds = obj.FindReactionsBySpecies(reactants, @IsReactant);
            rxns = obj.reactions(rxn_inds);
        end
        
        function [rxns, rxn_inds] = FindReactionsWithProducts(obj, products)
            rxn_inds = obj.FindReactionsBySpecies(products, @IsProduct);
            rxns = obj.reactions(rxn_inds);
        end
        
        function [rxns, rxn_inds] = FindReactionsWithSpecies(obj, species)
            rxn_inds = obj.FindReactionsBySpecies(species, @IsSpecie);
            rxns = obj.reactions(rxn_inds);
        end
        
        function [production_rxns, loss_rxns] = FindSteadyStateReactions(obj, ss_specie)
            production_inds = false(size(obj.reactions));
            loss_inds = false(size(obj.reactions));
            for i_rxn = 1:numel(obj.reactions)
                production_inds(i_rxn) = IsProduct(obj.reactions{i_rxn}, ss_specie);
                loss_inds(i_rxn) = IsReactant(obj.reactions{i_rxn}, ss_specie);
            end
            production_rxns = obj.reactions(production_inds);
            loss_rxns = obj.reactions(loss_inds);
        end
        
        function conc = CalculateSpeciesSteadyState(obj, ss_species, wrf_file, varargin)
            % CalculateSpeciesSteadyState(obj, SS_SPECIE, WRF_FILE)
            %   Calculates the steady state concentration of SS_SPECIE
            %   given the other species concentrations in WRF_FILE. There
            %   are several parameters:
            %
            %       'n_iter' - how many times to iterate through the
            %       solution. Default is 1, but if you are trying to solve
            %       for multiple steady-state species that depend on each
            %       other, increasing this *may* improve your result, since
            %       successive times through will use concentrations of
            %       steady state species from the previous iteration,
            %       instead of just skipping their contribution to
            %       production/loss.
            %
            %       'return_unit' - what unit to return the data in.
            %       Default is 'ppm', i.e. parts-per-million. This can be
            %       set to any unit that CONVERT_UNITS() knows how to
            %       convert to from 'ppp' (parts-per-part i.e. unscaled
            %       mixing ration) or 'ndens'/'number density', which will 
            %       return it in molec. cm^{-3}.
            %
            %       'start_count_stride' - a cell array containing the
            %       start, count, and stride arguments to pass to NCREAD
            %       when loading the necessary WRF-Chem species.
            
            % We'll calculate the steady state by first finding all the
            % production and loss reactions for that species. At steady
            % state, the production and loss rates must be equal, so for
            % species S, we'd have:
            %
            %   k_l1[C_1][S] + ... + k_ln[C_n][S] =
            %       k_p1[A_1][B_1] + ... k_pn[A_n][B_n]
            %
            % therefore [S] will be the sum of the rates of the production
            % reactions divided by the sum of the rate constants times
            % other species in the loss reactions.
            %
            % What we need to do is loop over both sets of reactions,
            % calculate their rate constants, and then calculate the
            % necessary sums.
            
            E = JLLErrors;
            p = advInputParser;
            p.addParameter('n_iter', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
            p.addParameter('start_count_stride', {});
            p.addParameter('return_unit', 'ppm');
            p.parse(varargin{:});
            pout = p.Results;
            
            n_iter = pout.n_iter;
            return_unit = pout.return_unit;
            start_count_stride = pout.start_count_stride;
            if isempty(start_count_stride)
                start = ones(1,4);
                count = inf(1,4);
                stride = ones(1,4);
            else
                start = start_count_stride{1};
                count = start_count_stride{2};
                stride = start_count_stride{3};
            end
            
            if ischar(ss_species)
                ss_species = {ss_species};
            elseif ~iscellstr(ss_species)
                E.badinput('SS_SPECIES must be a char array or cell array of chars');
            end
            
            % We'll need temperature and number density of air for the
            % non-photolysis reaction rates
            wrf_t = read_wrf_preproc(wrf_file, 'temperature');
            wrf_ndens = read_wrf_preproc(wrf_file, 'number density');
            wrf_size = size(wrf_t);
            % We'll need the lat/lon and date/time for photolysis rates
            wrf_lon = ncread(wrf_file, 'XLONG');
            wrf_lat = ncread(wrf_file, 'XLAT');
            wrf_datetime = date_from_wrf_filenames(wrf_file);
            
            % Prep a structure to cache species that we've already read in
            wrf_conc = struct();
            for i_iter = 1:n_iter
                % prep the output
                conc = make_empty_struct_from_cell(ss_species);
                
                for i_ss = 1:numel(ss_species)
                    [prod_rxn, loss_rxn] = obj.FindSteadyStateReactions(ss_species{i_ss});
                    total_prod = calculate_ss_sum(prod_rxn, ss_species(i_ss), wrf_conc);
                    total_loss = calculate_ss_sum(loss_rxn, ss_species(i_ss), wrf_conc);
                    conc.(ss_species{i_ss}) = total_prod ./ total_loss;
                end
                
                if i_iter < n_iter
                    wrf_conc = copy_structure_fields(conc, wrf_conc);
                end
            end
            
            if ~any(strcmpi(return_unit, {'ndens', 'number density'}))
                for i_ss = 1:numel(ss_species)
                    conc.(ss_species{i_ss}) = conc.(ss_species{i_ss}) ./ wrf_ndens;
                    conc.(ss_species{i_ss}) = convert_units(conc.(ss_species{i_ss}), 'ppp', return_unit);
                end
            end
            
            function total_rate = calculate_ss_sum(rxns, species_to_omit, wrf_conc)
                total_rate = zeros(wrf_size);
                for i_rxn = 1:numel(rxns)
                    this_rxn = rxns{i_rxn};
                    if ~this_rxn.is_photolysis
                        rate = this_rxn.rate_handle(wrf_t, wrf_ndens);
                    else
                        rate = KPP_OOP.interp_tuv_photolysis(this_rxn.rate_handle, wrf_datetime, wrf_lon, wrf_lat);
                    end
                    
                    
                    species_to_load = this_rxn.reactant_names(~ismember(this_rxn.reactant_names, species_to_omit));
                    
                    for i_spec = 1:numel(species_to_load)
                        this_species = species_to_load{i_spec};
                        try
                            wrf_conc = load_species(this_species, wrf_conc);
                        catch err
                            if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:unknownLocation')
                                fprintf('Problem calculating SS for %s:\n  "%s"\n  Skipping', char(this_rxn), err.message);
                                continue
                            else
                                rethrow(err)
                            end
                        end
                        
                        rate = rate .* wrf_conc.(this_species);
                    end
                    total_rate = total_rate + rate;
                end
            end
            
            function wrf_conc = load_species(species, wrf_conc)
                % We don't want to read species concentrations over and
                % over because that can be time consuming, so we'll cache
                % them in the wrf_conc structure and only load them if we
                % need to.
                if ~isfield(wrf_conc, species)
                    try
                        val = ncread(wrf_file, lower(species), start, count, stride);
                        unit = ncreadatt(wrf_file, lower(species), 'units');
                    catch err
                        if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:unknownLocation')
                            val = ncread(wrf_file, species, start, count, stride);
                            unit = ncreadatt(wrf_file, species, 'units');
                        else
                            rethrow(err)
                        end
                    end
                    
                    % We need to convert the concentrations in
                    % molecules/cm^3. First convert to unscaled mixing
                    % ratio then multiple by number density of air.
                    val = wrf_ndens .* convert_units(val, unit, 'ppp');
                    
                    wrf_conc.(species) = val;
                end
            end
        end
    end
    
    methods(Access = private)
        function rxn_inds = FindReactionsBySpecies(obj, species, test_handle)
            rxn_inds = false(size(obj.reactions));
            for i_rxn = 1:numel(rxn_inds)
                rxn_inds(i_rxn) = any(test_handle(obj.reactions{i_rxn}, species));
            end
        end
    end
    
    methods(Static, Access = private)
        function mdir = MyDir()
            mdir = fileparts(mfilename('fullpath'));
        end
        
        function reactions = ParseMechFile(mech_file)
            % This should eventually also parse the species file, if only
            % to make sure that no erroneous species are expected in the
            % equation file.
            reactions = cell(1000,1); % will clean up extra entries at the end;
            fid = fopen(mech_file);
            try
                tline = fgetl(fid);
                i_rxn = 1;
                while ischar(tline)
                    if strcmp(tline(1),'#')
                        % Skip comment lines
                        tline = fgetl(fid);
                        continue
                    end
                    fprintf('Parsing line: %s\n', tline);
                    reactions{i_rxn} = KPP_OOP.Reaction(tline);
                    i_rxn = i_rxn + 1;
                    tline = fgetl(fid);
                end
            catch err
                fclose(fid);
                rethrow(err);
            end
            fclose(fid);
            reactions(i_rxn:end) = [];
        end
        
    end
end

