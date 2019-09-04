function write_fits_to_ncdf(nc_dir, fit_type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    fit_type = 'behr';
end

switch lower(fit_type)
    case 'behr'
        fit_file_name = @misc_emissions_analysis.behr_fit_file_name;
        dow = {'TWRF', 'US'};
        suffix = '';
    case 'wrf'
        fit_file_name = @(win, dow) misc_emissions_analysis.wrf_fit_file_name(win, 'NO2_VCDS', dow);
        dow = {'TWRF'};
        suffix = '_wrf';
    case 'nasa'
        fit_file_name = @misc_emissions_analysis.nasa_fit_file_name;
        dow = {'TWRF', 'US'};
        suffix = '_nasa';
end

year_windows = arrayfun(@(y) (y-1):(y+1), 2006:2013, 'uniform', false);
for i_yr = 1:numel(year_windows)
    win = year_windows{i_yr};
    for i_dow = 1:numel(dow)
        fit_file = fit_file_name(win, dow{i_dow});
        fit = load(fit_file);
        nc_file_name = sprintf(fullfile(nc_dir, 'EMG_fits_%d-%d_%s%s.nc'), win(1), win(end), dow{i_dow}, suffix);
        delete(nc_file_name);
        write_one_fit(fit.locs, nc_file_name);
    end
end

end

function write_one_fit(locs, nc_file_name)

city_inds = misc_emissions_analysis.get_loc_inds_of_type('Cities');
locs = locs(city_inds);
n_ld = numel(locs(23).no2_sectors.x); % use LA, will definitely have the largest line densities
dims = struct('cities', numel(locs), 'nchars', 256, 'days', size(locs(1).WindDir,1),...
    'orbits', size(locs(1).WindDir,2), 'line_density', n_ld, 'fitting_parameters', 5);

atts_to_write = {'units', 'description'};

variables = struct('Location', struct('dims', {{'cities', 'nchars'}}, 'type', 'char'),...
    'Latitude', struct('dims', {{'cities'}}, 'units', 'degrees (south is negative)'),...
    'Longitude', struct('dims', {{'cities'}}, 'units', 'degrees (west is negative)'),...
    'WindDir', struct('dims', {{'cities', 'days', 'orbits'}}, 'units', 'degrees CCW from east', 'description', 'This gives the direction the wind is blowing towards; i.e. a 45 deg value means the wind blows southwest to northeast'),...
    'WindSpeed', struct('dims', {{'cities', 'days', 'orbits'}}, 'units', 'm s^-1'),...
    'WindUsedBool', struct('dims', {{'cities', 'days', 'orbits'}}, 'description', 'Set to 1 if the corresponding winds meet the criteria for this day to be included in the line densities', 'type', 'int8'),...
    'x', struct('dims', {{'cities', 'line_density'}}, 'name', 'linedens_x', 'units', 'km', 'description', 'Distance from the geographic center of the city'),...
    'linedens', struct('dims', {{'cities', 'line_density'}}, 'name', 'no2_line_density', 'units', 'mol km^-1'),...
    'ffit', struct('dims', {{'cities', 'fitting_parameters'}}, 'name', 'emg_fit_parameters', 'description', 'The fitting parameters for the EMG fits (a, x_0, mu_x, sigma_x, and B, in order)', 'conv', @ffit2vec),...
    'emgfit', struct('dims', {{'cities', 'line_density'}}, 'units', 'mol km^-1', 'description', 'The EMG fit on the same x-coordinates as the line density'),...
    'sd', struct('dims', {{'cities', 'fitting_parameters'}}, 'name', 'fit_param_std_dev', 'description', 'The standard deviations for the five fitting parameters'),...
    'percentsd', struct('dims', {{'cities', 'fitting_parameters'}}, 'name', 'fit_param_percent_std_dev', 'description', 'The percent standard deviations for the five fitting parameters'),...
    'ci95', struct('dims', {{'cities', 'fitting_parameters'}}, 'name', 'fit_param_conf_int', 'description', 'The 95% confidence intervals of the five fitting parameters'),...
    'percent_ci95', struct('dims', {{'cities', 'fitting_parameters'}}, 'name', 'fit_param_percent_conf_int', 'description', 'The 95% confidence intervals of the five fitting parameters as a percent of the value'),...
    'r2', struct('dims', {{'cities'}}, 'description', 'R2 of the EMG fit'),...
    'tau', struct('dims', {{'cities'}}, 'units', 'h', 'description', 'Lifetime derived from the EMG fit'),...
    'tau_uncert', struct('dims', {{'cities'}}, 'units', 'h', 'description', 'Uncertainty in the lifetime as a 95% confidence interval'),...
    'tau_sd', struct('dims', {{'cities'}}, 'units', 'h', 'description', 'Uncertainty in the lifetime as the standard deviation'),...
    'n_dofs', struct('dims', {{'cities'}}, 'description', 'Number of degrees of freedom in the lifetime fits'));
    
var_fields = fieldnames(variables);
for i_var = 1:numel(var_fields)
    vfn = var_fields{i_var};
    vinfo = variables.(vfn);
    
    nc_name = vfn;
    if isfield(vinfo, 'name')
        nc_name = vinfo.name;
    end
    
    nc_type = 'double';
    if isfield(vinfo, 'type')
        nc_type = vinfo.type;
    end
    
    [dim_cell, req_size] = make_dim_cell(vinfo.dims);
    var_array = init_array(nc_type, req_size);
    
    for i_loc = 1:numel(locs)
        % Get the variable and pad, if necessary
        this_var_data = find_substruct_field(locs(i_loc), vfn);
        % Imaginary values (sometimes happens in uncertainty) should be
        % fills
        if isnumeric(this_var_data)
            xx_imag = imag(this_var_data) ~= 0;
            this_var_data(xx_imag) = nan;
        end
        if isfield(vinfo, 'conv')
            this_var_data = vinfo.conv(this_var_data);
        end
        this_var_data = pad2size(this_var_data, req_size, nc_type);
        if isempty(this_var_data)
            continue
        end
        var_array(i_loc, :) = this_var_data(:);
    end
    
    
    nccreate(nc_file_name, nc_name, 'Dimensions', dim_cell, 'Datatype', nc_type);
    ncwrite(nc_file_name, nc_name, var_array);
    for i_att = 1:numel(atts_to_write)
        att = atts_to_write{i_att};
        if isfield(vinfo, att)
            ncwriteatt(nc_file_name, nc_name, att, vinfo.(att));
        end
    end
end

    function [dimcell, req_dims] = make_dim_cell(dimnames)
        dimcell = cell(1, numel(dimnames)*2);
        % every variable has the cities dimension, we need to check the
        % others
        req_dims = ones(1, numel(dimnames)-1);
        for i=1:numel(dimnames)
            j = (i-1)*2 + 1;
            dimcell{j} = dimnames{i};
            dimcell{j+1} = dims.(dimnames{i});
            if i > 1
                req_dims(i-1) = dimcell{j+1};
            end
        end
    end

    function arr = init_array(dtype, req_size)
        if isempty(req_size)
            req_size = 1;
        end
        switch dtype
            case 'char'
                arr = repmat(' ', [numel(locs), req_size]);
            otherwise
                arr = nan([numel(locs), req_size]);
        end
end

end



function arr = pad2size(arr, req_size, dtype)
switch dtype
    case 'char'
        arr = pad_char(arr, req_size);
    otherwise
        arr = pad_numeric(arr, req_size);
end
end

function arr = pad_numeric(arr, req_size)
if isrow(arr)
    arr = arr';
end
for i_dim = 1:numel(req_size)
    sz = size(arr);
    sz(i_dim) = req_size(i_dim) - sz(i_dim);
    padding = nan(sz);
    arr = cat(i_dim, arr, padding);
end
end

function str = pad_char(str, req_size)
padding = repmat(' ', 1, req_size - numel(str));
str = [str, padding];
end

function v = ffit2vec(ffit)
v = struct2array(ffit);
end