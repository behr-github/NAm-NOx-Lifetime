function [outputArg1,outputArg2] = write_avgs_to_ncdf(ncdir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

years = 2005:2014;

for iyr = 1:numel(years)
    yr = years(iyr);
    avg_file = misc_emissions_analysis.avg_file_name(yr, 'TWRF');
    avgs = load(avg_file);
    wkend_avg_file = misc_emissions_analysis.avg_file_name(yr, 'US');
    wkend_avgs = load(wkend_avg_file);
    [~, nc_filename] = fileparts(avg_file);
    nc_filename = fullfile(ncdir, [nc_filename, '.nc']);
    delete(nc_filename);
    write_one_file(avgs, wkend_avgs, nc_filename);
end

end

function write_one_file(wkday_avgs, wkend_avgs, nc_file)
wkday_avgs = wkday_avgs.daily;
wkend_avgs = wkend_avgs.daily;
lon = wkday_avgs.lon(1,:)';
lon_dim = {'lon', length(lon)};
lat = wkday_avgs.lat(:,1);
lat_dim = {'lat', length(lat)};
full_dims = [lat_dim, lon_dim];

nccreate(nc_file, 'lon', 'Dimensions', lon_dim, 'Datatype', 'double');
ncwrite(nc_file, 'lon', lon);
ncwriteatt(nc_file, 'lon', 'units', 'degrees_east');
ncwriteatt(nc_file, 'lon', 'description', 'Longitude coordinate for average grids');

nccreate(nc_file, 'lat', 'Dimensions', lat_dim, 'Datatype', 'double');
ncwrite(nc_file, 'lat', lat);
ncwriteatt(nc_file, 'lat', 'units', 'degrees_north');
ncwriteatt(nc_file, 'lat', 'description', 'Latitude coordinate for average grids');

nccreate(nc_file, 'weekday_no2', 'Dimensions', full_dims, 'Datatype', 'double');
ncwrite(nc_file, 'weekday_no2', wkday_avgs.no2);
ncwriteatt(nc_file, 'weekday_no2', 'units', 'molec. cm^{-3}');
ncwriteatt(nc_file, 'weekday_no2', 'description', 'BEHR NO2 VCDs averaged between 1 Apr and 30 Sept, weekdays (TWRF)');

nccreate(nc_file, 'weekday_weights', 'Dimensions', full_dims, 'Datatype', 'double');
ncwrite(nc_file, 'weekday_weights', wkday_avgs.weights);
ncwriteatt(nc_file, 'weekday_weights', 'units', 'unitless');
ncwriteatt(nc_file, 'weekday_weights', 'description', 'VCD weights for weekdays (TWRF)');

nccreate(nc_file, 'weekend_no2', 'Dimensions', full_dims, 'Datatype', 'double');
ncwrite(nc_file, 'weekend_no2', wkend_avgs.no2);
ncwriteatt(nc_file, 'weekend_no2', 'units', 'molec. cm^{-3}');
ncwriteatt(nc_file, 'weekend_no2', 'description', 'BEHR NO2 VCDs averaged between 1 Apr and 30 Sept, weekends (US)');

nccreate(nc_file, 'weekend_weights', 'Dimensions', full_dims, 'Datatype', 'double');
ncwrite(nc_file, 'weekend_weights', wkend_avgs.weights);
ncwriteatt(nc_file, 'weekend_weights', 'units', 'unitless');
ncwriteatt(nc_file, 'weekend_weights', 'description', 'VCD weights for weekends (US)');
end