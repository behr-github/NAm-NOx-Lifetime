function [this_rate] = interp_tuv_photolysis(j_handle, wrf_datetime, wrf_lon, wrf_lat, varargin)
%INTERP_TUV_PHOTOLYSIS Calculate TUV photolysis rates and interpolate to a grid
%   Calling TUV is slow, and can only be done for a single lat/lon at a
%   time, so if we need to call it for every point on a WRF grid, that
%   would take a long time. However, the photolysis rates shouldn't change
%   too much with lat/lon, at least around noon, or at least should change
%   in a way that we can reasonably interpolate between a few points. This
%   function will therefore call TUV for a limited number of points on the
%   input grids and interpolate the photolysis rates to the entire grid.
%
%   RATE = INTERP_TUV_PHOTOLYSIS( J_HANDLE, WRF_DATETIME, WRF_LON, WRF_LAT )
%       - J_HANDLE should be a photlysis rate handle from a
%       KPP_OOP.Reaction.
%       - WRF_DATETIME should be a date number or a date string implicitly
%       understood by Matlab that gives the date and UTC hour that the
%       photolysis rate should be calculated for.
%       - WRF_LON and WRF_LAT should be 2D arrays of longitudes and
%       latitudes that you want the J-values interpolated to.
%
%   RATE = INTERP_TUV_PHOTOLYSIS( ___, WRF_ALT ) 
%       WRF_ALT is an optional argument that gives the altitude that the
%       photolysis rate should be calculated for, in kilometers. Default is
%       0.5. Currently, this must be a scalar value. 
%
%   Additional parameters:
%
%       'nx', 'ny' - the number of points in the longitudinal and
%       latitudinal directions, respectively, that TUV rates should be
%       explicitly calculated for. These will be evenly spaced along the
%       WRF grid. Default is 3 for each. Increasing this will improve
%       accuracy, but at the cost of increased compute time.

p = advInputParser;
p.addOptional('wrf_alt', 0.5);
p.addParameter('nx', 3);
p.addParameter('ny', 3);

p.parse(varargin{:});
pout = p.Results;

wrf_alt = pout.wrf_alt;
nx = pout.nx;
ny = pout.ny;

[wrf_sparse_lon, wrf_sparse_lat] = make_sparse_latlon(wrf_lon, wrf_lat, nx, ny);
wrf_date = datestr(wrf_datetime, 'yyyy-mm-dd');
wrf_hour = hour(wrf_datetime);
j_rate = j_handle(wrf_date, wrf_hour, wrf_sparse_lon, wrf_sparse_lat, true, wrf_alt);

% Calling TUV is slow, but a quick test
% interpolating from a 3x2 matrix to a 7x3 one for
% 1 June 2005 20:00:00 gave me < 5% errors, so
% we'll use interpolation.
SInterp = scatteredInterpolant(wrf_sparse_lon, wrf_sparse_lat, j_rate);
this_rate = SInterp(wrf_lon, wrf_lat);

end

function [sparse_lon, sparse_lat] = make_sparse_latlon(lon, lat, nx, ny)
xx = round(linspace(1, size(lon,1), nx));
yy = round(linspace(1, size(lon,2), ny));
sparse_lon = lon(xx, yy);
sparse_lat = lat(xx, yy);
end