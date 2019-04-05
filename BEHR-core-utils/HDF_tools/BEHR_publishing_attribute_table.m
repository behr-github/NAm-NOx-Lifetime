function [ attr_table ] = BEHR_publishing_attribute_table( varargin )
%BEHR_PUBLISHING_ATTRIBUTE_TABLE Returns table of HDF attributes for BEHR fields
%   Detailed explanation goes here

E = JLLErrors;

as_struct = ismember('struct', varargin);

allowed_subsets = {'all', 'sp', 'behr', 'behr-insitu', 'pub', 'pub-insitu'};
xx = ismember(varargin, allowed_subsets);
if sum(xx) > 1
    E.badinput('Only one of %s may be an input', strjoin(allowed_subsets, ', '));
elseif sum(xx) < 1
    subset = 'all';
else
    subset = varargin{xx};
end 


% This cell array will have the variable name, unit, range, fill, product
% (SP or BEHR) and description in that order.
longfill = single(-1.267650600228229401496703205376e30);
shortfill = single(-32767);
pixcorfill = single(-1.00000001504747e+30);
behrfill = single(behr_fill_val());
behrflagfill = bitset(0, 32);
nofill = NaN;

% Variable name, unit, range, fill value, product, description

attr_table = {  'AmfStrat', 'unitless', [0, Inf], longfill, 'SP', 'Stratospheric AMF';...
                'AmfTrop', 'unitless', [0, Inf], longfill, 'SP', 'Tropospheric AMF (standard product)';...
                'Areaweight', 'km^-2', [0, Inf], behrfill, 'BEHR', 'Reciprocal of pixel area; use to weight a temporal average of grids to account for pixel representativeness';...
                'BEHRAMFTrop', 'unitless', [0, Inf], behrfill, 'BEHR', 'Tropospheric AMF (BEHR) for total NO2 (including ghost column), calculated with MODIS Albedo, GLOBE+WRF Terr. Pres., and 12 km NO2 profiles';...
                'BEHRAMFTropVisOnly', 'unitless', [0, Inf], behrfill, 'BEHR', 'Tropospheric AMF (BEHR) for visible NO2 only, calculated with MODIS Albedo, GLOBE+WRF Terr. Pres., and 12 km NO2 profiles';...
                'BEHRColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], behrfill, 'BEHR', 'Tropospheric NO2 VCD (BEHR) including estimated ghost column, calculated as SCD_trop / AMF_BEHR';...
                'BEHRColumnAmountNO2TropVisOnly', 'molec./cm^2', [0, Inf], behrfill, 'BEHR', 'Tropospheric NO2 VCD (BEHR) with only visible NO2 included, calculated as SCD_trop / AMF_BEHR_VIS';...
                'BEHRScatteringWeightsClear', 'unitless', [0, Inf], behrfill, 'BEHR', 'Clear-sky scattering weights derived from the MODIS albedo and GLOBE surface pressure. Includes NO2 cross section temperature correction.';...
                'BEHRScatteringWeightsCloudy', 'unitless', [0, Inf], behrfill, 'BEHR', 'Cloudy-sky scattering weights derived from the MODIS albedo and GLOBE surface pressure. Includes NO2 cross section temperature correction.';...
                'BEHRSurfacePressure', 'hPa', [0, 1013], behrfill, 'BEHR', 'Surface pressure used in the BEHR AMF calculation, is the WRF surface pressure adjusted with the hypsometric equation and the difference in WRF and GLOBE terrain heights.';...
                'WRFSurfacePressure', 'hPa', [0, Inf], behrfill, 'BEHR', 'Surface pressure from the WRF model.';...
                'BEHRAvgKernels', 'unitless', [0, Inf], behrfill, 'BEHR', 'Averaging kernels computed for the weighted average of cloudy and clear conditions';...
                'BEHRPressureLevels', 'hPa', [0, Inf], behrfill, 'BEHR', 'Pressure levels that correspond to the scattering weight, averaging kernel, and NO2 a priori vectors';...
                'BEHRQualityFlags', 'bit array flag', [0 (2^32)-1], behrflagfill, 'BEHR', 'Quality flags field summarizing NASA and BEHR flags. Do not use pixel if least significant bit != 0 (i.e. if value is odd)';...
                'CloudFraction', 'unitless', [0, 1], shortfill, 'SP', 'OMI geometric cloud fraction';...
                'CloudPressure', 'hPa', [0, Inf], shortfill, 'SP', 'OMI cloud top pressure';...
                'CloudRadianceFraction', 'unitless', [0, 1], shortfill, 'SP', 'OMI cloud radiance (top of atmosphere light fraction)';...
                'ColumnAmountNO2', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Total NO2 VCD';...
                'ColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Tropospheric NO2 VCD (standard product)';...
                'ColumnAmountNO2TropStd', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Standard deviation of SP NO2 tropospheric VCD';...
                'ColumnAmountNO2Strat', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Stratospheric NO2 VCD';...
                'GLOBETerrainHeight', 'm', [-Inf, Inf], behrfill, 'BEHR', 'Terrain height derived from GLOBE (1 x 1 km) topography averaged to OMI pixel';...
                'FoV75Area', 'km^2', [0, Inf], pixcorfill,'PIXCOR', 'Area of pixel based on the part of the Earth from which 75% of the observed radiance is received';...
                'FoV75CornerLatitude', 'deg', [-90, 90], pixcorfill, 'PIXCOR', 'Latitude corners of pixel based on the part of the Earth from which 75% of the observed radiance is received';...
                'FoV75CornerLongitude', 'deg', [-180, 180], pixcorfill, 'PIXCOR', 'Longitude corners of pixel based on the part of the Earth from which 75% of the observed radiance is received';...
                'Latitude', 'deg', [-90, 90], longfill, 'SP', 'Center latitude of pixels';...
                'Longitude', 'deg', [-180, 180], longfill, 'SP', 'Center longitude of pixels';...
                'MODISAlbedo', 'unitless', [0, 1], behrfill, 'BEHR', 'Surface reflectance derived from MCD43D07-09 products, avg. to OMI pixel';...
                'MODISCloud', 'unitless', [0, 1], behrfill, 'BEHR', 'MODIS MYD06_L2 5 x 5 km cloud fraction, avg. to OMI pixel';...
                'RelativeAzimuthAngle', 'deg', [0, 180], behrfill, 'BEHR', 'Calculated azimuth angle between sun and satellite';...
                'Row', 'unitless', [0, 59], nofill, 'SP', 'Across track row number, 0 based';...
                'SlantColumnAmountNO2', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Total NO2 SCD';...
                'SolarAzimuthAngle', 'deg', [-180, 180], longfill, 'SP', 'Solar azimuth angle';...
                'SolarZenithAngle', 'deg', [0 90], longfill, 'SP', 'Solar zenith angle';...
                'SpacecraftAltitude', 'm', [0, Inf], longfill, 'SP', 'Spacecraft (Aura satellite) altitude above WGS84 ellipsoid';...
                'SpacecraftLatitude', 'deg', [-90, 90], longfill, 'SP', 'Spacecraft (Aura satellite) latitude above WGS84 ellipsoid';...
                'SpacecraftLongitude', 'deg', [-180, 180], longfill, 'SP', 'Spacecraft (Aura satellite) longitude above WGS84 ellipsoid';...
                'Swath', 'unitless', [0, Inf], nofill, 'SP', 'Orbit number since OMI launch';...
                'TerrainHeight', 'm', [-Inf, Inf], shortfill, 'SP', 'Terrain height';...
                'TerrainPressure', 'hPa', [0, Inf], shortfill 'SP', 'Terrain pressure';...
                'TerrainReflectivity', 'unitless', [0, 1], shortfill, 'SP', 'Terrain albedo (OMI albedo product)';...
                'TiledArea', 'km^2', [0, Inf], pixcorfill, 'PIXCOR', 'Pixel area assuming that pixels are perfectly tiled (do not overlap)';...
                'TiledCornerLatitude', 'km^2', [-90, 90], pixcorfill, 'PIXCOR', 'Latitude corners of pixel assuming that pixels are perfectly tiled (do not overlap)';...
                'TiledCornerLongitude', 'km^2', [-180, 180], pixcorfill, 'PIXCOR', 'Longitude corners of pixel assuming that pixels are perfectly tiled (do not overlap)';...
                'Time', 's', [0, Inf], longfill, 'SP', 'Time at start of scan (TAI93: seconds since Jan 1, 1993)';...
                'ViewingAzimuthAngle', 'deg', [-180, 180], longfill, 'SP', 'Viewing azimuth angle';...
                'ViewingZenithAngle', 'deg', [0, 90], longfill, 'SP', 'Viewing zenith angle';...
                'XTrackQualityFlags', 'bit array flag', 'N/A', uint8(255), 'SP', 'Across track quality flag (for row anomaly)';...
                'VcdQualityFlags', 'bit array flag', 'N/A', uint16(65535), 'SP', 'Ground pixel quality flags';...
                'InSituAMF', 'unitless', [0, Inf], behrfill, 'BEHR-InSitu', 'AMF calculated using co-located in situ NO2 profile';...
                'BEHR_R_ColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-InSitu', 'BEHR Tropospheric NO2 VCD calculated with the in situ AMF';...
                'ProfileCount', 'unitless', [0, Inf], behrfill, 'BEHR-InSitu', 'Number of aircraft profiles averaged to create the in situ a priori NO2 profile';...
                'InSituFlags', 'bit array flag', 'N/A', nofill, 'BEHR-InSitu', 'In situ profile quality flag';...
                'BEHRColumnAmountNO2Trop_L3', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-L3','BEHR tropospheric NO2 VCDs filtered for quality and row anomaly';...
                'BEHRColumnAmountNO2Trop_L3MODISCloud', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-L3','BEHR tropospheric NO2 VCDs additionally filtered for MODIS Cloud < 20%';...
                'BEHRNO2apriori', 'parts-per-part', [-Inf, Inf], behrfill, 'BEHR', 'NO2 a priori profile used for each pixel. Pressure levels given in BEHRPressureLevels.';...
                'BEHRTropopausePressure','hPa',[0, Inf], behrfill, 'BEHR', 'Tropopause pressure derived from WRF temperature profiles; used in the BEHR AMF calculation';...
                };
            
attr_table = add_psm_weight_fields(attr_table);
attr_table = choose_subset(attr_table, subset);
            
if numel(unique(attr_table(:,1))) < size(attr_table, 1)
    E.callError('attr_def', 'One or more attributes is multiply defined in the attributes table');
end
            
if as_struct
    fields = {'unit', 'range', 'fillvalue', 'product', 'description'};
    S = struct;
    for a = 1:size(attr_table,1)
        S.(attr_table{a,1}) = cell2struct(attr_table(a,2:end), fields, 2);
    end
    attr_table = S;
end

%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%
    function attr_table = add_psm_weight_fields(attr_table)
        weights_vars = BEHR_publishing_gridded_fields.psm_weight_vars;
        table_lines = cell(numel(weights_vars), size(attr_table, 2));
        for i=1:numel(weights_vars)
            if ismember(weights_vars{i}, BEHR_publishing_gridded_fields.psm_weight_vars('insitu'))
                product = 'BEHR-InSitu';
            else
                product = 'BEHR';
            end

            table_lines(i,:) = [weights_vars(i), {'unitless', [0, Inf], behrfill, product, sprintf('Weight field for the %s field', BEHR_publishing_gridded_fields.all_psm_vars{i})}];
        end
        attr_table = cat(1, attr_table, table_lines);
    end

    
    function attr_table = choose_subset(attr_table, subset)
        if strcmpi(subset, 'all')
            return
        elseif strcmpi(subset, 'pub')
            products = {'SP', 'BEHR', 'PIXCOR'};
        elseif strcmpi(subset, 'pub-insitu')
            products = {'SP', 'BEHR', 'PIXCOR', 'BEHR-InSitu'};
        else
            products = {subset};
        end

        xx = false(size(attr_table,1));

        for a=1:size(attr_table,1)
            xx(a) = any(strcmpi(attr_table{a,5}, products));
        end

        attr_table = attr_table(xx,:);
    end

end

