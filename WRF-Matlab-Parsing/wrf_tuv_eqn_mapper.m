function [ tuv_eqn ] = wrf_tuv_eqn_mapper( wrf_j )
%[ TUV_EQN ] WRF_TUV_EQN_MAPPER( WRF_J )
%   Returns the equation from TUV that corresponds to the desired WRF-Chem
%   photolysis reaction. WRF_J must be the name of the WRF-Chem KPP
%   photolysis constant, that is, if in the .eqn file, it is listed as
%   'j(Pj_no2)' pass in j(Pj_no2) or Pj_no2
E = JLLErrors;
[s,e] = regexp(wrf_j, '\(.*\)');
if ~isempty(s)
    wrf_j = wrf_j(s+1:e-1);
end

switch wrf_j
    case 'Pj_no2'
        tuv_eqn = 'NO2 -> NO + O(3P)';
    case 'Pj_o33p'
        tuv_eqn = 'O3 -> O2 + O(3P)';
    case 'Pj_o31d'
        tuv_eqn = 'O3 -> O2 + O(1D)';
    case 'Pj_h2o2'
        tuv_eqn = 'H2O2 -> 2 OH';
    case 'Pj_no3o2'
        tuv_eqn = 'NO3 -> NO + O2';
    case 'Pj_no3o'
        tuv_eqn = 'NO3 -> NO2 + O(3P)';
    case 'Pj_hno2'
        tuv_eqn = 'HNO2 -> OH + NO';
    case 'Pj_hno3'
        tuv_eqn = 'HNO3 -> OH + NO2';
    case 'Pj_hno4'
        tuv_eqn = 'HNO4 -> HO2 + NO2';
    case 'Pj_ch2om'
        tuv_eqn = 'CH2O -> H2 + CO';
    case 'Pj_ch2or'
        tuv_eqn = 'CH2O -> H + HCO';
        issue_warning(wrf_j, tuv_eqn);
    case 'Pj_ch3coch3'
        tuv_eqn = 'CH3COCH3 -> CH3CO + CH3';
    otherwise
        E.notimplemented('No mapping defined for "%s"', wrf_j);
end

end

function issue_warning(wrf_j, tuv_eqn)
warning('I am uncertain if the the mapping "%s -> %s" is correct', wrf_j, tuv_eqn);
end