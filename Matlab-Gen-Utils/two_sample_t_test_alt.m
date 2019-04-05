function [t_calc, t_table, is_different] = two_sample_t_test_alt(val1, sigma1, n1, val2, sigma2, n2, varargin)
%TWO_SAMPLE_T_TEST_ALT Alternate test if values with different std. dev. differ 
%   [T_CALC, T_TABLE, IS_DIFFERENT] = TWO_SAMPLE_T_TEST_ALT( VAL1, SIGMA1,
%   N1, VAL2, SIGMA2, N2 ) computes the calculated t-value (T_CALC) for the
%   values (VAL1 and VAL2), their sample standard deviations (SIGMA1 and
%   SIGMA2), and the number of points that went into each measurement (N1
%   and N2). It also returns the tablulated t-value (T_TABLE) for the
%   effective degrees of freedom for these values and a two-sided 95%
%   confidence level, as well as a boolean (IS_DIFFERENT) that is true if
%   the two values are different at the 95% confidence levels.
%
%   Parameters:
%       'confidence' - the confidence level, expressed as a decimal. 0.95
%       by default, i.e. 95% confidence level.
%
%       'sided' - whether to use a one sided or two sided t-test. Default
%       is 'two', must be either 'one' or 'two'.
%
%   From Harris "Quantitative Chemical Analysis" Ch 4. (p. 77-78), when you
%   want to test if two values are statistically different, but the two
%   samples have different standard deviations, a good estimate of the
%   calculated t value and degrees of freedom is:
%
%       t_calc =          | v1 - v2 |
%               ----------------------------
%                sqrt(s1^2 / n1 + s2^2 / n2)
%
%       DoF    = (s1^2 / n1 + s2^2 / n2)^2
%               ---------------------------
%                (s1^2/n1)^2   (s2^2/n2)^2
%                ----------- + -----------
%                  n1 - 1        n2 - 1
%
%   where v1 and v2 are the values, s1 and s2 are the sample standard
%   deviations, and n1 and n2 are the number of measurements that go into
%   each value.

%%%%%%%%%%%%%%%%%
% INPUT PARSING %
%%%%%%%%%%%%%%%%%

E = JLLErrors;

p = inputParser;
p.addParameter('confidence', 0.95);
p.addParameter('sided', 'two');

p.parse(varargin{:});
pout = p.Results;

confidence_level = pout.confidence;
t_sidedness = pout.sided;

if ~isnumeric(val1) || ~isscalar(val1)
    E.badinput('VAL1 must be a numeric scalar')
end
if ~isnumeric(sigma1) || ~isscalar(sigma1)
    E.badinput('SIGMA1 must be a numeric scalar')
end

if ~isnumeric(n1) || ~isscalar(n1) || n1 < 1 || isnan(n1)
    E.badinput('N1 must be a positive, numeric scalar')
end


if ~isnumeric(val2) || ~isscalar(val2)
    E.badinput('VAL2 must be a numeric scalar')
end
if ~isnumeric(sigma2) || ~isscalar(sigma2)
    E.badinput('SIGMA2 must be a numeric scalar')
end
if ~isnumeric(n2) || ~isscalar(n2) || n2 < 1 || isnan(n2)
    E.badinput('N2 must be a positive, numeric scalar')
end

if ~isnumeric(confidence_level) || ~isscalar(confidence_level) || confidence_level < 0 || confidence_level > 1
   E.badinput('"confidence_level" must be a scalar number between 0 and 1');
end

allowed_sidedness = {'one','two'};
if ~ischar(t_sidedness) || ~ismember(t_sidedness, allowed_sidedness)
    E.badinput('"sided" must be one of: ''%s''', strjoin(allowed_sidedness, ''', '''));
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

denom = sqrt(sigma1.^2 ./ n1 + sigma2.^2 ./ n2);
t_calc = abs(val1 - val2) ./ denom;

dof_num = (sigma1.^2 ./ n1 + sigma2.^2 ./ n2).^2;
dof_denom = (sigma1.^2 ./ n1).^2 ./ (n1 - 1) + (sigma2.^2 ./ n2).^2 ./ (n2 - 1);
n_dofs = dof_num ./ dof_denom;

if strcmpi(t_sidedness, 'two')
    % tinv gives 1-sided t values. For a two sided t-test (i.e. one where
    % the values may be above or below each other) the same percent chance
    % of being different must be split onto either side. Concretely, for a
    % 95% two-sided confidence interval, we need to compare against the
    % 97.5% t-value returned by tinv.
    confidence_level = 1 - ( (1 - confidence_level)/2 );
end

t_table = tinv(confidence_level, n_dofs);

is_different = t_calc > t_table;

end

