function subplot_stretch(m,n,varargin)
%SUBPLOT_STRETCH Stretch a figure to better show subplots
%
%   SUBPLOT_STRETCH(M,N) For the current figure, with M-by-N subplots,
%   expand it so each subplot is roughly a normal figure size.
%
%   Parameters:
%       'figh' - a handle to the figure to stretch. Default is the current
%       figure.
%
%       'factor' - an extra factor to apply to both dimensions. Default is
%       1.
p = inputParser;
p.addParameter('figh', gcf);
p.addParameter('factor',1);

p.parse(varargin{:});
pout = p.Results;

fig = pout.figh;
factor = pout.factor;

fig.Position(3) = fig.Position(3)*factor*n;
fig.Position(4) = fig.Position(4)*factor*m;

end

