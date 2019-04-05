function [ fighandle, h ] = plot_error_envelope_y( x_in, y_lower, y_upper, varargin )
%plot_error_envelope_y(x_in, y_lower, y_upper, [opt] colorspec, [opt] figure handle): Uses the fill function to plot an error envelope for error in x value.
%   This plots an error envelope on your plot for error in the x_values.
%   This is set up to allow the use of standard deviation/std. error or
%   quartiles, hence the lower and upper bounds are input separately.  The
%   default color of the envelope is light grey, but the color can be
%   specified as an optional argument.  Most likely this will need to be
%   called before you plot your actual data, as the fill will probably
%   cover the line if plotted afterwards.
%
%   Returns the figure handle and handle to the fill object.

p = advInputParser;

p.addOptional('colorspec',[0.8 0.8 0.8], @(t) (isnumeric(t) || ischar(t)));
p.addParameter('alpha', 0.5);
p.addParameter('parent', gca);
p.parse(varargin{:});
pout = p.Results;

alpha = pout.alpha;
ax = pout.parent;

%Convert all to column vectors
x_in = x_in(:);
y_lower = y_lower(:);
y_upper = y_upper(:);

colspec = pout.colorspec;

Y = [y_lower; flipud(y_upper)];
X = [x_in; flipud(x_in)];

h=fill(ax,X,Y,colspec,'edgecolor','none','facealpha',alpha);
hold off

end

