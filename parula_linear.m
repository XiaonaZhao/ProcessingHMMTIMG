function map = parula_linear(m)
%PARULA_LINEAR   Turbo colormap.
%   TURBO(M) returns an M-by-3 matrix containing the turbo colormap, a
%   variant of the jet colormap that is more perceptually uniform.
%
%   See also JET, COLORMAP.

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

mycolorpoint = [
    0.6350, 0.0780, 0.1840;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4660, 0.6740, 0.1880;
    0.3010, 0.7450, 0.9330;
    0, 0.4470, 0.7410;
    ];

% mycolorposition=[1 9 24 40 56 64];
% mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:64,'linear','extrap');
% mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:64,'linear','extrap');
% mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:64,'linear','extrap');
% mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/256;
% mycolor = round(mycolor*10^4)/10^4; 

P = size(mycolorpoint,1);
map = interp1(1:P, mycolorpoint, linspace(1,P,m), 'linear');
    
end