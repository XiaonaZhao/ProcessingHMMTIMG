function mycolor = SublimeVivid(m)
%JSHINE   JShine colormap.
%   JSHINE(M) returns an M-by-3 matrix containing the JShine colormap, a
%   variant of the jet colormap that is more perceptually uniform.
%   https://uigradients.com/#JShine
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
    66,94,249;
    95,90,226;
    105,88,219;
    137,85,180;
    157,82,180;
    176,80,165;
    197,77,148;
    217,74,133;
    249,70,109;
    ];
mycolorposition=[1890 1589 1491 1170  964 774 551 353 29];
mm = 1920*ones(1, length(mycolorposition));
mycolorposition = ceil((mm- mycolorposition)/1920*64);

mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:64,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:64,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:64,'linear','extrap');
mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/256;
mycolor = round(mycolor*10^4)/10^4; %保留4位小数

P = size(mycolor,1);
mycolor = interp1(1:P, mycolor, linspace(1,P,m), 'linear');
end