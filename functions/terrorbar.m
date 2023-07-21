function h = terrorbar(x,y,err,varargin)
% h = terrorbar(x,y,err,varargin)
% error bar without the annoying horizontal line

h = errorbar(x,y,err,varargin{:});
set(h,'CapSize',0);

end
