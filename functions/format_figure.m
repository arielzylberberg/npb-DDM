function format_figure(hfig,varargin)
%format_figure(hfig,varargin)

hax = get(hfig,'children');

invert_colors = 0;
% FontSize = 6.8;
FontSize = 14;
% LineWidthAxes = 0.335;
LineWidthAxes = 0.5;
LineWidthPlot = 1;
BackColor = [1 1 1];
for i=1:length(varargin)
    if strcmp(varargin{i},'presentation')
        FontSize = 18;
        LineWidthAxes = 1;
%         hax = get(gcf,'Children');
        LineWidthPlot = 2.5;
%         for j = 1:length(hax)
%             set(hax(j),'LineWidth',2.5)
%         end
        
    elseif strcmp(varargin{i},'invert_colors')
        invert_colors = 1;
    elseif strcmp(varargin{i},'FontSize')
        FontSize = varargin{i+1};
    elseif strcmp(varargin{i},'BackgroundColor')
        BackColor = varargin{i+1};
    elseif strcmp(varargin{i},'LineWidthPlot')
        LineWidthPlot = varargin{i+1};    
    elseif strcmp(varargin{i},'LineWidthAxes')
        LineWidthAxes = varargin{i+1};    
    end
end

all_text = findall(hfig,'Type','text');
set(all_text,'FontSize',FontSize)

for i=1:length(hax)
    set(hax(i),'FontSize',FontSize,'LineWidth',LineWidthAxes,'box','off')
    try
        set(get(hax(i),'Children'),'LineWidth',LineWidthPlot)
    catch
        
    end
end

set(gcf,'color',BackColor);

if invert_colors==1
    a = findall(hfig);
    w = findobj(a,'Color','w');
    b = findobj(a,'Color','k');
    set(w,'Color','k');
    set(b,'Color','w');
    
    for j=1:length(hax)
        set(hax(j),'Ycolor','w')
        set(hax(j),'Xcolor','w')
        set(hax(j),'Color','none') % Para fondo
        %transparente
    end
    
    set(hfig,'Color','none') % Para fondo transparente
    set(hfig,'InvertHardcopy','off')
end