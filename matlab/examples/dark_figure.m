function h = dark_figure(n)
    if(nargin==0)
        h=figure();
    else
        h=figure(n);
    end
    
    set(gcf,'color',[0 0 0]);
    set(gcf,'defaultaxescolor',[0 0 0]);
    set(gcf,'defaultaxesXcolor',[1 1 1]);
    set(gcf,'defaultaxesYcolor',[1 1 1]);
    set(gcf,'defaultaxesZcolor',[1 1 1]);    
    set(gcf,'defaulttextcolor',[1 1 1]);
    set(gcf,'inverthardcopy','off');
    set(gcf,'defaultaxesGridColor',[0.85 0.85 0.85]);
end
