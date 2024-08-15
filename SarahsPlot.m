function SarahsPlot(MS,xarg,yarg,PlotColor,PlotTitle)
    figure; hold on;
    plot([MS.(xarg)],[MS.(yarg)],'.','Color',PlotColor);
    title(PlotTitle);
end