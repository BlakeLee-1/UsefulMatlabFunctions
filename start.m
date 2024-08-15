set(0,'defaultAxesFontSize',20);
set(0,'defaultaxeslinewidth',1.2);
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on');
set(groot,'defaultLineLineWidth',2.5);
set(0, 'DefaultAxesBox', 'on');
set(0,'defaultAxesFontSize',20);
set(groot,'defaultFigurePosition', [400,250,900,750])
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))

set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))