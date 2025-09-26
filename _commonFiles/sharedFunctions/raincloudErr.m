function h = raincloudErr(x,col,xpos,boxwidth,markersize,jit,fig)

% set fig as current fig
set(0, 'CurrentFigure', fig)
box on;

% jitter for raindrops
a = boxwidth;
b = -a;
jitter = ((b-a).*rand(length(x),1) + a)*jit;

% info for making boxplot
Y = quantile(x,[0.25 0.75 0.5 0 1]);

% 'box' of 'boxplot'
hold on;
h{1} = rectangle('Position',[xpos-(boxwidth*0.5) Y(1) boxwidth Y(2)-Y(1)]);
set(h{1},'EdgeColor','k');
set(h{1},'FaceColor',col);
set(h{1},'LineWidth',1);

% median line
ml_size = 0.35; %0.35
h{2} = line([xpos-(boxwidth*ml_size) xpos+(boxwidth*ml_size)],[Y(3) Y(3)],'col','k','LineWidth',3);

% whiskers
h{3} = line([xpos xpos],[Y(2) Y(5)],'col','k','LineWidth',1);
h{4} = line([xpos-(boxwidth*0.15) xpos+(boxwidth*0.15)],[Y(5) Y(5)],'col','k','LineWidth',1);
h{5} = line([xpos xpos],[Y(1) Y(4)],'col','k','LineWidth',1);
h{6} = line([xpos-(boxwidth*0.15) xpos+(boxwidth*0.15)],[Y(4) Y(4)],'col','k','LineWidth',1);

% raindrops
h{7} = scatter(jitter+xpos,x);
h{7}.SizeData = markersize;
h{7}.Marker = 'o';
h{7}.MarkerFaceColor = col*0.8;
h{7}.MarkerEdgeColor = col*0.8;
h{7}.LineWidth = 0.5;
end