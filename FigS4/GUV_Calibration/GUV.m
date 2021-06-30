T=readtable('SurfaceDensity_GUVs.xlsx','Sheet','Fit');

c = table2array(T(1,2:6));
f = table2array(T(2,2:6));

plot(c,f,'o','MarkerSize',5,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[0,0,0])
xlabel('Dye Density (molecules/{\mu}m^2)')
ylabel('Fluorescence (a.u.)')

fi = fit(c',f','poly1');

h = gca;
h.FontSize = 16;
set(h,'box','off')
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
set(gcf,'color','w')

hold on

fun_l = @(x) fi.p1*x + fi.p2;
plot(0:1:max(c),fun_l(0:1:max(c)),'LineWidth',2,'Color','k')
