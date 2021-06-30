%mCherry Calibration
%Right scope, 5% Excitation (?)
%New is right scope 10% excitation
%h = subplot(2,2,1);
x = [0,148.25,296.5,593];%nM
c = [104.731,1637.452,3385.455,7291.853];
std = [0.083862983,41.87679552,33.99020306,31.81235372];





x = [0,39.25,78.5,157,314];
c = [111.7883,1043.502,2374.978,6581.292,11929.428];
std = [0.123581282,5.468664279,53.41028375,123.7847072,145.822707];

%Theoretical cyto conc from the fluor
c_t  =  (c-111)/475;
ratio = c_t ./ (x/1000);
ratio = mean(ratio(2:end));
x_t = (x/1000) * ratio;


f = @(x) 565.6/1000 *x;

plot(x_t,c,'o','MarkerSize',5,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[0,0,0])

hold on
errorbar(x_t,c,3*std,'ok')
plot(x_t,f(x_t*1000),'Color','k','LineWidth',2);
h = gca;
h.FontSize = 16;
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
xlabel('Concentration (uM)')
ylabel('Fluorescence (a.u.)')
title('mCherry Calibration')
ylim([0,12500])
box off
set(gca,'tickdir','out')
%%
figure
%Lucifer Yellow Calibration
h = subplot(2,2,2);
x = [0,1250,2500,5000];%nM
c=[108.73,2363.86,4461.716,8917.78];
%f = @(x) 1.7599*x+108.73;
f = @(x) 0.91796 *x + 111;

plot(x,c,'o','MarkerSize',5,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[0,0,0])
hold on
plot(x,f(x),'Color','k','LineWidth',2);
h = gca;
h.FontSize = 16;
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
xlabel('Concentration (nM)')
ylabel('Fluorescence (a.u.)')
title('Lucifer Yellow Calibration')
xlim([0,5100])
box off
set(gca,'tickdir','out')

%Fluorescence over time
cd = readmatrix('/Users/ikuz/Box Sync/ModelingBcLOV/CytosolicConcentration/25um_sample_100ms_exposure_10%_cyan.csv');
scaling_ratio= 1.7599/0.0754;
h = subplot(2,2,3);
t = cd(:,1);
plot(t,(cd(:,2)-111)*scaling_ratio,'k','LineWidth',2)
h = gca;
h.FontSize = 16;
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
xlabel('Time (s)')
ylabel('Fluorescence (a.u.)')
title('Cellular Fluorescence')
xlim([0,1620])
ylim([0,600])

set(gcf,'color','w')
box off
set(gca,'tickdir','out')
%%
figure

bc_c = readmatrix('/Users/ikuz/Box Sync/ModelingBcLOV/CytosolicConcentration/ConfocalConcentrations.csv');
bc_c(bc_c > 17) = [];
%xlim([0,25])
%xh = subplot(2,2,4);
edges = 0:15;
histogram(bc_c,edges,'FaceColor','none','EdgeColor','k','LineWidth',2)
h = gca;
h.FontSize = 16;
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
xlabel('Cell Concentration ({\mu}M)')
ylabel('Frequency')
title('Concentration Distribution')

box off
set(gca,'tickdir','out')
