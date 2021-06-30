%% Figures
%From Nuclear Membrane Dynamics and Reassembly in Living Cells: 
%Targeting of an Inner Nuclear Membrane Protein in Interphase and Mitosis
frap = readmatrix('/Users/ikuz/Box Sync/ModelingBcLOV/MembraneFrap/12.28.20/Mem-7-ROI-2.csv');
t = frap(:,1);
f = frap(:,2);
a = frap(:,3);
a = a / max(a);
b = frap(:,4);
f = (f - b);
%f(22:end) = f(22:end) ./ a(22:end);
f = f/max(f);
subplot(1,2,1)
plot(t,f,'o','MarkerSize',5,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[0,0,0])
ylabel('Fluorescence (a.u.)')
xlabel('Time (s)')
h = gca;
h.FontSize = 18;
set(h,'box','off')
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
set(gcf,'color','w')
xlim([0,max(t)])
ylim([0.2,1.01])
ft1 = fittype('a*(1-exp(-b*x))');
ft2 = fittype('a*(1-exp(-b*x))+c*(1-exp(-d*x))');
ft1 = fittype('a*(1-(2.071^2*(2.071^2+4*pi*b*(x)).^-1).^(1/2))');
ft2 = fittype('a*(1-(2.071^2*(2.071^2+4*pi*b*(x)).^-1).^(1/2))+c*(1-(2.071^2*(2.071^2+4*pi*d*(x)).^-1).^(1/2))');
start=22;
f1 = fit(t(start+25:end)-t(start),f(start+25:end)-f(start),ft1,'StartPoint',[1,0.1],'Lower',[0,0,0,0],'Upper',[1e3,10,1e3,0.05])
f2 = fit(t(start:end)-t(start),f(start:end)-f(start),ft2,'StartPoint',[1,0.1,1,0.01],'Lower',[0,0,0,0],'Upper',[1e3,10,1e3,0.05]);
%fun1 = @(t) f1.a * (1-exp(-f1.b*t))+f(start);
%fun2 = @(t) f2.a * (1-exp(-f2.b*t)) + f2.c * (1-exp(-f2.d*t)) +f(start);
fun1 = @(t) f1.a*(1-(2.146^2*(2.146^2+4*pi*f1.b*(t)).^-1).^(1/2)) +f(start);
fun2 = @(t) f2.a*(1-(2.146^2*(2.146^2+4*pi*f2.b*(t)).^-1).^(1/2))+f2.c*(1-(2.146^2*(2.146^2+4*pi*f2.d*(t)).^-1).^(1/2)) +f(start);
hold on
plot(t(start:end),fun1(t(start:end)-t(start)),'k','LineWidth',2)

xi = linspace(-5,5,1000);
t = linspace(0,200,2001);
k_on_l = 31901;

k_off_l = 0.025;

k_off_d = 1/44.4;

k_off_p = 1/500000;
D = 8.932;
format long
Dm = 0;

k_on_d = 1126.1;

[aa,ab,ac,ad,ae,af,ag] = ndgrid(k_on_l,k_off_l,k_off_d,k_off_p,D,Dm,k_on_d);
v = [aa(:),ab(:),ac(:),ad(:),ae(:),af(:),ag(:)];

kd = v(1,2)/v(1,1)*1e6;
    s_c = .4;
    s_m = s_c/(s_c+kd);
    bleach_depth = 0.615;

    %ic_fun = @(x) [s_c;0;6000*s_m*(abs(x)>=1);0];
    ic_fun = @(x) [s_c*1e-6;0;6000*s_m*(abs(x)>=1.0355) + 6000*s_m*(1-bleach_depth)*(abs(x)<1.0355);6000*s_m*bleach_depth*(abs(x)<1.0355)];
    idx = find(abs(xi) < 1.0355);

    k_on_l = v(1,1);
    k_off_l = v(1,2);
    k_off_d = v(1,3);
    k_off_p = v(1,4);
    D = v(1,5);
    Dm = v(1,6);
    k_on_d = v(1,7);

    fun = @(x,t,u,dudx) rf_fun_alt5(x,t,u,dudx,k_on_l,k_off_l,D,Dm);
    sol_o = pdepe(0,fun,ic_fun,@bc_fun,xi,t);
    sol = squeeze(sol_o(:,:,3));

    model_res = mean(sol(:,idx),2)/6000/s_m;
top = mean(f(end-5:end));
bottom = model_res(1);
model_res = (model_res - bottom);
model_res = model_res ./ max(model_res);
model_res = model_res * (top - bottom) + bottom;

plot(t(11:end)+6.632,model_res(1:end-10),'k','LineWidth',2)
%plot(t(start:end),fun2(t(start:end)-t(start)),'k','LineWidth',2)
title('Membrane FRAP')
set(gca,'tickdir','out')
%%
%D via https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3731631/
frap = readmatrix('./CytosolicFRAP/Cyto-3-sample.csv');
t = frap(:,1);
f = frap(:,2);
a = frap(:,3);
a = a / max(a);
b = frap(:,4);
f = f/max(f);
%f = (f - b);
%f(22:end) = f(22:end) ./ a(22:end);
subplot(1,2,2)
plot(t,f,'o','MarkerSize',5,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[0,0,0])
ylabel('Fluorescence (a.u.)')
xlabel('Time (s)')
h = gca;
h.FontSize = 18;
set(h,'box','off')
h.TickLabelInterpreter = 'none';
h.LineWidth=2;
set(gcf,'color','w')
xlim([0,max(t)])
ylim([0.64,1.01])
ft1 = fittype('a*(1-exp(-b*x))+c');
ft2 = fittype('a*(1-exp(-b*x))+c*(1-exp(-d*x))');
start=27;
%f1 = fit(t(start:end)-t(start),f(start:end),ft1,'StartPoint',[1,11,0.6])
%f2 = fit(t(start:end)-t(start),f(start:end)-f(start),ft2,'StartPoint',[0.5,10,0.5,1],'Upper',[1,100,1,100])
ft1 = fittype('a*(1-(2.5^2*(2.5^2+4*pi*b*(x)).^-1).^(1/2))');
f1 = fit(t(start:end)-t(start),f(start:end)-f(start),ft1,'StartPoint',[1,0.1],'Lower',[0,0,0,0],'Upper',[0.4,8.932])

%fun1 = @(t) f1.a * (1-exp(-f1.b*t))+f1.c;
%fun2 = @(t) f2.a * (1-exp(-f2.b*t)) + f2.c * (1-exp(-f2.d*t)) +f(start);
fun1 = @(t) f1.a*(1-(2.5^2*(2.5^2+4*pi*f1.b*(t)).^-1).^(1/2)) +f(start);
hold on
plot(t(start:end),fun1(t(start:end)-t(start)),'k','LineWidth',2)
%plot(t(start:end),fun2(t(start:end)-t(start)),'k','LineWidth',2)
title('Cytosolic FRAP')
set(gca,'tickdir','out')