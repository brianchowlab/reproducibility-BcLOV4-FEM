load('TIRF_m_intrp.mat')

profile = squeeze(m_intrp(:,301,:) - m_intrp(:,301,1));

x = -30:0.1:30;
clipped_x = x(101:502);
clipped = profile(101:502,:);

wid = [];
wid_m = [];
cen = [];
amp=[];
max_height = [];
resamp = min(clipped_x):0.01:max(clipped_x);

for i=1:size(profile,2)
    f = clipped(:,i);
    %halfmax = (min(f) +max(f))/2;
    fr = interp1(clipped_x,f,resamp);
    halfmax = max(f)/2;
    fifthmax = 33.31/10;
    idx1 = find(fr>= halfmax,1,'first');
    idx2 = find(fr>= halfmax,1,'last');
    fwhm  = resamp(idx2) - resamp(idx1);
    idx1_m = find(fr>= fifthmax,1,'first');
    idx2_m = find(fr>= fifthmax,1,'last');
    fwhm_m  = resamp(idx2_m) - resamp(idx1_m);
    if size(fwhm_m,2)==0
       fwhm_m = NaN; 
    end

    f_g = fit(clipped_x',f,'gauss1','Start',[300,0,5],'Lower',[max(f)*1.2,-1,0],'Upper',[1000,1,100]);
    wid = [wid,fwhm / 2.355];
    wid_m = [wid_m,fwhm_m];
    cen = [cen,f_g.b1];
    amp = [amp,f_g.a1];
    max_height = [max_height,max(f)];
    if mod(i,60) == 0
        plot(f_g,clipped_x,f)
    end
end
wid(1) = 0;%0.75/sqrt(2*log(2));
wid_vcell = wid;

f = figure
t = 0:0.1:60;
plot(t,max_height,'k','LineWidth',2)
xlabel('Time (s)')
ylabel('Amplitude (molecules/um^2)')
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
xlim([0,60])

hold on
load('vcell_TIRF_result.mat')
t = 0:0.5:100;
plot(t,max_height_vcell,'--k','LineWidth',2)
f.Position = [1000 500 700 600];
ax = gca;
ax.TickLength = [0.025, 0.025];
