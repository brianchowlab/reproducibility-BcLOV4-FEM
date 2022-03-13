%T=readtable('./12.21.20/spatial-1-sim.xlsx','ReadVariableNames',0);
T=readtable('./12.21.20/Spatial-1-correct.xlsx','Sheet','1','ReadVariableNames',0);
%T=readtable('./12.22.20/Spatial-5.xlsx','ReadVariableNames',0);
%T=readtable('./12.22.20/Spatial-6.xlsx','ReadVariableNames',0);

t = table2array(T(2,1))-table2array(T(1,1));
x = table2array(T(:,2));
f = table2array(T(:,3:end));
t = t*(0:size(f,2)-1);
%%
%t1 = 0:0.3:0.3*121;
%t2 = 0.3*121+1:1: 0.3*121 + 200;
%t = [t1,t2];
clip_l = 1;%18;%24
clip_r = 1;%18;%8
pre = 20;%1
xt = x(clip_l:end-clip_r+1) - x(clip_l);
xi = 0:0.1:max(xt);
t_post = t(pre+1:end) - t(pre+1);
ti = 0:0.2:max(t_post);
f_pre = trimmean(f(clip_l:end-clip_r+1,1:pre),10,2);
f_post = f(clip_l:end-clip_r+1,(pre+1):end) - repmat(f_pre,[1,size(t,2)-pre]);
g = gausswin(5);
f_post = filter(g,1,f_post')';
f_post_i = zeros(size(xi,2),size(f_post,2));
for i = 1:size(f_post,2)
    f_post_i(:,i) = interp1(xt',f_post(:,i),xi);
end
f_post_ii = zeros(size(xi,2),size(ti,2));
for i = 1:size(xi,2)
    f_post_ii(i,:) = interp1(t_post',f_post_i(i,:),ti);
end
f_post_ii = f_post_ii / prctile(f_post_ii(:),99.7);

r=301;
h=heatmap(f_post_ii(:,1:r))
set(h.NodeChildren(3), 'XTickLabelRotation', 0);
h.GridVisible = 'off';
h.Colormap = parula;
%h.ColorLimits = [0,prctile(f_post(:),99.9)];
h.ColorLimits = [0,1];
h.XLabel = 'Time (s)';
h.YLabel = 'Position ({\mu}m)';
cdl = h.XDisplayLabels; 
%xd = num2cell(t_post);
xd = repmat(NaN,size(cdl,1), size(cdl,2));
xd(1:20:end) = string(ti(1:20:r))';
h.XDisplayLabels = xd;
cdl = h.YDisplayLabels; 
%xd = num2cell(t_post);
xd = repmat(NaN,size(cdl,1), size(cdl,2));
xd(1:20:end) = string(xi(1:20:end))';
h.YDisplayLabels = xd;
h = gca;
h.FontSize = 18;
samp = 7;
%%
%Bleach -0.5627
temp = f_post_ii(:,:);
temp = temp - temp(:,1);
for i = 1:size(temp,2)
   temp(:,i) = temp(:,i) + 0.6 * i/  size(temp,2);
end
temp = temp(1:end-1,1:end-1);
temp = splitapply(@mean,temp,ceil((1:size(temp))/2)');
temp = splitapply(@mean,temp',ceil((1:size(temp'))/2)')';
temp = temp ./ max(temp(:,21));

%1
plot(xi(1:2:end-1),temp(:,3))
hold on
%4
plot(xi(1:2:end-1),temp(:,10))
%20
plot(xi(1:2:end-1),temp(:,75))
xlim([0,15])
ylim([0,1])
xticks([0 5 10 15])
%%
%figure
%plot(xt,splitapply(@mean,f_post',ceil((1:size(f_post'))'/samp))')
%legend
f_avg  = splitapply(@mean,f_post',ceil((1:size(f_post'))'/samp))';
%splitapply(@mean,f_post',ceil((1:size(f_post'))'/3))'
wid = [];
cen = [];
amp = [];
for i = 1:size(f_avg,2)
    f_g = fit(xt,f_post(:,i),'gauss1','Lower',[0,0,0],'Upper',[1000,1000,1000]);
    wid = [wid,f_g.c1];
    cen = [cen,f_g.b1];
    amp = [amp,f_g.a1];
end

figure
subplot(1,3,1)
plot(t_post(1:samp:end),wid)
ylabel('Gaussian STD')

subplot(1,3,2)
plot(t_post(1:samp:end),cen)
ylabel('Gaussian Mean')

subplot(1,3,3)
plot(t_post(1:samp:end),amp)
legend('STD','Mean','Amp')
xlabel('Time (s)');
ylabel('Gaussian Amp')
set(gca,'FontSize',18)
%%
f_post_scaled = f_post - min(f_post);
f_post_scaled = f_post_scaled./ max(f_post);
h=heatmap(f_post_scaled)
h.GridVisible = 'off';
h.Colormap = jet;
h.ColorLimits = [0,1];
h.XLabel = 'Time (s)';
cdl = h.XDisplayLabels; 
%xd = num2cell(t_post);
xd = repmat(NaN,size(cdl,1), size(cdl,2));
xd(1:20:end) = string(t_post(1:20:end))';
h.XDisplayLabels = xd;
cdl = h.YDisplayLabels; 
%xd = num2cell(t_post);
xd = repmat(NaN,size(cdl,1), size(cdl,2));
xd(1:20:end) = string(xt(1:20:end))';
h.YDisplayLabels = xd;
h = gca;
h.FontSize = 16;
figure
plot(splitapply(@mean,f_post_scaled',ceil((1:size(f_post_scaled'))'/100))')