%% OFF FIT FINAL FOR PAPER
%Fittable: 1,2,4,6,7,8,9,10,11,12
figure
i = [1,2,4,6,7,8,9,10,11,12];
data = [];
cis = [];
cyto = [];
for j = 1:size(i,2)
    filename = ['Cell-',num2str(i(j)),'-Cyto.csv'];
    a = readmatrix('../Analysis.xlsx','Sheet','Off');
    V = a(:,9);
    SA = a(:,10);
    frac = a(:,11);
    frac = 0.33*ones(size(frac));
    d = readmatrix(filename);
    t = 0:(size(d,1)-1);
    f = d(:,2);
    c = (f - 104)/452.7271;
    cyto = [cyto,c(1:200)];
    ci = max(c);
    cis = [cis,ci];
    m = ci - c;
    m = m*602.2*V(i(j))/SA(i(j))*frac(i(j));
    t = t(1:200);
    m = m(1:200);
    data = [data,m];
end

cyto = cyto - cyto(1,:);
cyto = cyto ./ max(cyto(end-5:end,:));
ci = bootci(1000,@mean,cyto')';
temp = cyto;
cyto = mean(cyto') - min(ci(:,1));
temp = temp - min(ci(:,1));
ci = ci - min(ci(:,1));
scale = max(cyto(:,end-5:end));
cyto = cyto ./ scale;
ci = ci ./ scale;
patch([t fliplr(t)], [ci(:,1)' fliplr(ci(:,2)')], 'b','FaceAlpha',0.2,'EdgeAlpha',0)
hold on
plot(t,cyto,'k','LineWidth',1)

fun = @(off,to,x) (1-exp(-(x-off)/to));
off = fit(t(1:end-5)',cyto(6:end)',fun,'StartPoint',[5.5,44])
%plot(off,t(1:end-39)',cyto(40:end)')
plot(t,fun(5.65,off.to,t),'--r');
ylim([0,1.05])
set(gca,'FontSize',24)
set(gca,'tickdir','out')
set(gca,'linew',1.5)
xlabel('Time (s)')
ylabel('Scaled Fluorescence (a.u.)')

c=1:size(temp,2);
store = [];
for i = 1:100
    temp2 = datasample(c,max(c));
    temp2 = mean(temp(:,temp2)')';
    temp2 = temp2 ./ max(temp2(end-5:end,:));
    off = fit(t(1:end-5)',temp2(6:end),fun,'StartPoint',[5.5,44]);
    store = [store;off.to];
end
quantile(store,0.025)
quantile(store,0.975)