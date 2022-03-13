in = readmatrix('In_alt.csv');
in = in(:,2);
out = readmatrix('Out_alt.csv');
out = out(:,2);
div = in ./ out;
div = div / div(1);
%plot(div)
%t = 0:2.08:2.08*98;
plot(0:200,div(1:201))
plot(0:200,(in(1:201)-125.35)/505.52)
hold on
plot(0:200,(out(1:201)-125.35)/505.52)
%plot(t,div(1:99))
%plot(t,in(1:99)/452.7271)
hold on
%plot(t,out(1:99)/452.7271)
%ylim([300,1000])
xlim([0,200])
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off