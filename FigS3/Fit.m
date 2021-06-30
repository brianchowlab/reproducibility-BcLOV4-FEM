%Fittable: 1,2,3,4,5,6,7,10,11,12,13,14,18,19,20,22,23,25
%Fittable 2.28.21-Data
%1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,23,24

konv = [30000:200:40000];


koffv = [0.03:0.001:0.04];


%32000,0.038
[X,Y] = meshgrid(konv,koffv);
params = [X(:),Y(:)];
errors = [];

for k=1:size(params,1)
    k
    c = [1,2,3,4,5,6,7,8,10,11,12,13,15,16,17,18,19,23,24,26,28,29,30];%25
    %c = 1;
    %c = [1:8,10:30];
    a = readmatrix('~/Documents/Analysis.xlsx','Sheet','2.28.21-100%');
    %a = readmatrix('~/Documents/Analysis.xlsx','Sheet','2.28.21-SinglePulse');

    V = a(:,9);
    SA = a(:,10);
    frac = a(:,11);
    C = [];
    CAll = [];
    res = zeros(size(c,2),5001)';
    
    for i = 1:size(c,2)
        %i
        filename = ['./2.28.21-Data/Cell-',num2str(c(i)),'-Cyto-1.csv'];

        frac_cytoplasm = frac(c(i));

        if c(i) < 16
            d = readmatrix(filename);
            t = 0:50;
            a = d(1:51,2);
            C1 = (a-104) / 452.7271;
            filename = ['./2.28.21-Data/Cell-',num2str(c(i)),'-Cyto-2.csv'];
            d = readmatrix(filename);
            t = 0:50;
            a = d(1:51,2);
            C2 = (a-104) / 452.7271;
            filename = ['./2.28.21-Data/Cell-',num2str(c(i)),'-Cyto-3.csv'];
            d = readmatrix(filename);
            t = 0:50;
            a = d(1:51,2);
            C0 = (max(a)-104)/452.7271*1e-6;
            C3 = (a-104) / 452.7271;
            CA = mean([C1,C2,C3]')';
        else
            d = readmatrix(filename);
            t = 0:50;
            a = d(1:51,2);
            C0 = (max(a)-104)/452.7271*1e-6;
            CA = (a-104) / 452.7271;
        end
        C = [C,C0];
        CAll = [CAll,CA];



        [y] = zerodim_model(params(k,1),params(k,2),20,C0,1,0.1,SA(c(i))/V(c(i))/frac(c(i))/602.2e6);
        res(:,i) = y(:,1) + y(:,2);
        %out = fit(t',m/median(m(40:end)),fun,'StartPoint',1,'Lower',0);
        %kt = [kt,out.k];

    end
end
res = res(1:100:end,:) * 1e6;
plot(CAll)
hold on
plot(res)
MSE = (CAll - res).^2;
MSE_s = sum(MSE);
%Sort 19,20,5,23,22,9,16,1,21,15,17,10,6,4,11,8,18,12,13,7,14,2,3