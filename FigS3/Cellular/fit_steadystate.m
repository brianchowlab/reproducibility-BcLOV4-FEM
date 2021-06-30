function [error,y,m0,Cf] = fit_steadystate(data,KD_l,S,KD_d,SA,V,C0)
    m0 = S * C0 ./ (KD_d + C0);
    C0 = C0*1e-6;
    KD_l = KD_l * .15e-6;
    a = -SA./V/602.2e6;
    b = C0 + S * SA./V/602.2e6 + m0.*SA./V/602.2e6+KD_l;
    c =  -S*C0 - S*m0.*SA./V/602.2e6;
    y = (-b + sqrt(b.^2 - 4.*a.*c))./(2*a);
    
    %max_dark_mem = C0*1e6 * 452.7271 / 1.8448;

    Cf = (1-(y-m0).*SA./(C0*602.2e6)./V).*C0*1e6;
    %check = sum(((S-y).*(1-(y-m0).*SA./(C0*602.2e6)./V).*C0 - KD_l*y).^2)
    error = sum(C0.*(data + m0 - y).^2);
end