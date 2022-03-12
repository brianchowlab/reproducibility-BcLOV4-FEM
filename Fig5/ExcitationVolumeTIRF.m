function [il_3D,coords] = ExcitationVolumeTIRF(mask_il,z,params)

    %Define a 1um x 1um ROI that is rastered each 100 nm
    mask_il = zeros(100,100);
    mask_il(46:55,46:55) = 1;
    mask = mask_il';
    res = 0.1;
    sampling_density = 0.01;
    alpha = asind(params.NA/params.immersion_n);
    
    within_aperture = @(x,y,z,xo,yo,zo) sqrt((x-xo).^2 + (y-yo).^2) ./ abs(z-zo) < tand(alpha);
    
    %TIRF specific stuff%
    d =  @(lambda,n1,n2,theta) lambda/4/pi * (n2^2 * sind(theta)^2 - n1^2)^(-1/2);
    wave_z = @(lambda,n1,n2,theta,z) exp(-z/d(lambda,n1,n2,theta));
    lambda = .445;
    n1 = 1.33;
    n2 = 1.52;
    theta = 79;
    w0 = sqrt(2) * 0.325 * lambda / (sqrt(2) * 1.4^0.91);
    profile = @(r) exp(-2*r.^2./(w0.^2));
    
    [y,x] = find(mask==1) ;
    x0 = min(x); y0 = min(y) ;
    x1 = max(x); y1 = max(y) ;
    on = mask(y0:y1,x0:x1);
    
    %Define raster scanning account for difference between sampling density
    %and raster density
    on_raster = zeros(101,101);
    on_raster(6:10:end,6:10:end) = 1;
    
    left = tand(alpha) * min(z) - size(on,2)/2 * res;
    right = tand(alpha) * max(z) + size(on,2)/2 * res;
    X = 0:sampling_density:right;
    X = [fliplr(-X),X(2:end)];
    Y = X;

    [X,Y,Z] = meshgrid(X,Y,z);
    coords.X = X;
    coords.Y = Y;
    coords.Z = Z;
    R = sqrt(X.^2 + Y.^2);
    
    %Set irradiance here
    irradiance_integ = profile(R).*wave_z(lambda,n1,n2,theta,Z);
    irradiance_integ(find(Z<0)) = 0;
    C = zeros(size(irradiance_integ));
    for i = 1:size(z,2)
        C(:,:,i) = conv2(irradiance_integ(:,:,i),on_raster,'same');
    end
    irradiance_integ = C / max(C(:));
    il_3D = irradiance_integ;
end