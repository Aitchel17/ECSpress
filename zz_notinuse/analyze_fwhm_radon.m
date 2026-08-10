function [lineprofile,sinogram] = analyze_fwhm_radon(maskstack,vertices,background)
    arguments
        maskstack (:,:,:) {mustBeNumeric}
        vertices (:,:) {mustBeNumeric}
        background = 0 
    end
    vertices = vertices(1:2,1:2);
    lvec = vertices(2,:) - vertices(1,:);
    ltheta = rad2deg(atan2(lvec(2),lvec(1)));
    if ltheta<0
        ltheta = 360+ltheta;
    end
    maskstack(isnan(maskstack)) = background; % convert nan value to 0 for projection analysis
    sinogram = radon(maskstack(:,:,2380),0:1:360);
    lineprofile = sinogram(:,floor(ltheta)+1);
    figure(Name='lineprofile',NumberTitle='off')
    plot(lineprofile)
end