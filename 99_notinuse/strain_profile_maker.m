% strain_profile_maker
clear all
close all
x=1:100;
noise_level=.1;
reference_vector=max(0,sign(sin(2*pi*x/25)-.7))% this is the profile that we try to match, this generates periodic spikes

spatial_tau=25;
spatial_shift=5;
for n=1:length(x)
    test_x(n)=n+spatial_shift*exp(-(n)/spatial_tau);
end
test_vector=max(0,interp1(test_x,reference_vector,x,'nearest')+noise_level*randn(size(reference_vector)));%=max(0,spline(test_x,reference_vector,reference_x));
