function [shift_curve,calc_shift,my_fit,the_tau]=ProfileStrainFit(reference_profile,test_profile);
% reference_profile,test_profile - vectors of data.  this code assumes
% that the test profile is shifted to the right and the effects 'decay'
% away the higher the index
% shift curve - the displacement between the reference and test profile 
% my_fit - fitting the strain profile w/ a decaying exponential curve
% the_tau - spatial memeory of the strain
figure(111)
subplot(311)
hold off
plot(reference_profile)
hold on
plot(test_profile)
legend('reference','test')
title('profiles')
[dist,ix,iy] = dtw(reference_profile,test_profile,'absolute');

subplot(312)
hold off
plot(reference_profile(ix))
hold on
plot(test_profile(iy),'o')
legend('reference','test')
title('shifted profiles')
shift_curve_raw=iy-ix;
% this 'shift curve' ends up being padded with zeros  on the left - this
% corrects this
calc_shift=find(shift_curve_raw>=max(shift_curve_raw),1); % find the earlised peak in the shift curve, and fit he decay from there
shift_curve=shift_curve_raw(calc_shift:end)-min(shift_curve_raw(calc_shift:end)); 

% fit the displacement shift with an exponential function
my_fit=fit((1:length(shift_curve))',shift_curve,'exp1', 'StartPoint', [shift_curve(1) -1/10])%give reasonble start points,  tau ~ 10pixels
the_tau=1/my_fit.b;

%plot the shift_curve and fit for sanity checking
subplot(313)
hold off
plot(shift_curve)
hold on
plot(my_fit)
ylabel('shift')
xlabel('pixels')
figure
plot(shift_curve_raw)