close all
clear all
load('/Users/pjd17/Documents/Data/brain_distortion_struct-1.mat')%change this to match your file location
for k=1:length(brain_distortion_struct.normalized_kymograph);
    smoothed_data(:,k)=smooth(brain_distortion_struct.normalized_kymograph(:,k),3);%not the optimal smoothing but ok for quick and dirty
end
[Gx, Gy] = imgradientxy(smoothed_data,'sobel');%enhance edges with Sobel filter

half_ave_range=3;
the_range=120:255;
the_reference=1000:1020;
reference_vector=mean(Gy(the_range,1008:1010)');%this is the referene configuration everything is compared to
shift_profiles=zeros(size(Gy,2),2*length(the_range));
for k = (half_ave_range+1):(size(Gy,2)-half_ave_range)
    test_vector=mean(Gy(the_range,(k-half_ave_range):(k+half_ave_range))');%average of 3 time points
    [shift_curve,calc_shift,my_fit,the_goodness,the_tau,ix,iy]=myProfileStrainFit(reference_vector,test_vector,0);%function is at end of file
    shift_profiles(k,1:length(shift_curve))=shift_curve;
    fitted_tau(k)=the_tau;
    the_amp(k)=my_fit.a;
    goodness(k)=the_goodness.adjrsquare;
end
%%
figure(1)
subplot(211)
imagesc(smoothed_data)
title('raw data')
ylim([min(the_range) max(the_range)]);

subplot(212)
imagesc(Gy)
title('y gradient')
ylim([min(the_range) max(the_range)]);
figure(2)
subplot(211)
imagesc(shift_profiles')
ylim([0 the_range(end)-the_range(1)])
colorbar
title('shifts')
xlabel('time')
ylabel('pixels')
subplot(212)
hold off
plot(max(shift_profiles'))
hold on
plot(medfilt1(max(shift_profiles'),11))
title('maximum shift')
ylabel('shift, pixels')
xlabel('time')
legend('maximum shift','maximum shift, medfilt 11')

figure(3)
for k=(half_ave_range+1):(size(Gy,2)-half_ave_range)
    corrected_shift_profiles(k,:)=mean(shift_profiles((k-half_ave_range):(k+half_ave_range),:));
end
plot(max(corrected_shift_profiles'))
hold off
plot(medfilt1(max(corrected_shift_profiles'),11))
figure(111)
hold off
imagesc(smoothed_data)%Gy)
ylim([100 255])
hold on
shift_profile_index=[2:2:20 25 35 45 60 80 100 120 125];
for s=1:length(shift_profile_index)
plot(smoothdata((corrected_shift_profiles(:,shift_profile_index(s))'),'gaussian',3)+shift_profile_index(s)+min(the_range),'w', 'LineWidth',2)
end
plot(max(corrected_shift_profiles(:,1:20)')+min(the_range),'k', 'linewidth', 3)%plot the maximum shif tin the first 20 pixels

figure(4) %plots some shift curves
hold off
plot(mean(corrected_shift_profiles(530:540,:)))
hold on
plot(mean(corrected_shift_profiles(415:420,:)))
plot(mean(corrected_shift_profiles(1530:1540,:)))
plot(mean(corrected_shift_profiles(5795:5805,:)))
plot(mean(corrected_shift_profiles(1030:1040,:)))
legend('530:540','415:420','1530:1540','5795:5805','1030:1040')
xlim([0 200])
xlabel('location')
ylabel('displacement')


%%
function [shift_curve,calc_shift,my_fit,the_goodness,the_tau,ix,iy]=myProfileStrainFit(reference_profile,test_profile,plotting_on);
% reference_profile,test_profile - vectors of data.  this code assumes
% that the test profile is shifted to the right and the effects 'decay'
% away the higher the index
% shift curve - the displacement between the reference and test profile
% my_fit - fitting the strain profile w/ a decaying exponential curve
% the_tau - spatial memeory of the strain

[dist,ix,iy] = dtw((reference_profile),(test_profile),15,'euclidean');

shift_curve_raw=iy-ix;
% this 'shift curve' ends up being padded with zeros  on the left - this
% corrects this
calc_shift=find(shift_curve_raw>=max(shift_curve_raw),1); % find the earlised peak in the shift curve, and fit he decay from there
shift_curve=shift_curve_raw;%(calc_shift:end)%-min(shift_curve_raw(calc_shift:end));

% fit the displacement shift with an exponential function
[my_fit, the_goodness]=fit((1:length(shift_curve))',shift_curve,'exp1', 'StartPoint', [0 -1/10], 'Upper', [25 -1/100],'Lower', [0 -1]);%give reasonble start points,  tau ~ 10pixels
the_tau=1/my_fit.b;
%plot the shift_curve and fit for sanity checking
if (plotting_on==1)
    figure(111)
    subplot(311)
    hold off
    plot(reference_profile)
    hold on
    plot(test_profile)
    legend('reference','test')
    title('profiles')
    subplot(312)
    hold off
    plot(reference_profile(ix))
    hold on
    plot(test_profile(iy),'o')
    legend('reference','test')
    title('shifted profiles')
    subplot(313)
    hold off
    plot(shift_curve)
    hold on
    plot(my_fit)
    ylabel('shift')
    xlabel('pixels')
    figure
    plot(shift_curve_raw)
else
end
end