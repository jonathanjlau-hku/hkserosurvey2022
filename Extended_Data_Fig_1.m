clear all;
close all;

% Purpose: Plot age effect factors under Extended Data Figure 1
% Note: Figure will be created in a separate Matlab window named "Figure 1"

load("Extended_data_Figure_1.mat")


figure(1);
clf;
hold on;

h_7delay = fill([age_range fliplr(age_range)], [dat_7delay(1,:) fliplr(dat_7delay(end,:))], [0.902 0.294 0.208]);
transparency = 0.2;
set(h_7delay, "FaceAlpha", transparency, "lineStyle", "none");

h_14delay = fill([age_range fliplr(age_range)], [dat_14delay(1,:) fliplr(dat_14delay(end,:))], [0.302 0.733 1]);
transparency = 0.2;
set(h_14delay, "FaceAlpha", transparency, "lineStyle", "none");

h_21delay = fill([age_range fliplr(age_range)], [dat_21delay(1,:) fliplr(dat_21delay(end,:))], "g");
transparency = 0.2;
set(h_21delay, "FaceAlpha", transparency, "lineStyle", "none");


plot(age_range, dat_7delay(round(size(dat_7delay,1)/2),:),"Color",[0.902 0.294 0.208]);
plot(age_range, dat_14delay(round(size(dat_14delay,1)/2),:),"Color",[0.302 0.733 1]);
plot(age_range, dat_21delay(round(size(dat_21delay,1)/2),:),"Color","g");


axis tight;
xlabel("Age");
ylabel("Relative FOI vs 35 year old subject");
legend("7 day delay","14 day delay","21 day delay");

exportgraphics(gcf,'extended_data_fig_1.tiff','Resolution',300);