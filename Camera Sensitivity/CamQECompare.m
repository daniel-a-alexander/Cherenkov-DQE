
%%
colors3 = linspecer(3);


%%
gen3_wl = gen3(:,1);
gen3_blue_qe = gen3(:,2);
gen3_red_qe = gen3(:,3);
gen3_blue_sens = gen3(:,5);
gen3_red_sens = gen3(:,6);
wl = gen3_wl;


gen2_wl = gen2(:,1);
gen2_sens = gen2(:,2);
gen2_sens = interp1(gen2_wl,gen2_sens,wl);

gen2_qe = 100*gen2(:,4); % percentage
gen2_qe = interp1(gen2_wl,gen2_qe,wl);

%%

load('qesens.mat');

%% Plot QE

figure;

ax1 = subplot(2,1,1);
hold on;
plot(wl, gen3_blue_sens, '-.', 'Color', colors3(1,:), 'Linewidth',1.5);
plot(wl, gen3_red_sens, '-.', 'Color', colors3(2,:), 'Linewidth',1.5);
plot(wl, gen2_sens, '-.', 'Color', colors3(3,:), 'Linewidth',1.5);
ylabel('Sensitivty [mA/W]', 'FontSize', 16);
set(gca,'XTickLabel',[]);
axis([400,800,0,200]);
ax1.FontSize = 16;

p(1) = area(0,0, 'FaceColor', colors3(1,:));
p(2) = area(0,0, 'FaceColor', colors3(2,:));
p(3) = area(0,0, 'FaceColor', colors3(3,:));
legend(p, {'Gen3 Blue','Gen3 Red','Gen2+'}, 'Location', 'NorthWest');
grid on;

ax2 = subplot(2,1,2);
hold on;
plot(wl, gen3_blue_qe, '-', 'Color', colors3(1,:), 'Linewidth',1.5);
plot(wl, gen3_red_qe, '-', 'Color', colors3(2,:), 'Linewidth',1.5);
plot(wl, gen2_qe, '-', 'Color', colors3(3,:), 'Linewidth',1.5);
ylabel('Q.E. [%]', 'FontSize', 16);
xlabel('Wavelength [nm]', 'FontSize', 16);
axis([400,800,0,40]);
ax2.FontSize = 16;
grid on;





