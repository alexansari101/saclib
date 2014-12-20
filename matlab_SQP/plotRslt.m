global dt tf
lnwidth = 2;

figure(1)
subplot(1,2,1)
plot([0+dt:dt:tf],x(1:3:end),'b', 'LineWidth',lnwidth)
hold on
plot([0+dt:dt:tf],x(2:3:end),'Color', [.4 0 .6], 'LineWidth',lnwidth)
title('theta and theta_dot vs time')
xlabel('time')
ylabel('theta (blue) & theta_dot (magenta)')
hold on
subplot(1,2,2)
plot([0+dt:dt:tf],x(3:3:end),'b', 'LineWidth',lnwidth)
title('acceleration control vs time')
xlabel('time')
ylabel('cart acceleration')