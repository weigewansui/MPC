clear
close

data = importdata('accel_28234927.txt');


accel_x = data(:,1);
accel_y = data(:,2);
accel_z = data(:,3);

Time = 1:length(accel_x);
Time = Time*0.01612903;
figure(3)
hold on
plot(Time, accel_x);
plot(Time, accel_y);
plot(Time, accel_z);
xlim([0 60])
xlabel('Time [s]','FontSize', 20)
ylabel('Acceleration Commands [m/s^2]','FontSize', 20)
h_legend=legend('x direction','y direction','z direction');
set(h_legend,'FontSize',14);