clear
close

data = importdata('accel_correct.txt');
accel_x_true = data(:,1);

Time_true = 1:length(accel_x_true);
Time_true = Time_true*0.01612903;


data = importdata('accel_wrong.txt');
accel_x_w = data(:,1);

Time_w = 1:length(accel_x_w);
Time_w = Time_w*0.01612903;


figure(3)
hold on
plot(Time_true, accel_x_true,'r.-');
plot(Time_w, accel_x_w);
xlim([0 30])
xlabel('Time [s]','FontSize', 20)
ylabel('Acceleration Commands in x [m/s^2]','FontSize', 20)
h_legend=legend('with restart','without restart');
set(h_legend,'FontSize',14);