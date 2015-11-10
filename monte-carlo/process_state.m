clear
close

data = importdata('1/vicon_state_24171730.txt');
pos_x = data(:,2)/1000;
pos_y = data(:,3)/1000;
pos_z = data(:,4)/1000 + 17;

vel_x = data(:,5);
vel_y = data(:,6);
vel_z = data(:,7);

Time = 1:length(pos_x);
Time = Time*0.01612903;
boundary = 0.2*pos_z;

hold on
plot(Time, pos_x)
plot(Time, pos_y)
plot(Time, pos_z)
plot(Time, boundary, 'r-.')
plot(Time, -boundary, 'r-.')
xlabel('Time [s]')
ylabel('Position [m]')
legend('x', 'y', 'z')

xlim([0 60])
grid on



figure(2)
hold on
plot(Time, vel_x)
plot(Time, vel_y)
plot(Time, vel_z)
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('V_x', 'V_y', 'V_z')

xlim([0 60])