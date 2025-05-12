clear; clc;

%% Script 1: Actual Data and Fitted Curve
dt_actual = 10;
t_actual = 0:dt_actual:50*dt_actual;
T_actual = [ 28.0;
    28.0;
    28.0;
    28.0;
    28.1;
    28.1;
    28.1;
    28.2;
    28.2;
    28.2;
    28.2;
    28.2;
    28.2;
    28.3;
    28.3;
    28.3;
    28.3;
    28.4;
    28.5;
    28.6;
    28.8;
    29.1;
    29.1;
    29.2;
    29.3;
    29.5;
    29.5;
    29.6;
    29.7;
    29.7;
    29.8;
    29.8;
    29.9;
    29.9;
    29.9;
    29.9;
    30.0;
    30.1;
    30.1;
    30.1;
    30.1;
    30.2;
    30.2;
    30.3;
    30.3;
    30.3;
    30.3;
    30.3;
    30.4;
    30.5;
    30.6];

p_actual = polyfit(t_actual, T_actual, 1);
z_actual = polyval(p_actual, t_actual);

%% Script 2: Simulation
T_b0 = 90.5; % Initial temperature of each steel ball (°C)
T_w0 = 28; % Initial temperature of water (°C)
m_b = 3.98*1e-3; % Mass of each steel ball (kg)
m_w = 350*1e-3; % Mass of water in the bucket (kg)
c_b = 490; % Specific heat capacity of steel (J/kg°C)
c_w = 4186; % Specific heat capacity of water (J/kg°C)
r_b = 0.005; % Radius of a steel ball (m)
A_b = 4*pi*r_b^2; % Surface area of the steel ball (m²)
rho_b = m_b/(4*pi*(r_b^3)/3); % Density of steel (kg/m³)
V_b = (4/3)*pi*r_b^3; % Volume of the steel ball (m³)
h = 329.74; % average Convective heat transfer coefficient (W/m²°C)

dt_sim = 10; % time duration between dropping two consecutive balls
t_final = dt_sim*50; % Total simulation time (s)
time_sim = 0:dt_sim:t_final; % Time vector for simulation
T_w = zeros(size(time_sim)); 
T_w(1) = T_w0;
balls = [];
Qtot = 0;

for t_idx = 2:length(time_sim)
    current_time = time_sim(t_idx);

    if mod(current_time, dt_sim) == 0
        balls = [balls; T_b0];
    end

    Q_total = 0;

    for i = 1:length(balls)
        T_ball = balls(i);
        T_new = T_w(t_idx-1) + (T_ball - T_w(t_idx-1)) * exp(- (h*A_b)/(rho_b*c_b*V_b) * dt_sim);
        Q_i = m_b * c_b * (T_ball - T_new);
        Q_total = Q_total + Q_i;
        balls(i) = T_new;
    end

    Qtot = Qtot + Q_total;
    delta_Tw = Q_total / (m_w * c_w);
    T_w(t_idx) = T_w(t_idx-1) + delta_Tw;
end

Tavg_sim = Qtot / (t_final * A_b * h * 50);
Tavg_actual = (m_w * c_w * 2.5) / (500 * 50 * A_b * h);

%% Plot both first figures on same plot
figure;
hold on;
plot(t_actual, T_actual, '-x', 'Color', '#7E2F8E', 'LineWidth', 2);
plot(t_actual, z_actual, 'm', 'LineWidth', 2);
plot(time_sim, T_w, 'r-x', 'LineWidth', 2);
xlabel('Time (s)')
ylabel('Water Temperature (°C)')
title('Actual and Simulated Water Temperature Rise')
legend('Actual data', 'Fitted curve', 'Simulated data', 'Location', 'best');
grid on;
hold off;

%% Second simulation for ball temperatures
figure;
hold on;

balls_data = {};
T_w = zeros(size(time_sim)); 
T_w(1) = T_w0;
balls = [];
ball_times = [];

for t_idx = 2:length(time_sim)
    current_time = time_sim(t_idx);

    if mod(current_time, dt_sim) == 0
        balls = [balls; T_b0];
        balls_data{end+1} = T_b0;
        ball_times(end+1) = t_idx;
    end

    Q_total = 0;

    for i = 1:length(balls)
        T_ball = balls(i);
        T_new = T_w(t_idx-1) + (T_ball - T_w(t_idx-1)) * exp(- (h*A_b)/(rho_b*c_b*V_b) * dt_sim);
        Q_i = m_b * c_b * (T_ball - T_new);
        Q_total = Q_total + Q_i;
        balls(i) = T_new;

        if abs(T_new - T_w(t_idx-1)) > 0.5
            balls_data{i}(end+1) = T_new;
        end
    end

    delta_Tw = Q_total / (m_w * c_w);
    T_w(t_idx) = T_w(t_idx-1) + delta_Tw;
end

p1 = plot(NaN, NaN, 'b');
for i = 1:length(balls_data)
    bt = balls_data{i};
    t_start = ball_times(i);
    t_end = t_start + length(bt) - 1;

    if t_end > length(time_sim)
        t_end = length(time_sim);
        bt = bt(1:t_end - t_start + 1);
    end

    plot(time_sim(t_start:t_end), bt, 'b');
end

p2 = plot(time_sim, T_w, 'r-x', 'LineWidth', 2);

xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Temperature of Water and Steel Balls Over Time');
legend([p1 p2], {'Steel Balls', 'Water'}, 'Location', 'best');
grid on;
hold off;
