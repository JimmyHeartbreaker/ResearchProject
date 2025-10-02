function plot_control(T,N,U_opt)

    time = linspace(0,T,N);
    U_opt_reshaped = reshape(U_opt,2,N);
    figure;
    plot(time, U_opt_reshaped(1,:), 'r', 'LineWidth',2); hold on;
    plot(time, U_opt_reshaped(2,:), 'b', 'LineWidth',2);
    xlabel('Time [s]');
    ylabel('Control (Thrust)');
    legend('Ux','Uy');
    title('SCP Convexified 2D Minimum-Fuel Control');
    grid on;
end