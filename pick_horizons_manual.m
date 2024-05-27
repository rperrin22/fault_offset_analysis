function [pts] = pick_horizons_manual(index,Df,trace_vec,time_vec,Data,H076,H090,mat_min,mat_max,fault_xvec_1,fault_xvec_2,fault_yvec)

figure('units','normalized','outerposition',[0 0 1 1]);

imagesc(trace_vec,time_vec,Data)
colormap bone
set(gca,'clim',[-5,5])
hold on
plot([Df.pt1_trace(index), Df.pt2_trace(index)],[Df.pt1_time(index), Df.pt2_time(index)],'r','linewidth',2)
plot(fault_xvec_1,fault_yvec,'k:','linewidth',1.5)
plot(fault_xvec_2,fault_yvec,'b:','Linewidth',1.5)
plot(H076.trace,H076.time,'c','linewidth',1.5)
plot(H090.trace,H090.time,'g','linewidth',1.5)
xlabel('CDP Number')
ylabel('Time (s)')
grid on
set(gca,'xlim',[mat_min,mat_max])
title(sprintf('Fault %d',index))
set(gca,'ylim',[min(fault_yvec)*0.98,max(fault_yvec)*1.02])

% Initialize array to store calculated lag
pts = [];

% Initialize arrays to store coordinates
x_coords = [];
y_coords = [];

% Loop to interactively select points
while true
    % Wait for user to click on the image
    [x, y, button] = ginput(1);

    % Exit loop if right-click or Enter is pressed
    if button == 3 || button == 13
        break;
    end

    % Store coordinates
    x_coords = [x_coords, x];
    y_coords = [y_coords, y];

    % Plot the selected point
    hold on;
    plot(x, y, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');

    % If two points are selected, plot a line between them
    if length(x_coords) == 2
        plot(x_coords, y_coords, 'b-', 'LineWidth', 2);
        pts(end+1,1) = diff(y_coords);
        pts(end,2) = mean(y_coords);
        x_coords = [];
        y_coords = [];
    end
end
saveas(gcf,sprintf('L7_f%d_horizoncheck.png',index));

end