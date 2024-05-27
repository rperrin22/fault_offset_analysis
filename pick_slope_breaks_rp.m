function [xtemp,ytemp,slopes,slopelocs] = pick_slope_breaks_rp(ind_1,ind_x_1,ind_2,ind_x_2,ind_3,ind_x_3,pts)

figure;
ind_plot_1 = ind_1*1000;
plot(ind_plot_1,ind_x_1)
xlabel('Lag (ms)')
ylabel('Time (s)')
set(gca,'ydir','reverse')
set(gca,'xlim',[-15, 15])
grid on
hold on
ind_plot_2 = ind_2*1000;
plot(ind_plot_2,ind_x_2)
ind_plot_3 = ind_3*1000;
plot(ind_plot_3,ind_x_3)

if nargin > 6
    scatter(-pts(:,1)*1000,pts(:,2),10,'ko','filled');
end

xtemp=[];
ytemp=[];

while true
    % Use ginput to get x and y coordinates
    [xtemp(end+1), ytemp(end+1), button] = ginput(1); %#ok<AGROW>


    % Plot the point
    plot(xtemp(end), ytemp(end), 'ro'); % Change 'ro' to customize the appearance of the plotted points

    % Exit loop if 'q' is pressed
    if strcmp(char(button), 'q')
        break;
    end

end

slopes = zeros(length(xtemp)-1,1);
slopelocs = zeros(length(xtemp)-1,2);
for count = 1:length(slopes)
    tempdx = abs(xtemp(count+1) - xtemp(count));
    tempdy = abs(ytemp(count+1) - ytemp(count));
    slopes(count) = tempdx/tempdy;

    slopelocs(count,1) = mean([xtemp(count+1);xtemp(count)]);
    slopelocs(count,2) = mean([ytemp(count+1);ytemp(count)]);
end

end