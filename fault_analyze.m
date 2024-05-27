function fault_analyze(index,Df,H076,H090,segy_filename,Data,trace_vec,time_vec)

%   index=8;

% get fault center x-position
fault_x_pos = mean([Df.pt1_trace(index); Df.pt2_trace(index)]);
mat_min = round(fault_x_pos - 500);
mat_max = round(fault_x_pos + 500);
mat_min_2 = round(fault_x_pos - 50);
mat_max_2 = round(fault_x_pos + 50);

[Data2,SegyTraceHeader2,SegyHeade2]=ReadSegy(segy_filename,'minmax','cdp',mat_min_2,mat_max_2); %#ok<ASGLU>
trace_vec_2 = mat_min_2:1:mat_max_2;
time_vec_2 = SegyHeade2.time;

% make position mesh
[XX2,YY2] = meshgrid(trace_vec_2,time_vec_2);

Fseis = scatteredInterpolant(XX2(:),YY2(:),Data2(:));


fault_xvec_1 = linspace(Df.pt1_trace(index)-10,Df.pt2_trace(index)-15,1000);
fault_xvec_2 = linspace(Df.pt1_trace(index)+10,Df.pt2_trace(index)+15,1000);
fault_yvec = linspace(Df.pt1_time(index),Df.pt2_time(index),1000);
dy = fault_yvec(2) - fault_yvec(1);

fault_amp_1 = Fseis(fault_xvec_1,fault_yvec);
fault_amp_2 = Fseis(fault_xvec_2,fault_yvec);

ind_1 = zeros(750,1);
ind_x_1 = zeros(750,1);

ind_2 = zeros(850,1);
ind_x_2 = zeros(850,1);

ind_3 = zeros(950,1);
ind_x_3 = zeros(950,1);

for count = 1:750
    [C,lags] = xcorr(fault_amp_1(count:count+250),fault_amp_2(count:count+250));
    ind_1(count) = lags(find(C == max(C))) * dy; %#ok<FNDSB>
    ind_x_1(count) = fault_yvec(count);
end

for count = 1:850
    [C,lags] = xcorr(fault_amp_1(count:count+150),fault_amp_2(count:count+150));
    ind_2(count) = lags(find(C == max(C))) * dy; %#ok<FNDSB>
    ind_x_2(count) = fault_yvec(count);
end

for count = 1:950
    [C,lags] = xcorr(fault_amp_1(count:count+50),fault_amp_2(count:count+50));
    ind_3(count) = lags(find(C == max(C))) * dy; %#ok<FNDSB>
    ind_x_3(count) = fault_yvec(count);
end

[xtemp,ytemp,slopes,slopelocs] = pick_slope_breaks_rp(ind_1,ind_x_1);
% close all
figure;
ax(1) = subplot(1,5,1:3);
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

ax(2) = subplot(1,5,4);
plot(fault_amp_1,fault_yvec,'k')
hold on
plot(fault_amp_2,fault_yvec,'b')
set(gca,'ydir','reverse')
xlabel('Seismic Amplitude')
ylabel('Time (s)')
grid on
set(gca,'xlim',[-5, 5])

ax(3) = subplot(1,5,5);
ind_plot_1 = ind_1*1000;
ind_plot_2 = ind_2*1000;
ind_plot_3 = ind_3*1000;
plot(ind_plot_1,ind_x_1)
xlabel('Lag (ms)')
ylabel('Time (s)')
set(gca,'ydir','reverse')
set(gca,'xlim',[-10, 10])
grid on
hold on
plot(ind_plot_2,ind_x_2)
plot(ind_plot_3,ind_x_3)
plot(xtemp,ytemp,'r+:','linewidth',1.5);
for count = 1:length(slopes)
    text(slopelocs(count,1),slopelocs(count,2),sprintf('%.2f ms/s',slopes(count)));
end


linkaxes(ax,'y')
set(gca,'ylim',[min(fault_yvec)*0.98,max(fault_yvec)*1.02])

end