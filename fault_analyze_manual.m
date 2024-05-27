function [pts,slope_out] = fault_analyze_manual(Line_num,index,Df,H076,H090,segy_filename,Data,trace_vec,time_vec)

%   index=8;

% get fault center x-position
fault_x_pos = mean([Df.pt1_trace(index); Df.pt2_trace(index)]);
mat_min = round(fault_x_pos - 500);
mat_max = round(fault_x_pos + 500);
mat_min_2 = round(fault_x_pos - 50);
mat_max_2 = round(fault_x_pos + 50);

% read in part of the segy surrounding the fault
[Data2,SegyTraceHeader2,SegyHeade2]=ReadSegy(segy_filename,'minmax','cdp',mat_min_2,mat_max_2); %#ok<ASGLU>
trace_vec_2 = mat_min_2:1:mat_max_2;
time_vec_2 = SegyHeade2.time;

% make position mesh and interpolant
[XX2,YY2] = meshgrid(trace_vec_2,time_vec_2);
Fseis = scatteredInterpolant(XX2(:),YY2(:),Data2(:));

%--------old version------------------------------
% fault_xvec = linspace(Df.pt1_trace(index),Df.pt2_trace(index)-15,1000);
% fault_xvec_1 = linspace(Df.pt1_trace(index)-10,Df.pt2_trace(index)-15,1000);
% fault_xvec_2 = linspace(Df.pt1_trace(index)+10,Df.pt2_trace(index)+15,1000);
% fault_yvec = linspace(Df.pt1_time(index),Df.pt2_time(index),1000);
% dy = fault_yvec(2) - fault_yvec(1);
%--------------------------------------------------

%--------new test version------------------------
dy = 0.00025; % 0.25 ms for this prototype
fault_yvec = Df.pt1_time(index):dy:Df.pt2_time(index);
fault_xvec = interp1([Df.pt1_time(index),Df.pt2_time(index)],[Df.pt1_trace(index),Df.pt2_trace(index)],fault_yvec,'linear');
fault_xvec_1 = fault_xvec - 15;
fault_xvec_2 = fault_xvec + 15;
%--------------------------------------------------

% interpolate amplitudes on either side of the fault
fault_amp_1 = Fseis(fault_xvec_1,fault_yvec);
fault_amp_2 = Fseis(fault_xvec_2,fault_yvec);

% calculate lags
fault_length = length(fault_yvec);
window_1_length = 0.02; % s
window_2_length = 0.03; % s
window_3_length = 0.04; % s
window_1 = window_1_length/dy; % samples
window_2 = window_2_length/dy; % samples
window_3 = window_3_length/dy; % samples

% output lengths
lag_length_1 = fault_length - window_1;
lag_length_2 = fault_length - window_2;
lag_length_3 = fault_length - window_3;

ind_1 = zeros(lag_length_1,1);
ind_x_1 = zeros(lag_length_1,1);
ind_2 = zeros(lag_length_2,1);
ind_x_2 = zeros(lag_length_2,1);
ind_3 = zeros(lag_length_3,1);
ind_x_3 = zeros(lag_length_3,1);

for count = 1:lag_length_1
    [C,lags] = xcorr(fault_amp_1(count:count+window_1),fault_amp_2(count:count+window_1));
    ind_1(count) = lags(find(C == max(C))) * dy; %#ok<FNDSB>
    ind_x_1(count) = fault_yvec(count);
end

for count = 1:lag_length_2
    [C,lags] = xcorr(fault_amp_1(count:count+window_2),fault_amp_2(count:count+window_2));
    ind_2(count) = lags(find(C == max(C))) * dy; %#ok<FNDSB>
    ind_x_2(count) = fault_yvec(count);
end

for count = 1:lag_length_3
    [C,lags] = xcorr(fault_amp_1(count:count+window_3),fault_amp_2(count:count+window_3));
    ind_3(count) = lags(find(C == max(C))) * dy; %#ok<FNDSB>
    ind_x_3(count) = fault_yvec(count);
end


[pts] = pick_horizons_manual(index,Df,trace_vec,time_vec,Data,H076,H090,mat_min,mat_max,fault_xvec_1,fault_xvec_2,fault_yvec);
saveas(gcf,sprintf('Horizon_ctrl_L%d_F%d.png',Line_num,index));

[xtemp,ytemp,slopes,slopelocs] = pick_slope_breaks_rp(ind_1,ind_x_1,ind_2,ind_x_2,ind_3,ind_x_3,pts);
saveas(gcf,sprintf('Slope_breaks_L%d_F%d.png',Line_num,index));

% close all
figure('units','normalized','outerposition',[0 0 1 1]);
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
ind_plot_1 = ind_1*1000; % to ms
ind_plot_2 = ind_2*1000; % to ms
ind_plot_3 = ind_3*1000; % to ms
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

slope_out = [];

for count = 1:length(slopes)
    text(slopelocs(count,1),slopelocs(count,2),sprintf('%.2f ms/s',slopes(count)));
    slope_out_temp_y = linspace(ytemp(count),ytemp(count+1),10)';
    slope_out_temp_x = interp1(fault_yvec,fault_xvec,slope_out_temp_y);
    slope_out_temp_z = ones(size(slope_out_temp_y)) * slopes(count);
    slope_out = [slope_out; slope_out_temp_x, slope_out_temp_y, slope_out_temp_z];
end

scatter(-pts(:,1)*1000,pts(:,2),10,'ko','filled')

linkaxes(ax,'y')
set(gca,'ylim',[min(fault_yvec)*0.98,max(fault_yvec)*1.02])
saveas(gcf,sprintf('L%d_F%d_output.png',Line_num,index));

end