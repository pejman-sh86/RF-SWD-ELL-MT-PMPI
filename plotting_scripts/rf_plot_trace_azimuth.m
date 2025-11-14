function []=rf_plot_trace_azimuth(waveform,nmx,col);

nrf  = size(waveform,1);
nsmp = size(waveform,2);
maxamp = max(max(waveform));
waveform = waveform/maxamp;

for irf=1:nrf;
  plot(waveform(irf,:) + (irf/2),col);
  hold on;box on;
end;
set(gca,'XLim',[0,nmx],'YLim',[0 nrf/2+maxamp])

return;
