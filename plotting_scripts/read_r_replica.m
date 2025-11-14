function [F] = read_r_replica(filename,NBAND,NFP);

ip = fopen(filename,'rt');   % 'rt' means read text
if (ip < 0)
    error('could not open file');% just abort if error
end;

%------------------------------------------------­------
%   find length of longest line
%------------------------------------------------­------

mx=0;                    % record length of longest string
cnt=0;                   % record number of strings
s = fgetl(ip);           % get a line
while (ischar(s))        % while not end of file
   cnt = cnt+1;
   %if(cnt==1);NBAND = length(str2num(s));end;       % Number freq bands
   %if(cnt==2);NFP = length(str2num(s));end;         % Number forward parameters
   if (length(s) > mx)                     % keep record of longest line
     mx = length(s);
   end;
   s = fgetl(ip);      % get next line
end;

frewind(ip);% rewind the file to the beginning

% create an empty matrix of appropriate size
tab=char(zeros(cnt,mx));% fill with ASCII zeros
% load the strings for real
cnt=0;
s = fgetl(ip);
while (ischar(s))
  cnt = cnt+1;
  tab(cnt,1:length(s)) = s;% slot into table
  s = fgetl(ip);
end;
% close the file and return
fclose(ip);

F(1).NDPF   = str2num(tab(1,:));
%F(1).minlim = str2num(tab(2,:));
%F(1).maxlim = str2num(tab(3,:));
%F(1).mtrue  = str2num(tab(4,:));
for ifreq = 2:NBAND+1

    s2n = str2num(tab(ifreq,:));
    F(ifreq-1).dat = s2n;

end
for ifreq = NBAND+2:2*NBAND+1

    s2n = str2num(tab(ifreq,:));
    F(ifreq-1-NBAND).ang = s2n;

end

%figure(1);
%for i=1:NBAND;
%   subplot(3,4,i);hold on;box on;
%   plot(F(i).ang,F(i).dat,'+k');
%   set(gca,'XLim',[0 90],'YLim',[0 40]);
%end;

return;
