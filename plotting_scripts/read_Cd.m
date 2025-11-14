function [F] = read_Cd(filename);

ip = fopen(filename,'rt');   % 'rt' means read text
if (ip < 0)
    error('could not open file');% just abort if error
end;

%------------------------------------------------­------
%   find length of longest line
%------------------------------------------------­------

NBAND = 7;
NDPF = [51, 94, 110, 118, 121, 121, 120];
mx=0;                    % record length of longest string
cnt=0;                   % record number of strings
s = fgetl(ip);           % get a line
while (ischar(s))        % while not end of file
   cnt = cnt+1;
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

jr = 1;
for ifreq = 1:NBAND
   for ir = 1:NDPF(ifreq)
      s2n = str2num(tab(jr,:));
      F(ifreq).icsave(ir,:) = s2n;
      jr = jr+1;
   end
end
for ifreq = 1:NBAND
   for ir = 1:NDPF(ifreq)
      s2n = str2num(tab(jr,:));
      F(ifreq).csave(ir,:) = s2n;
      jr = jr+1;
   end
end

return;
