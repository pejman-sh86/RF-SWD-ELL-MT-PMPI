function [R,ang,freq] = read_oases_r(filename);

ip = fopen(filename,'rt');   % 'rt' means read text
if (ip < 0)
    error('could not open file');% just abort if error
end;


%------------------------------------------------­------
%   find length of longest line
%------------------------------------------------­------

max=0;          % record length of longest string
cnt=0;          % record number of strings
s = fgetl(ip);              % get a line
while (ischar(s))     % while not end of file
   cnt = cnt+1;
   if (length(s) > max)    % keep record of longest
     max = length(s);
   end;
   s = fgetl(ip);      % get next line
end;

frewind(ip);% rewind the file to the beginning

% create an empty matrix of appropriate size
tab=char(zeros(cnt,max));% fill with ASCII zeros
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

s2n = str2num(tab(1,:));
fbeg = s2n(1);
fend = s2n(2);
nfreq = s2n(3);
freq = [fbeg:(fend-fbeg)/(nfreq-1):fend];

s1 = strrep(tab(2,:),'# Frequency, # of angles','');

s2n = str2num(s1);
nang = s2n(2);
clear s2n;

irow = 3;
for ifreq = 1:nfreq

    s2n = str2num(tab(irow:irow+nang-1,:));
    ang = s2n(:,1);
    abs_R(:,ifreq) = s2n(:,2);
    fi_R(:,ifreq) = s2n(:,3);
    irow = irow + nang + 1;

end

R = abs_R .* exp(i*fi_R);
save bla R abs_R ang;

return;
