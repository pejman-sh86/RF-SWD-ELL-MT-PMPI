
filebase = 'HON';
include_rv = 0;
include_swd = 0;
include_ell = 0;
include_mt = 1;
%%
repfile     = strcat(filebase,'_replica.dat');
dobsfile    = strcat(filebase,'_dobs.dat');
%%
dpred = [];
dobs = [];
%%

if include_rv == -1
   predrf = importdata(strcat(filebase, '_mappred.dat')); % 1xNTIME array
   %NTIME = length(predrf);
   datarf = importdata(strcat(filebase, '_RF.txt')); % 2xNTIME array
   obsrf = datarf(1,:); % 1xNTIME array

   dpred = [dpred predrf];
   dobs = [dobs obsrf];
end

if include_swd == 1
    predswd = importdata(strcat(filebase, '_mappredSWD.dat')); % 1xNDAT_SWD array
    dataswd = importdata(strcat(filebase, '_SWD.dat')); % NDAT_SWDx2 array
    obsswd = dataswd(:,2)'; % 1xNDAT_SWD array

    dpred = [dpred predswd];
    dobs = [dobs obsswd];
end

if include_ell == 1
    predell = importdata(strcat(filebase, '_mappredELL.dat')); % 1xNDAT_ELL array
    dataell = importdata(strcat(filebase, '_ELL.dat')); % NDAT_ELLx2 array
    obsell = dataell(:,2)'; % 1xNDAT_ELL array

    dpred = [dpred predell];
    dobs = [dobs obsell];
end

if include_mt == 1
    predmt = importdata(strcat(filebase, '_mappredMT.dat')); % 1x2NDAT_MT array
    datamt = importdata(strcat(filebase, '_MT.dat')); % NDAT_MTx3 array
    obsmt = reshape(datamt(:,2:3),size(predmt)); % 1x2NDAT_MT array

    dpred = [dpred predmt];
    dobs = [dobs obsmt];
end

save(repfile, 'dpred', '-ascii')
save(dobsfile, 'dobs', '-ascii')
