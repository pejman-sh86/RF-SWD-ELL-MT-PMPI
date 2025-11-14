function [IMAP,ICOV,NVMX,ILOC,IAR,IEXCHANGE,IDBLCPL,...
         NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
         minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
         mindelay,maxdelay,minMw,minstr,mindip,minrak,...
         maxMw,maxstr,maxdip,maxrak]=cmt_read_parfile(parfile);

fid = fopen(parfile);
s = fgets(fid); %% IMAP:       Predict data for MAP file and exit.
IMAP = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IDBLCPL
IDBLCPL = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% ICOV:       Selects likelihood function (0=implicit, 1=)
ICOV = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NVMX:       Max number of nodes.  
NVMX = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IVRUP:      Vrup switch, 1 if sampling Vr on grid.
ILOC  = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IAR:        Switches on autoregressive error model
IAR = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IEXCHANGE:  Switches on exchange moves (parallel tempering)
IEXCHANGE = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NPTCHAINS1: No. chains at T = 1.
NPTCHAINS1 = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% dTlog:      PT chain spacing.
dTlog = sscanf(s, '%f'); % convert to number
s = fgets(fid); %% ICHAINTHIN: Amount of chain thinning (keep low, more PT chains are better)
ICHAINTHIN = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NKEEP:      Buffer length for saving sample.
NKEEP = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IADAPT:     Adapt the step size as function of acceptance.
IADAPT = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NBUF:
NBUF = sscanf(s, '%d');
s = fgets(fid); %% TCHCKPT:
for i=1:NVMX;
  s = fgets(fid); %% minMw:
  minMw(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxMw:
  maxMw(i) = sscanf(s, '%f');
  s = fgets(fid); %% minstr:
  minstr(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxstr:
  maxstr(i) = sscanf(s, '%f');
  s = fgets(fid); %% mindip:
  mindip(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxdip:
  maxdip(i) = sscanf(s, '%f');
  s = fgets(fid); %% minrak:
  minrak(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxrak:
  maxrak(i) = sscanf(s, '%f');
  s = fgets(fid); %% minlat:
  minlat(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxlat:
  maxlat(i) = sscanf(s, '%f');
  s = fgets(fid); %% minlon:
  minlon(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxlon:
  maxlon(i) = sscanf(s, '%f');
  s = fgets(fid); %% mindepth:
  mindepth(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxdepth:
  maxdepth(i) = sscanf(s, '%f');
  s = fgets(fid); %% mindelay:
  mindelay(i) = sscanf(s, '%f');
  s = fgets(fid); %% maxdelay:
  maxdelay(i) = sscanf(s, '%f');
end;
return;
