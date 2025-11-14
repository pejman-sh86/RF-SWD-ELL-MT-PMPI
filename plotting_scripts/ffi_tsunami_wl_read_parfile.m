function [IMAP,ICOV,NVMX,NPV,NMISC,IVRUP,IAR,IEXCHANGE,...
         NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,...
         NBUF,MAXDISP,MINDISP]=ffi_tsunami_wl_read_parfile(parfile);

fid = fopen(parfile);

s = fgets(fid); %% IMAP:       Predict data for MAP file and exit.
IMAP = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IWLRAN:       
IWLRAN = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% ISNGLCORE:       
ISNGLCORE = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% ICOV:       Selects likelihood function (0=implicit, 1=)
ICOV = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NVMX:       Max number of nodes.  
NVMX = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NPV:        Number of parameters per node.
NPV = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NMISC:      Number of misc parameters (hyp loc and rup vel).
NMISC = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IVRUP:      Vrup switch, 1 if sampling Vr on grid.
IVRUP = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IAR:        Switches on autoregressive error model
IAR = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IEXCHANGE:  Switches on exchange moves (parallel tempering)
IEXCHANGE = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NPTCHAINS1: No. chains at T = 1.
NPTCHAINS1 = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% dTlog:      PT chain spacing.
dTlog = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% ICHAINTHIN: Amount of chain thinning (keep low, more PT chains are better)
ICHAINTHIN = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NKEEP:      Buffer length for saving sample.
NKEEP = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% IADAPT:     Adapt the step size as function of acceptance.
IADAPT = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NBUF:       Buffer length for adapting step size
NBUF = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% TCHECKPT:       Buffer length for adapting step size
TCHCKPT = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% NRESEED: No. of re-seeding sweeps (NRESEED*NTHREAD times; 0 turns off re-seeding)
NRESEED = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% MAXDISP:    Max sea-surface displacement
MAXDISP = sscanf(s, '%d'); % convert to number
s = fgets(fid); %% MINDISP:    Min sea-surface displacement
MINDISP = sscanf(s, '%d'); % convert to number

return;
