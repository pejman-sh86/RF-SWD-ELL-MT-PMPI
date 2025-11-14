function Digitize_map(action)

%  Digitize map has been developed for a french coastal company
%  This software allows you to:
%               digitize maps and get [x-y-z] coordinates (left click)
%               load digitized points
%               enter the depth for each digitized point
%               enter the depth for each digitized line
%               zoom on the map (Tools menu)
%               delete the last digitized points (click on the wheel)
%               save each line/group of point in a file (right click)
%  If a problem occurs, don't worry! the last serie of digitized points are
%  saved in a matrice called X_digitized
%  You have to save the line digitized at the depth m before to digitize another one at the depth n
%
% Note: When you enter the depth, nothing happens when you press enter on the keyboard
% Press 2 times the key 'tab' and 1 time the key 'space' instead
%
%  Only problem: You will have to work out the transposition of the x-y coordinates:
%  The origin is at the top/left corner of the image
%  This program is a basic one, no error messages are considered
%
%   Acknowledgements:
%   This was developed based on the function "digitize2"
%   which is available from the MATLAB Central File Exchange.
%
%   Author:
%   Arnaud Gizolme (arnaud_gizolme@yahoo.fr)
%   01-October-2008


clear all
close all
screensize = get(0,'ScreenSize');
figpos = screensize*0.8 + [0.1*screensize(3),0.1*screensize(4),0,0];
hFig = figure('Name','DIGITIZE MAP: Map Digitization Tool', ...
    'NumberTitle','off','Visible','off',... 
    'HandleVisibility','on','Position',figpos, ...
    'BusyAction','Queue','Interruptible','off', ...
    'Color', 0.8*[1,1,1],'NextPlot','Add', ...
    'DoubleBuffer','On','IntegerHandle','off');
hLoad = uimenu('Label','Digitazing');

uimenu(hLoad,'Label','Image File','Callback',@LoadImageSession);
uimenu(hLoad,'Label','Digitize points','Callback',@digitize23);
uimenu(hLoad,'Label','Digitize lines','Callback',@digitize24);
uimenu(hLoad,'Label','Load points','Callback',@LoadPoints);
set(hFig,'Visible','On')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LoadImageSession(varargin)
hFig = get(0,'CurrentFigure');
[iname,ipath] = uigetfile('*.*','Select an Image-File:');
fullname = [ipath,iname];
X = imread(fullname);
axes('Position',[0,0,1,1])
ih = imagesc(X);
hold on
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LoadPoints(varargin)
hFig = get(0,'CurrentFigure');
[iname,ipath] = uigetfile('*.*','Select an Image-File:');
fullname = [ipath,iname];
hold on
PP=load (fullname);
if isnumeric(PP)
    X=PP;
else
    names=fieldnames(PP);
    X=getfield(PP,char(names(1)));
end
plot(X(:,1),X(:,2),'g.','MarkerSize',17)
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X=digitize23(varargin)

% Commence Data Acquisition from image
msgStr{1} = 'Click with LEFT mouse button to ACQUIRE';
msgStr{2} = ' ';
msgStr{3} = 'Click with RIGHT mouse button to QUIT';
titleStr = 'Ready for data acquisition';
uiwait(msgbox(msgStr,titleStr,'warn','modal'));
drawnow

numberformat = '%6.2f';
nXY = [];
while 1,
     fprintf(['\n INFO >> Click with RIGHT mouse button to QUIT \n\n']);
     n = 1;
     disp(sprintf('\n %s \n',' Index          X            Y        Z'))
     
% %%%%%%%%%%%%%% DATA ACQUISITION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     while 1
	    [x,y, buttonNumber] = ginput(1);                       
	    if buttonNumber == 1, 
          % Prompt user for depth
           prompt={'Enter the depth (Z value) at the selected point'};
           def={'1'};
           dlgTitle='DIGITIZE: user input required';
           lineNo=1;
           answer=inputdlg(prompt,dlgTitle,lineNo,def);
           if isempty(answer); continue ; end
           z= str2num(char(answer{:}));

	       H(n)=plot(x,y,'r.');
	       xpt(n) = x;
	       ypt(n) = y;
           zpt(n) = z;
           
	       disp(sprintf(' %4d         %f      %f',n, x, y,z))
	       nXY(n,:) = [n x y z];
           n = n+1;
        elseif buttonNumber == 2,
           query = questdlg('Delete last point ?', 'DIGITIZE MAP: delete point', 'YES', 'NO', 'NO');
           drawnow
	       switch upper(query),
		      case 'YES',
                 n=n-1;
			     nXY(n,:) = [];
                 xpt(n) = []; ypt(n) = []; zpt(n) = [];
                 delete(H (n))
                 disp(' Last point deleted')
		      case 'NO',
           end % switch query
           
	    else  %if buttonNumber == 3, 
	       query = questdlg('STOP digitizing and QUIT ?', 'DIGITIZE: confirmation', 'YES', 'NO', 'NO');
	       drawnow
	       switch upper(query),
		      case 'YES',
			     disp(sprintf('\n'))
			     break
		      case 'NO',
                 end % switch query
	    end
        X_digitized=[xpt' ypt' zpt'];
        save X_digitized X_digitized
     end

     if nargout  == 0,
	    % Save data to file
	    [writefname, writepname] = uiputfile('*.dat','Save data as');
        if writefname(1)==0; break;end
	    writepfname = fullfile(writepname, writefname);
	    writedata = [xpt' ypt' zpt'];
	    fid = fopen(writepfname,'w');
	    fprintf(fid,'%g %g %g\n',writedata');
	    fclose(fid);
	    disp(sprintf('\n'))
     elseif nargout == 1,
	    outputdata = [xpt' ypt' zpt'];
	    varargout{1} = outputdata;
     end
     break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X=digitize24(varargin)

% Commence Data Acquisition from image
msgStr{1} = 'Click with LEFT mouse button to ACQUIRE';
msgStr{2} = ' ';
msgStr{3} = 'Click with RIGHT mouse button to QUIT';
titleStr = 'Ready for data acquisition';
uiwait(msgbox(msgStr,titleStr,'warn','modal'));
drawnow

numberformat = '%6.2f';
nXY = [];
ng = 0;
while 1,
     fprintf(['\n INFO >> Click with RIGHT mouse button to QUIT \n\n']);
     n = 1;
     disp(sprintf('\n %s \n',' Index          X            Y        Z'))
     % Prompt user for depth
       prompt={'Enter the depth (Z value) at the selected point'};
       def={'1'};
       dlgTitle='DIGITIZE MAP: user input required';
       lineNo=1;
       answer=inputdlg(prompt,dlgTitle,lineNo,def);
       if isempty(answer); continue ; end
       z= str2num(char(answer{:}));
     
% %%%%%%%%%%%%%% DATA ACQUISITION LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     while 1
	    [x,y, buttonNumber] = ginput(1);                       
	    if buttonNumber == 1, 
           H(n)=plot(x,y,'r.');
	       xpt(n) = x;
	       ypt(n) = y;
           zpt(n) = z;
	       disp(sprintf(' %4d         %f      %f',n, x, y,z))
	       nXY(n,:) = [n x y z];
           n = n+1;
        elseif buttonNumber == 2,
           query = questdlg('Delete last point ?', 'DIGITIZE MAP: delete point', 'YES', 'NO', 'NO');
           drawnow
	       switch upper(query),
		      case 'YES',
                 n=n-1;
			     nXY(n,:) = [];
                 xpt(n) = []; ypt(n) = []; zpt(n) = [];
                 delete(H (n))
                 disp(' Last point deleted')
		      case 'NO',
           end % switch query
	    else  %if  
	       query = questdlg('STOP digitizing and QUIT ?', 'DIGITIZE MAP: confirmation', 'YES', 'NO', 'NO');
	       drawnow
	       switch upper(query),
		      case 'YES',
			     disp(sprintf('\n'))
			     break
		      case 'NO',
           end % switch query
	    end
        X_digitized=[xpt' ypt' zpt'];
        save X_digitized X_digitized
     end

     if nargout  == 0,
	    % Save data to file
	    [writefname, writepname] = uiputfile('*.dat','Save data as');
        if writefname(1)==0; break;end
	    writepfname = fullfile(writepname, writefname);
	    writedata = [xpt' ypt' zpt'];
	    fid = fopen(writepfname,'w');
	    fprintf(fid,'%g %g %g\n',writedata');
	    fclose(fid);
	    disp(sprintf('\n'))
     elseif nargout == 1,
	    outputdata = [xpt' ypt' zpt'];
	    varargout{1} = outputdata;
     end
     break
end   