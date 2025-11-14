function  [out_regs] = zinput(arg1)
% ZINPUT - Graphical input from mouse with zoom
%   
%   [OUT_REGS] = ZINPUT(N) gets N points or regions from the
%   current axes and returns the X- and Y-ranges in a length Nx4
%   matrix OUT_REGS. 
%   The cursor can be positioned using a mouse. Data points are entered by
%   pressing the right mouse button, the region of the current axis
%   are selected with the middle button, left button is for
%   zooming,  single cklick zooms in, click-and drag zooms to
%   region (and doubble-click should zoom out - feature pending).
%   Any key on the keyboard except carriage return zooms out to the
%   orignal axis except carriage return, which terminates the input
%   before N points are entered.
%   
%   [OUT_REGS] = ZINPUT gathers an unlimited number of points until the
%   return key is pressed.
%   
%   See also GINPUT

% Copyright Bjorn Gustavsson 20050314

ax0 = axis;
out_regs = [];
c = computer;
if ~strcmp(c(1:2),'PC') 
  tp = get(0,'TerminalProtocol');
else
  tp = 'micro';
end

if ~strcmp(tp,'none') & ~strcmp(tp,'x') & ~strcmp(tp,'micro') & 0,
  % I dont know about this so better make short-cut and blindly try
  % what works for X in all environments - sorry about that.
else
  
  fig = gcf;
  figure(gcf);
  
  if nargin == 0
    how_many = inf;
    b = [];
  else
    how_many = arg1;
    b = [];
    if  isstr(how_many) ...
          | size(how_many,1) ~= 1 | size(how_many,2) ~= 1 ...
          | ~(fix(how_many) == how_many) ...
          | how_many < 0
      error('Requires a positive integer.');
    end
    if how_many == 0
      ptr_fig = 0;
      while(ptr_fig ~= fig)
        ptr_fig = get(0,'PointerWindow');
      end
      scrn_pt = get(0,'PointerLocation');
      loc = get(fig,'Position');
      pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
      out1 = pt(1); y = pt(2);
    elseif how_many < 0
      error('Argument must be a positive integer.');
    end
  end
  
  % Remove figure button functions
  state = uisuspend(fig);
  pointer = get(gcf,'pointer');
  set(gcf,'pointer','fullcrosshair');
  
  fig_units = get(fig,'units');
  char = 0;
  while  size(out_regs,1) < how_many
    % Use no-side effect WAITFORBUTTONPRESS
    waserr = 0;
    try
      keydown = wfbp;
    catch
      waserr = 1;
    end
    if(waserr == 1)
      if(ishandle(fig))
        set(fig,'units',fig_units);
        uirestore(state);
        error('Interrupted');
      else
        error('Interrupted by figure deletion');
      end
    end
    
    ptr_fig = get(0,'CurrentFigure');
    if(ptr_fig == fig)
      if keydown
        axis(ax0);
        char = get(fig, 'CurrentCharacter');
      else
        char = get(fig, 'CurrentCharacter');
        button = abs(get(fig, 'CurrentCharacter'));
        
        pnt = get(gcf,'currentpoint');
        xy1 = get(gca,'currentpoint');
        rbbox([pnt 0 0],pnt)
        xy2 = get(gca,'currentpoint');
        selection_type = get(gcf,'selectiontype');
        
        if all(xy1==xy2)
          ax1 = axis;
          
          dx = abs(ax1(2)-ax1(1));
          dy = abs(ax1(4)-ax1(3));
          
          xmin = max(ax1(1),xy1(1)-dx/4);
          xmax = min(ax1(2),xy1(1)+dx/4);
          ymin = max(ax1(3),xy1(1,2)-dy/4);
          ymax = min(ax1(4),xy1(1,2)+dy/4);
          zoom2ax = [xmin xmax ymin ymax];
        else
          %%% zoom to selected rectangle
          zoom2ax = [sort([xy1(1,1) xy2(1,1)]) sort([xy1(1,2) xy2(1,2)])];
        end
        
        switch selection_type
         case 'normal'
          axis(zoom2ax)
          reg = [];
         case 'alt'
          reg = axis;
          axis(ax0)
         case 'extend'
          reg = xy1(1,[1 1 2 2]);
          axis(ax0)
         otherwise %%% open (doubleclick in linux)
          axis(ax0);
          reg = [];
        end
      end
      pt = get(gca, 'CurrentPoint');
      
      if (char == 'r')
        rmi = size(out_regs,1);
        if rmi > 0
          out_regs(rmi,:) = [];
        end
        set(fig, 'CurrentCharacter','q')
        char = 'q';
        reg = [];
      end
      if(char == 13) % & how_many ~= 0)
                     % if the return key was pressed, char will == 13,
                     % and that's our signal to break out of here whether
                     % or not we have collected all the requested data
                     % points.  
                     % If this was an early breakout, don't include
                     % the <Return> key info in the return arrays.
                     % We will no longer count it if it's the last input.
        break;
      end
      
      if ~isempty(reg)
        out_regs = [out_regs;reg];
      end
    end
    %[size(out_regs,1), how_many, size(out_regs,1) < how_many]
  end
  
  uirestore(state);
  set(gcf,'pointer','arrow');
  set(fig,'units',fig_units);
  set(fig, 'CurrentCharacter','q');
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = [];

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
  h=findall(fig,'type','uimenu','accel','C');   % Disabling ^C for edit menu so the only ^C is for
  set(h,'accel','');                            % interrupting the function.
  keydown = waitforbuttonpress;
  current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
  if~isempty(current_char) & (keydown == 1)           % If the character was generated by the 
    if(current_char == 3)                       % current keypress AND is ^C, set 'waserr'to 1
      waserr = 1;                             % so that it errors out. 
    end
  end
  
  set(h,'accel','C');                                 % Set back the accelerator for edit menu.
catch
  waserr = 1;
end
drawnow;
if(waserr == 1)
  set(h,'accel','C');                                % Set back the accelerator if it errored out.
  error('Interrupted');
end

selection_type = get(gcf,'selectiontype');
if strcmp(selection_type,'open')
  axis(ax0)
end
if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

