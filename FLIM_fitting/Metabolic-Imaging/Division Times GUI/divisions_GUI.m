%GUI for automated detection and correction of changepoints in embryo development time series
% Versions:
% 2018-01-02: Going backward in this GUI currently clears all the
% subsequent frames. I corrected that.

function divisions_GUI
clear all; close all; clc;


hgui = figure('position',[300 100 800 600]);
set(hgui,'NumberTitle','off');
set(hgui, 'Name','Embryo GUI');
% create structure of handles
myhandles = guihandles(hgui); 
% initialize structure variables
myhandles.s1a1dir = '';         %s1a1/ directory path
myhandles.filedir = {};         %cell containing full path to directories of images
myhandles.nfiles = struct([]);  %structure to store filenames of image tifs in directory
myhandles.ind = 0;              %index of current subdirectory
myhandles.ndir = 0;             %total number of directories loaded by load image button 
myhandles.fcell = {};           %feature cell 
myhandles.batch = 'MLsets_';           %name of prefix
myhandles.chngpnts = {};        %cell array for storing changepoints {name, [ind1,ind2,indb]} where
                                %name is in format MLsets#_Pos#_Cell#
                                %ind1 is the changepoint for 1-cell to
                                %2-cell
                                %ind2 is the changepoint for 2-cell to
                                %multi-cell
                                %indb is the changepoint for morula to
                                %blastocyst
myhandles.chngpnts_accum = {};  %cell array to accumulate chngpnts, to be saved to file
myhandles.toggle=0;             %toggle to modify ind1, ind2 or indb via the GUI
myhandles.PanelImHands = [];
setappdata(0,'hgui',gcf);
% Save the structure
guidata(hgui,myhandles)


% CREATE BUTTONS:
% Get the initial directory of images to loop over
h_load_image_file = uicontrol('Parent',hgui,'style','pushbutton','units','pixels','position',[650 380 100 50],'String','Select Sample Dir','Callback',@load_image_file_callback);

% Pos, mask, and embryo number display
h_EmbNumLab_txt = uicontrol('Parent',hgui,'Style','text','Position',[625 270 50 20],'String','Embryo:');
h_EmbNum_txt = uicontrol('Parent',hgui,'style','edit','units','pixels','position',[680 273 30 20],'String','','Callback',@prefix_EmbNum_txt_callback);
h_EmbSlash_txt = uicontrol('Parent',hgui,'Style','text','Position',[710 270 15 20],'String','\');
h_Embtot_txt = uicontrol('Parent',hgui,'style','text','units','pixels','position',[725 270 30 20],'String','','Callback',@prefix_Embtot_txt_callback);

h_PosLab_txt = uicontrol('Parent',hgui,'Style','text','Position',[625 330 50 20],'String','Pos:');
h_PosNum_txt = uicontrol('Parent',hgui,'style','text','units','pixels','position',[680 330 50 20],'String','');
h_MaskLab_txt = uicontrol('Parent',hgui,'Style','text','Position',[625 300 50 20],'String','Mask:');
h_MaskNum_txt = uicontrol('Parent',hgui,'style','text','units','pixels','position',[680 300 50 20],'String','');


% Indices and radio buttons.
h_toggle_b.c(1) = uicontrol('style','radiobutton','units','pixels','position',[645,240,40,20],'Tag','1','Callback',@toggleb_callback);
h_toggle_b.c(2) = uicontrol('style','radiobutton','units','pixels','position',[685,240,40,20],'Tag','2','Callback',@toggleb_callback);
h_toggle_b.c(3) = uicontrol('style','radiobutton','units','pixels','position',[725,240,40,20],'Tag','3','Callback',@toggleb_callback);
h_ind1_txt = uicontrol('Parent',hgui,'Style','text','Position',[635 220 50 20],'String','ind1:');
h_ind1_edit = uicontrol('Parent',hgui,'style','edit','units','pixels','position',[645 205 30 20],'String','','Callback',@prefix_ind1_callback);
h_ind2_txt = uicontrol('Parent',hgui,'Style','text','Position',[675 220 50 20],'String','ind2:');
h_ind2_edit = uicontrol('Parent',hgui,'style','edit','units','pixels','position',[685 205 30 20],'String','','Callback',@prefix_ind2_callback);
h_indb_txt = uicontrol('Parent',hgui,'Style','text','Position',[715 220 50 20],'String','indb:');
h_indb_edit = uicontrol('Parent',hgui,'style','edit','units','pixels','position',[725 205 30 20],'String','','Callback',@prefix_indb_callback);
h_previous_file = uicontrol('Parent',hgui,'style','pushbutton','units','pixels','position',[650 145 100 50],'String','Previous','Callback',@load_previous_file_callback);
h_next_file = uicontrol('Parent',hgui,'style','pushbutton','units','pixels','position',[650 80 100 50],'String','Next','Callback',@load_next_file_callback);



function load_image_file_callback(hgui,event)
   %Select and Change to s1_a1 directory
   s1a1dir = uigetdir();
   cd(s1a1dir)
   % Get the data storage structure using guidata in the local function
   myhandles = guidata(gcbo);
   myhandles.s1a1dir = s1a1dir;
   assignin('base', 's1a1dir', s1a1dir);
   %look for features file and load if it exists
   featuresfile = dir('*features.mat');
   if ~isempty(featuresfile)
        fcell={};
        load(featuresfile(1).name);        
        % Update the value of fcell
        myhandles.fcell = fcell;
        % Assign to workspace variable 'fcell' 
        assignin('base', 'fcell', fcell);
        % Save the change made to the structure
        guidata(gcbo,myhandles) 
   else
       fprintf('Error: no features file found\n');
   end
   cd('sorted_sdts')
   batch = myhandles.batch;
   names = [batch '*'];
   d = dir(names);
   nd = size(d,1); %number of position folders
   ndir = 0; %initialize number of mask directories
   filedir = cell(1,ndir);  %declare master filedir cell
   %loop over all position and mask subfolders to build master filedir
   for i=1:nd
       dirname = fullfile(s1a1dir, 'sorted_sdts', d(i).name);
       cd(dirname);
       d1=dir();
       d1 = d1([d1.isdir]);  %get list of folder names
       d1 = d1(3:end);      %get rid of . and ..
       nd1 = size(d1,1);
       ndir=ndir+nd1; 
       for j=1:nd1
            subnames = fullfile(s1a1dir, 'sorted_sdts', d(i).name,d1(j).name);
            filedir = [filedir {subnames}];
       end
   end
   myhandles.filedir = filedir;
   assignin('base', 'filedir', filedir);
   myhandles.ndir = ndir;
   assignin('base', 'ndir', ndir);
   cd(s1a1dir)
   %look for chngpnts file and load if it exists
   chngpntsfile = dir('*changepoints.mat');  
   TotNumEmbs = size(filedir,2);
   if ~isempty(chngpntsfile)
        load(chngpntsfile(1).name);        
        % Update the value of chngpnts_accum
        myhandles.chngpnts_accum = chngpnts_accum;
        % Update the valur of chngpnts_accum 
        assignin('base', 'chngpnts_accum', chngpnts_accum); 
        %Update the value of chngpnts
        chngpnts{1} = chngpnts_accum{end,1};
        chngpnts{2}(1) = chngpnts_accum{end,2}(1);
        chngpnts{2}(2) = chngpnts_accum{end,2}(2);  
        chngpnts{2}(3) = chngpnts_accum{end,2}(3);
        myhandles.chngpnts = chngpnts;
        %Broadcast to workspace
        assignin('base', 'chngpnts', chngpnts);
        %update index fields
        set(h_ind1_edit,'String',chngpnts{2}(1));
        set(h_ind2_edit,'String',chngpnts{2}(2));
        set(h_indb_edit,'String',chngpnts{2}(3));
        %Update the value of ind
        ind = find(strncmp(chngpnts_accum{end,1},chngpnts_accum(:,1),numel(chngpnts_accum{end,1})));
        myhandles.ind = ind;
        assignin('base','ind',ind);
        guidata(gcbo,myhandles) 
        draw_panel(hgui)   %draw the main panel of GUI with embyro time points
        %Update the previous value of chngpnts_accum
%         chngpnts_accum{end,1} = [];
%         chngpnts_accum{end,2} = [];
%         D1 = chngpnts_accum'; 
%         chngpnts_accum = reshape(D1(~cellfun(@isempty,D1)),2,[])'; %ugly; removes last row from chngpnts_accum
        myhandles.chngpnts_accum = chngpnts_accum;
        guidata(gcbo,myhandles)
   else
        ind = 1;
        myhandles.ind = ind;
        assignin('base', 'ind', ind);
        %analysis for first position, mask
        nfiles = dir(filedir{ind});          %only ind = 1 should be called initially 
        nfiles = nfiles(~[nfiles.isdir]);    %get only the image .tif files
        nsize = size(nfiles,1);
        for i=1:nsize
            nfiles(i).name = fullfile(filedir{ind}, nfiles(i).name);
        end;
        myhandles.toggle=0;
        assignin('base', 'toggle', myhandles.toggle);
        myhandles.nfiles=nfiles;
        assignin('base', 'nfiles', nfiles);
        name = nfiles(1).name;  %only the first filename should be called here
        slash = strfind(name,filesep);
        dot = strfind(name,'.');
        name = name((slash(end)+1):(dot(1)-1));  %extract the filename
        under = strfind(name,'_');
        name = name(1:under(2)-1);               %remove the frame number part from name
        filename = [myhandles.batch name];       %filename format MLsets_Pos#_cell#
        chngpnts = changepoints(myhandles.fcell,filename); %calculate the changepoints for the time series of filename
        %Update the value of chngpnts
        myhandles.chngpnts = chngpnts;
        %Broadcast to workspace
        assignin('base', 'chngpnts', chngpnts);
        %update index fields
        set(h_ind1_edit,'String',chngpnts{2}(1));
        set(h_ind2_edit,'String',chngpnts{2}(2));
        set(h_indb_edit,'String',chngpnts{2}(3));
        myhandles.chngpnts_accum = chngpnts;
        guidata(gcbo,myhandles) 
        draw_panel(hgui)   %draw the main panel of GUI with embyro time points 
   end
   
   set(h_EmbNum_txt,'String',ind);
   set(h_Embtot_txt,'String',TotNumEmbs);
   
end

%pushbutton_callback is executed when GUI buttons defined in draw_panel function are pressed
function pushbutton_callback(hgui,event,x)
   myhandles = guidata(gcbo);
   ind = myhandles.ind;
   nfiles=myhandles.nfiles;
   chngpnts=myhandles.chngpnts;
   switch myhandles.toggle
       case 1
           chngpnts{2}(1) = x;
           set(h_ind1_edit,'String',chngpnts{2}(1));
       case 2
           chngpnts{2}(2) = x;
           set(h_ind2_edit,'String',chngpnts{2}(2));
       case 3
           chngpnts{2}(3) = x;
           set(h_indb_edit,'String',chngpnts{2}(3));
       otherwise
           set(h_ind1_edit,'String',chngpnts{2}(1));
           set(h_ind2_edit,'String',chngpnts{2}(2));
           set(h_indb_edit,'String',chngpnts{2}(3));
   end
   myhandles.chngpnts = chngpnts; 
   myhandles.chngpnts_accum(ind,:) = chngpnts;
   % Broadcast to workspace 
   assignin('base', 'chngpnts', chngpnts);
   guidata(gcbo,myhandles) 
   draw_panel(hgui);
end


function toggleb_callback(hgui,event)
    % Get the structure using guidata in the local function
    myhandles = guidata(gcbo);
    if myhandles.toggle ~= 0 
        set(h_toggle_b.c(myhandles.toggle),'Value',0); %unclick previous radiobutton if clicked, ensures only one radio is checked
    end
    vals = get(h_toggle_b.c,'Value');
    button_state = find([vals{:}] == 1);
    if ~isempty(button_state)
        myhandles.toggle = button_state;
    else 
        myhandles.toggle = 0;
    end    
    % Broadcast to workspace 
    assignin('base', 'toggle', button_state);
    guidata(gcbo,myhandles)
end

function prefix_EmbNum_txt_callback(hgui,event)
    myhandles = guidata(gcbo);
    chngpnts_accum = myhandles.chngpnts_accum;
    ind = str2num(get(h_EmbNum_txt,'String'));
    % If we're already at the first frame, just don't do anything else.
    if ind>0
        if myhandles.toggle ~= 0
            set(h_toggle_b.c(myhandles.toggle),'Value',0); %unclick previous radiobutton
        end
        myhandles.toggle = 0;
        %Broadcast to workspace
        assignin('base', 'toggle', 0);
        chngpnts = myhandles.chngpnts_accum(ind,:);
        myhandles.chngpnts = chngpnts;
        
        %restore previous index values
        set(h_ind1_edit,'String',chngpnts{2}(1));
        set(h_ind2_edit,'String',chngpnts{2}(2));
        set(h_indb_edit,'String',chngpnts{2}(3));
        
        % Broadcast to workspace
        assignin('base', 'ind', ind);
        assignin('base', 'chngpnts_accum', chngpnts_accum);
        myhandles.ind = ind;
        guidata(gcbo,myhandles)
        set(h_EmbNum_txt,'String',ind);
        %Draw Panel
        draw_panel(hgui)
    else
        set(h_EmbNum_txt,'String',num2str(myhandles.ind));
    end
        
end

function prefix_ind1_callback(hgui,event)
    % Get the structure using guidata in the local function
    myhandles = guidata(gcbo);
    ind = myhandles.ind;
    x = get(h_ind1_edit,'String'); %ind1 being Tag of ind1 edit box
    if isempty(x)
        fprintf('Error: Enter Text first\n');
    else
        % Update the value of ind1
        chngpnts = myhandles.chngpnts;
        chngpnts{2}(1) = str2num(x);
        myhandles.chngpnts = chngpnts;
        myhandles.chngpnts_accum(ind,:) = chngpnts;
    	% Broadcast to workspace 
        assignin('base', 'ind1', chngpnts{2}(1));
        assignin('base', 'chngpnts', chngpnts);
        guidata(gcbo,myhandles)
        %update panel
        draw_panel(hgui)
    end
end


function prefix_ind2_callback(hgui,event)
    % Get the structure using guidata in the local function
    myhandles = guidata(gcbo);
    ind = myhandles.ind;
    x = get(h_ind2_edit,'String'); %ind2 being Tag of ind2 edit box
    if isempty(x)
        fprintf('Error: Enter Text first\n');
    else
        % Update the value of ind2
        chngpnts = myhandles.chngpnts;
        chngpnts{2}(2) = str2num(x);
        myhandles.chngpnts = chngpnts;
        myhandles.chngpnts_accum(ind,:) = chngpnts;
    	% Broadcast to workspace 
        assignin('base', 'ind2', chngpnts{2}(2));
        assignin('base', 'chngpnts', chngpnts);
        guidata(gcbo,myhandles)
        %update panel
        draw_panel(hgui)
    end
end


function prefix_indb_callback(hgui,event)
    %Get the structure using guidata in the local function
    myhandles = guidata(gcbo);
    ind = myhandles.ind;
    x = get(h_indb_edit,'String'); %indb being Tag of indb edit box
    if isempty(x)
        fprintf('Error: Enter Text first\n');
    else
        %Update the value of indb
        chngpnts = myhandles.chngpnts;
        chngpnts{2}(3) = str2num(x);
        myhandles.chngpnts = chngpnts;
        myhandles.chngpnts_accum(ind,:) = chngpnts;
    	%Broadcast to workspace 
        assignin('base', 'indb', chngpnts{2}(3));
        assignin('base', 'chngpnts', chngpnts);
        guidata(gcbo,myhandles) 
        %Update panel
        draw_panel(hgui)
    end
end

function load_previous_file_callback(hgui,event)
    myhandles = guidata(gcbo);
    chngpnts_accum = myhandles.chngpnts_accum;
    ind = myhandles.ind - 1;
    % If we're already at the first frame, just don't do anything else.
    if ind>0
        
        if myhandles.toggle ~= 0
            set(h_toggle_b.c(myhandles.toggle),'Value',0); %unclick previous radiobutton
        end
        myhandles.toggle = 0;
        %Broadcast to workspace
        assignin('base', 'toggle', 0);
        chngpnts = myhandles.chngpnts_accum(ind,:);
        myhandles.chngpnts = chngpnts;
        
        %restore previous index values
        set(h_ind1_edit,'String',chngpnts{2}(1));
        set(h_ind2_edit,'String',chngpnts{2}(2));
        set(h_indb_edit,'String',chngpnts{2}(3));
        
%         
%         if ~isempty(chngpnts_accum)
%             chngpnts{1} = chngpnts_accum{end,1};
%             chngpnts{2}(1) = chngpnts_accum{end,2}(1);
%             chngpnts{2}(2) = chngpnts_accum{end,2}(2);
%             chngpnts{2}(3) = chngpnts_accum{end,2}(3);
%             myhandles.chngpnts = chngpnts;
%             
%             %         %Update the previous value of chngpnts_accum
%             %         chngpnts_accum{end,1} = [];
%             %         chngpnts_accum{end,2} = [];
%             %         D1 = chngpnts_accum';
%             %         chngpnts_accum = reshape(D1(~cellfun(@isempty,D1)),2,[])'; %ugly; removes last row from chngpnts_accum
%             myhandles.chngpnts_accum = chngpnts_accum;
%         else
%             ind = myhandles.ind;
%         end
        % Broadcast to workspace
        assignin('base', 'ind', ind);
        assignin('base', 'chngpnts_accum', chngpnts_accum);
        myhandles.ind = ind;
        guidata(gcbo,myhandles)
        set(h_EmbNum_txt,'String',ind);
        %Draw Panel
        draw_panel(hgui)
        
    end
end
    

function load_next_file_callback(hgui,event)
    myhandles = guidata(gcbo);
    chngpnts_accum = myhandles.chngpnts_accum;
    ind = myhandles.ind + 1;

    if ind == myhandles.ndir + 1  %case of last image file
        % Change: do nothing
        %         % Update the value of chngpnts
        %         chngpnts = myhandles.chngpnts;
        %         % Update the value of chngpnts_accum
        %         chngpnts_accum = myhandles.chngpnts_accum;
        %         chngpnts_accum = [chngpnts_accum;chngpnts];
        %         myhandles.chngpnts_accum = chngpnts_accum;
        %         %update ind
        %         myhandles.ind = ind;
        %         % Broadcast to workspace
        %         assignin('base', 'chngpnts', chngpnts);
        %         assignin('base', 'chngpnts_accum', chngpnts_accum);
        %         assignin('base', 'ind', ind);
        %         %save chngpnts_accum to s1a1/
        %         s1a1dir = myhandles.s1a1dir;
        %         save(fullfile(s1a1dir,'embryo_data_changepoints.mat'),'chngpnts_accum');
        %         guidata(gcbo,myhandles)
        
    elseif ind <= myhandles.ndir
        
        
        % If change points were previously specified and we are reviewing
        % them, load them from chngpnts_accum
        if ind<=size(chngpnts_accum,1)
            chngpnts = chngpnts_accum(ind,:);
        else % otherwise, calculate the changepoints using the machine learning magic
            % Update the value of chngpnts
            chngpnts = myhandles.chngpnts;
            % Update the value of chngpnts_accum
            chngpnts_accum = myhandles.chngpnts_accum;
            chngpnts_accum(ind,:) = chngpnts;
            myhandles.chngpnts_accum = chngpnts_accum;
            % Broadcast to workspace
            assignin('base', 'chngpnts', chngpnts);
            assignin('base', 'chngpnts_accum', chngpnts_accum);
            if myhandles.toggle ~= 0
                set(h_toggle_b.c(myhandles.toggle),'Value',0); %unclick previous radiobutton
            end
            myhandles.toggle = 0;
            
            filedir = myhandles.filedir;
            nfiles = dir(filedir{ind});
            nfiles = nfiles(~[nfiles.isdir]);  %get only the .tif files
            nsize = size(nfiles,1);
            %Broadcast to workspace
            myhandles.toggle=0;
            assignin('base', 'toggle', 0);
            myhandles.nfiles=nfiles;
            assignin('base', 'nfiles', nfiles);
            for i=1:nsize
                nfiles(i).name = fullfile(filedir{ind}, nfiles(i).name);
            end
            k = size(nfiles);
            n_row = ceil(sqrt(k(1)));
            n_col = floor(k(1)/n_row);
            name = nfiles(1).name;  %only the first filename is needed to lookup in fcell
            slash = strfind(name,filesep);
            dot = strfind(name,'.');
            name = name((slash(end)+1):(dot(1)-1));  %extract the filename
            under = strfind(name,'_');
            name = name(1:under(2)-1);               %remove the frame number part from name
            filename = [myhandles.batch name];       %format Batch#_Pos#_cell#
            chngpnts = changepoints(myhandles.fcell,filename); %calculate the changepoints for the time series of filename
        end
        %Update the value of chngpnts and chngpnts_accum
        myhandles.chngpnts = chngpnts;
        myhandles.chngpnts_accum(ind,:) = chngpnts;
        %update ind
        myhandles.ind = ind;
        %Broadcast to workspace
        assignin('base', 'chngpnts', chngpnts);
        assignin('base', 'ind', ind);
        %update index fields
        set(h_ind1_edit,'String',chngpnts{2}(1));
        set(h_ind2_edit,'String',chngpnts{2}(2));
        set(h_indb_edit,'String',chngpnts{2}(3));
        guidata(gcbo,myhandles)
        draw_panel(hgui)
        set(h_EmbNum_txt,'String',ind);
    end
    
end

% Save the current chngpnts_accum handle variable file named
% "embryo_data_changepoints.mat" in s1a1/
h_save_file = uicontrol('Parent',hgui,'style','pushbutton','units','pixels','position',[650 20 100 50],'String','Save','Callback',@save_file_callback);

function save_file_callback(hgui,event) 
    myhandles = guidata(gcbo);
    chngpnts_accum = myhandles.chngpnts_accum;
    s1a1dir = myhandles.s1a1dir;
    save(fullfile(s1a1dir,'embryo_data_changepoints.mat'),'chngpnts_accum');
end    

%Draw time series button panel
function draw_panel(hgui)
   myhandles = guidata(gcbo);
   ind = myhandles.ind;
   filedir = myhandles.filedir;
   nfiles = dir(filedir{ind});            
   nfiles = nfiles(~[nfiles.isdir]);  %get only the .tif files
   nsize = size(nfiles,1);
%    Delete all the panel images before redrawing
   for i = 1:length(myhandles.PanelImHands)
        delete(myhandles.PanelImHands{i});
   end
   for i=1:nsize
       nfiles(i).name = fullfile(filedir{ind}, nfiles(i).name);
   end
   chngpnts=myhandles.chngpnts;
   k = size(nfiles);
   % Update Pos and mask number displays
   currfile = nfiles(1).name;
   slashes = strfind(currfile,'\');
   currfile = currfile((slashes(end)+1):end);
   Posind = strfind(currfile,'Pos');
   Maskind = strfind(currfile,'mask');
   dashes = strfind(currfile,'_');
   PosNum = currfile(Posind+3:dashes(1)-1);
   MaskNum = currfile(Maskind+4:dashes(2)-1);
   set(h_PosNum_txt,'string',PosNum);
   set(h_MaskNum_txt,'string',MaskNum);
   n_row = ceil(sqrt(k(1)));
   n_col = floor(k(1)/n_row);
   row_count = 0;
   col_count = 0;
   for j = 1:1:k(1)
      name = nfiles(j).name;
      img = imread([nfiles(j).name]);
      img = imadjust(img);
      [xdim ydim] = size(img);
      xwidth = floor(xdim/2);
      ywidth = floor(ydim/2);
      if j == chngpnts{2}(1) | j == chngpnts{2}(2) | j == chngpnts{2}(3)
        timg = insertText(img,[round(.68*xdim),-5],j,'FontSize',32,'BoxOpacity',0.8,'TextColor','red');
      else
        timg = insertText(img,[round(.68*xdim),-5],j,'FontSize',32,'BoxOpacity',0.5,'TextColor','blue');
      end
      I2=imresize(timg, [xwidth ywidth]);
      varname=genvarname('h',who);
      foo = mod(j,n_row);
      if foo == 0
          col_count=col_count + 1;
      end
      PanelImHands{j} = uicontrol('style','pushbutton','units','pixels','position',[foo*xwidth n_col*ywidth - col_count*ywidth xwidth+10 ywidth+10],'cdata',I2,'Callback',{@pushbutton_callback,j});
   end
    myhandles.PanelImHands = PanelImHands;
    guidata(gcf,myhandles)
end

function [ chngpnts ] = changepoints( fcell,filename )
%Computes the changepoints for the transition from 1- to 2-cell, 2- to multi-cell
%and morula to blastocyst
%INPUT: fcell is {name,[1-eigratio,intensity,area,numhcpnts]}
%       filename is name of the fcell entry expected in the format
%       "batch#_pos#_cell#"
%OUTPUT: chngpnts is a cell {name,[ind1,ind2,indb]} where name is the filename and ind1,ind2 and indb
%        are the indices of the 1-cell,2-cell and blastocyst transitions 

    TF = strncmp(filename,fcell(:,1),numel(fcell{1,1})-4); %get features for cell matching filename
    foo = fcell(TF,2);
    a = [];
    eigr = [];
    numhcpnts = []; 
    for m=1:size(foo)
         a = [a,foo{m}(3)];
         eigr = [eigr,foo{m}(1)];
         numhcpnts = [numhcpnts,foo{m}(4)];
    end
    %standardize
    a = (a - mean(a))/std(a);
    eigr = (eigr - mean(eigr))/std(eigr);
    numhcpnts = (numhcpnts - mean(numhcpnts))/std(numhcpnts);
    %differentiate
    da = diff(a);
    de = diff(eigr);
    dr = diff(numhcpnts);
    %1.9-sigma criterion for 1-cell and 2-cell changepoint
    f = find(abs(de) > 1.9*std(de)); %changepoint is the first time abs(eigratio change) is > 1.9-sigma
    if length(f) >= 2 
        ind1 = f(1)+1;
        ind2 = f(2)+1; %changepoint is the first time jump is > 1.9sigma
    elseif length(f) == 1
        ind1 = 0;
        ind2 = f(1)+1;
    else
        ind1 = 0;
        ind2 = 0;
    end

    %0.55-sigma criterion for blastocyst changepoint
    f = find(abs(da) > 0.55*std(da)); %changepoint is the first time area change is > 0.55-sigma
    if ~isempty(f) 
        indb = f(1)+1;
    else
        indb = 0;
    end
    foo = [ind1,ind2,indb];
    %rearrange indices until cnodition ind1<=ind2<=indb is true
    while ~(foo(1)<=foo(2) & foo(2) <= foo(3))
        inds = randperm(3);
        foo = foo(inds);
    end
    chngpnts={filename,foo};
end




function prefix_FrNum_edit_callback(hObject, eventdata, handles)
% hObject    handle to FrDispTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrDispTxt as text
%        str2double(get(hObject,'String')) returns contents of FrDispTxt as a double

end



end