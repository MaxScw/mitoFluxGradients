function varargout = ROIsByFrameGUI(varargin)
% Simple GUI to cycle through all the frames of an acquisition and draw
% manual ROIs for each frame. Written to investigate whether there are
% systematic differences between the ICM and trophectoderm of blastocysts

% ROISBYFRAMEGUI MATLAB code for ROIsByFrameGUI.fig
%      ROISBYFRAMEGUI, by itself, creates a new ROISBYFRAMEGUI or raises the existing
%      singleton*.
%
%      H = ROISBYFRAMEGUI returns the handle to a new ROISBYFRAMEGUI or the handle to
%      the existing singleton*.
%
%      ROISBYFRAMEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROISBYFRAMEGUI.M with the given input arguments.
%
%      ROISBYFRAMEGUI('Property','Value',...) creates a new ROISBYFRAMEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROIsByFrameGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROIsByFrameGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROIsByFrameGUI

% Last Modified by GUIDE v2.5 17-Oct-2019 16:28:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROIsByFrameGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ROIsByFrameGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ROIsByFrameGUI is made visible.
function ROIsByFrameGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROIsByFrameGUI (see VARARGIN)

% Choose default command line output for ROIsByFrameGUI
handles.output = hObject;

% Select acquisition folder. 
handles.acqpath = [uigetdir '\'];
% handles.acqpath = 'C:\Dropbox\data\s1_a1_zscan';
slashes = strfind(handles.acqpath,'\'); acq = handles.acqpath(slashes(end-1)+1:end-1);

% Load indices
load([handles.acqpath '\multiD_indices.mat']);
handles.nameinds = nameinds;

% Determine which channels are present
uchans = unique(nameinds(:,4));
if ~isempty(find(strcmp(uchans,'NADH'))) Chind(1) = 1; end
if ~isempty(find(strcmp(uchans,'FAD'))) Chind(2) = 2; end
if ~isempty(find(strcmp(uchans,'UserChan'))) Chind(3) = 3; end
Chind(Chind==0)=[];
ChLabs = {'NADH','FAD','UserChan'};
if Chind(1) 
    handles.NADHbutton.Visible = 'on';
    handles.NADHbutton.Value = 1;
else
    handles.NADHbutton.Visible = 'off'; 
    handles.FADbutton.Value = 1;
end
if Chind(2) 
    handles.FADbutton.Visible = 'on';
else
    handles.FADbutton.Visible = 'off'; 
end

% Determine number of XY positions in acquisition
sdtpath = [handles.acqpath 'sorted_sdts\'];
Dpos = dir(sdtpath); Dpos(1:2)=[]; Dpos(~[Dpos.isdir])=[];
remove = [];
for i = 1:length(Dpos)
    if ~strcmp(Dpos(i).name(1:3),'Pos')
        remove = [remove i];
    end
end
Dpos(remove)=[];
handles.Poss = 0:(length(Dpos)-1);
handles.PosTxt.String = ['Pos # (' num2str(handles.Poss(1)) '-' num2str(handles.Poss(end)) ')'];

% Get Tiffs dir and frame numbers (for each Pos)
for pos = handles.Poss
    Dtiffs{pos+1} = dir([handles.acqpath '\sorted_sdts\IntTiffs_Pos' num2str(pos) '_' acq '\*.tif']);
end
handles.Dtiffs = Dtiffs;

% make an array of frame numbers
for pos = handles.Poss
    for i = 1:length(Dtiffs{pos+1}) handles.frs(i,pos+1) = str2num(Dtiffs{pos+1}(i).name(3:end-4)); end
end
% Also initialize the ManROIs cell array to have the same length
handles.ManROIs = cell(size(handles.frs,1),length(handles.Poss),2);

% Set frame text box
handles.FrTxt.String = num2str(handles.frs(1,1));

% Set slider limits to correspond with frame numbers
handles.FrSlider.Min = min(handles.frs(:,1));
handles.FrSlider.Max = max(handles.frs(:,1));
step = 1/(max(handles.frs(:,1))-min(handles.frs(:,1)));
handles.FrSlider.SliderStep = [step 5*step]; 
handles.FrSlider.Value = handles.frs(1,1);

% Load and display the first frame present
im0 = imread([Dtiffs{1}(1).folder '\' Dtiffs{1}(1).name]);
% Reshape tiff image to load channels into the 3rd dim with
% consistent indexing - NADH=1, FAD=2
% Also get number of scans
xdim = size(im0,2); ydim = size(im0,1);
handles.ydim = ydim;
im = reshape(im0,ydim,ydim,[]);

axes(handles.axes1);
set(handles.figure1,'Position',[50 10  180.6000 57.9231])
if handles.NADHbutton.Value
    nameinds = handles.nameinds;
    numscans = double(nameinds{strcmp(nameinds(:,3),'0')'&strcmp(nameinds(:,4),'NADH')'&([nameinds{:,6}]==handles.frs(1,1)),9});
    im = double(im(:,:,1))./numscans;
    handles.DispMin.String = num2str(min(min(im)));
    handles.DispMax.String = num2str(max(max(im)));
    imshow(im,[]);
elseif handles.FADbutton.Value
    nameinds = handles.nameinds;
    numscans = double(nameinds{strcmp(nameinds(:,3),'0')'&strcmp(nameinds(:,4),'FAD')'&([nameinds{:,6}]==handles.frs(1,1)),9});
    im = double(im(:,:,2))./numscans;
    handles.DispMin.String = num2str(min(min(im)));
    handles.DispMax.String = num2str(max(max(im)));
    imshow(im,[]);
end

set(handles.figure1,'Name','ROIs by frame')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ROIsByFrameGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ROIsByFrameGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function FrSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.FrSlider.Value = round(handles.FrSlider.Value);
fr = handles.FrSlider.Value;
pos = str2num(handles.PosEdit.String);

% If slider value is for a missing frame, don't do anything
if ~isempty(find(handles.frs(:,pos+1)==fr))

handles.FrTxt.String = num2str(fr);

% Replot new image
% Load and display the first frame present
im0 = imread([handles.Dtiffs{pos+1}(fr).folder '\' handles.Dtiffs{pos+1}(fr).name]);
% Reshape tiff image to load channels into the 3rd dim with
% consistent indexing - NADH=1, FAD=2
% Also get number of scans
xdim = size(im0,2); ydim = size(im0,1);
im = reshape(im0,ydim,ydim,[]);

axes(handles.axes1);
if handles.NADHbutton.Value
    % Adjust by number of FOV scans
    nameinds = handles.nameinds;
    numscans = double(nameinds{strcmp(nameinds(:,3),num2str(pos))'&strcmp(nameinds(:,4),'NADH')'&([nameinds{:,6}]==fr),9});
    im = double(im(:,:,1))./numscans;
elseif handles.FADbutton.Value
    nameinds = handles.nameinds;
    numscans = double(nameinds{strcmp(nameinds(:,3),num2str(pos))'&strcmp(nameinds(:,4),'FAD')'&([nameinds{:,6}]==fr),9});
    im = double(im(:,:,2))./numscans;
end

% Display image
if handles.AutoCont.Value
    imshow(im,[]);
    handles.DispMin.String = num2str(min(min(im)));
    handles.DispMax.String = num2str(max(max(im)));
else
    imshow(im,[str2num(handles.DispMin.String) str2num(handles.DispMax.String)]);
end

% Show previous polygon, if it exists.
try
    hold(handles.axes1,'on')
    if handles.NADHbutton.Value
        [y,x] = find(bwperim(handles.ManROIs{fr,pos+1,1}));
        plot(x,y,'b.')
    elseif handles.FADbutton.Value
        [y,x] = find(bwperim(handles.ManROIs{fr,pos+1,2}));
        plot(x,y,'g.')
    end
    hold(handles.axes1,'off')
catch
end

    
end
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function FrSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function FrTxt_Callback(hObject, eventdata, handles)
% hObject    handle to FrTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrTxt as text
%        str2double(get(hObject,'String')) returns contents of FrTxt as a double
handles.FrSlider.Value = str2num(handles.FrTxt.String);
FrSlider_Callback(handles.FrSlider, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function FrTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ROIpolygon.
function ROIpolygon_Callback(hObject, eventdata, handles)
% hObject    handle to ROIpolygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = str2num(handles.PosEdit.String);
axes(handles.axes1)
hp = impoly;
active_region = getPosition(hp);
polymask = poly2mask(active_region(:,1), active_region(:,2),handles.ydim,handles.ydim);

if handles.NADHbutton.Value
    handles.ManROIs{handles.FrSlider.Value,pos+1,1} = polymask;
elseif handles.FADbutton.Value
    handles.ManROIs{handles.FrSlider.Value,pos+1,2} = polymask;
end

% AutoSave
ManROIs = handles.ManROIs;
frs = handles.frs;
save([handles.acqpath '\ManROIs_' handles.SaveSuffix.String '.mat'],'ManROIs','frs')

% Update handles structure
guidata(hObject, handles);
FrSlider_Callback(handles.FrSlider, eventdata, handles)


function DispMin_Callback(hObject, eventdata, handles)
% hObject    handle to DispMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DispMin as text
%        str2double(get(hObject,'String')) returns contents of DispMin as a double
handles.AutoCont.Value = 0;
FrSlider_Callback(handles.FrSlider, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function DispMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DispMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DispMax_Callback(hObject, eventdata, handles)
% hObject    handle to DispMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DispMax as text
%        str2double(get(hObject,'String')) returns contents of DispMax as a double
handles.AutoCont.Value = 0;
FrSlider_Callback(handles.FrSlider, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function DispMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DispMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoCont.
function AutoCont_Callback(hObject, eventdata, handles)
% hObject    handle to AutoCont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutoCont
FrSlider_Callback(handles.FrSlider, eventdata, handles)


% --- Executes on button press in NADHbutton.
function NADHbutton_Callback(hObject, eventdata, handles)
% hObject    handle to NADHbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NADHbutton
FrSlider_Callback(handles.FrSlider, eventdata, handles)


% --- Executes on button press in FADbutton.
function FADbutton_Callback(hObject, eventdata, handles)
% hObject    handle to FADbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FADbutton
FrSlider_Callback(handles.FrSlider, eventdata, handles)



function SaveSuffix_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SaveSuffix as text
%        str2double(get(hObject,'String')) returns contents of SaveSuffix as a double


% --- Executes during object creation, after setting all properties.
function SaveSuffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveSuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadROIs.
function LoadROIs_Callback(hObject, eventdata, handles)
% hObject    handle to LoadROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    load([handles.acqpath '\ManROIs_' handles.SaveSuffix.String '.mat'])
    handles.ManROIs = ManROIs;
catch
end
% Update handles structure
guidata(hObject, handles);
FrSlider_Callback(handles.FrSlider, eventdata, handles)


% --- Executes on button press in ResetROIs.
function ResetROIs_Callback(hObject, eventdata, handles)
% hObject    handle to ResetROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Delete active file, but make a backup in case somebody hits by mistake
copyfile([handles.acqpath '\ManROIs_' handles.SaveSuffix.String '.mat'],[handles.acqpath '\ManROIs_' handles.SaveSuffix.String '_BU.mat'])
delete([handles.acqpath '\ManROIs_' handles.SaveSuffix.String '.mat'])
handles.ManROIs = cell(size(handles.frs,1),length(handles.Poss),2);
FrSlider_Callback(handles.FrSlider, eventdata, handles)


% --- Executes on button press in SaveTroph.
function SaveTroph_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTroph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save the trophectoderm, just the inverse of the ICM
for pos = handles.Poss
    for ch = 1:2
        for i = 1:size(handles.frs,1)
            % If no ICM specified, assume troph ROI is whole image
            if isempty(handles.ManROIs{i,pos+1,ch})
                ManROIs{i,pos+1,ch} = ones(handles.ydim,handles.ydim);
            else
                ManROIs{i,pos+1,ch} = ~handles.ManROIs{i,pos+1,ch};
            end
        end
    end
end
frs = handles.frs;
save([handles.acqpath '\ManROIs_troph.mat'],'ManROIs','frs')



function PosEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PosEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PosEdit as text
%        str2double(get(hObject,'String')) returns contents of PosEdit as a double
if str2num(handles.PosEdit.String)>max(handles.Poss)
    handles.PosEdit.String = num2str(max(handles.Poss));
end
if str2num(handles.PosEdit.String)<0 handles.PosEdit.String = '0'; end

FrSlider_Callback(handles.FrSlider, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function PosEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PosEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function PosTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PosTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
