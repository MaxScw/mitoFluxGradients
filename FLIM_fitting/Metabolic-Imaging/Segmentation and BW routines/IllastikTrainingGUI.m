function varargout = IllastikTrainingGUI(varargin)
%ILLASTIKTRAININGGUI MATLAB code file for IllastikTrainingGUI.fig
%      ILLASTIKTRAININGGUI, by itself, creates a new ILLASTIKTRAININGGUI or raises the existing
%      singleton*.
%
%      H = ILLASTIKTRAININGGUI returns the handle to a new ILLASTIKTRAININGGUI or the handle to
%      the existing singleton*.
%
%      ILLASTIKTRAININGGUI('Property','Value',...) creates a new ILLASTIKTRAININGGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to IllastikTrainingGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      ILLASTIKTRAININGGUI('CALLBACK') and ILLASTIKTRAININGGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in ILLASTIKTRAININGGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IllastikTrainingGUI

% Last Modified by GUIDE v2.5 12-Jul-2019 14:28:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IllastikTrainingGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @IllastikTrainingGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before IllastikTrainingGUI is made visible.
function IllastikTrainingGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for IllastikTrainingGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IllastikTrainingGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IllastikTrainingGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = handles.NADH.Value;
varargout{2} = handles.FAD.Value;
varargout{3} = handles.Product.Value;
varargout{4} = handles.Other.Value;
varargout{5} = handles.checkbox1.Value;
varargout{6} = handles.SepMsks.Value;
varargout{7} = str2num(handles.Threshlow.String);
varargout{8} = str2num(handles.Threshhigh.String);
delete(handles.figure1);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in Go.
function Go_Callback(hObject, eventdata, handles)
% hObject    handle to Go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargout{1} = handles.NADH.Value;
% varargout{2} = handles.FAD.Value;
% varargout{3} = handles.Product.Value;
% varargout{4} = handles.checkbox1.Value;
% close(handles.figure1);
uiresume(handles.figure1)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function Threshlow_Callback(hObject, eventdata, handles)
% hObject    handle to Threshlowtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshlowtxt as text
%        str2double(get(hObject,'String')) returns contents of Threshlowtxt as a double


% --- Executes during object creation, after setting all properties.
function Threshlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshlowtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Threshlowtxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshlowtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Threshhigh_Callback(hObject, eventdata, handles)
% hObject    handle to Threshhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threshhigh as text
%        str2double(get(hObject,'String')) returns contents of Threshhigh as a double


% --- Executes during object creation, after setting all properties.
function Threshhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Other.
function Other_Callback(hObject, eventdata, handles)
% hObject    handle to Other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Other


% --- Executes on button press in Product.
function Product_Callback(hObject, eventdata, handles)
% hObject    handle to Product (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Product


% --- Executes on button press in SepMsks.
function SepMsks_Callback(hObject, eventdata, handles)
% hObject    handle to SepMsks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SepMsks
