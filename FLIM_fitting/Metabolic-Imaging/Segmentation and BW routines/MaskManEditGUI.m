function varargout = MaskManEditGUI(varargin)
% MASKMANEDITGUI MATLAB code for MaskManEditGUI.fig
%      MASKMANEDITGUI, by itself, creates a new MASKMANEDITGUI or raises the existing
%      singleton*.
%
%      H = MASKMANEDITGUI returns the handle to a new MASKMANEDITGUI or the handle to
%      the existing singleton*.
%
%      MASKMANEDITGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASKMANEDITGUI.M with the given input arguments.
%
%      MASKMANEDITGUI('Property','Value',...) creates a new MASKMANEDITGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MaskManEditGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MaskManEditGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MaskManEditGUI

% Last Modified by GUIDE v2.5 05-Oct-2018 12:41:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MaskManEditGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MaskManEditGUI_OutputFcn, ...
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


% --- Executes just before MaskManEditGUI is made visible.
function MaskManEditGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MaskManEditGUI (see VARARGIN)

% Choose default command line output for MaskManEditGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MaskManEditGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MaskManEditGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in finishsave.
function finishsave_Callback(hObject, eventdata, handles)
% hObject    handle to finishsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Join.
function Join_Callback(hObject, eventdata, handles)
% hObject    handle to Join (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Join
disp('HELLOOOOOO')
handles.asdf =1;
guidata(hObject, handles);


% --- Executes on button press in loadmap.
function loadmap_Callback(hObject, eventdata, handles)
% hObject    handle to loadmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
