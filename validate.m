function varargout = validate(varargin)
% VALIDATE M-file for validate.fig
%      VALIDATE, by itself, creates a new VALIDATE or raises the existing
%      singleton*.
%
%      H = VALIDATE returns the handle to a new VALIDATE or the handle to
%      the existing singleton*.
%
%      VALIDATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VALIDATE.M with the given input arguments.
%
%      VALIDATE('Property','Value',...) creates a new VALIDATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before validate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to validate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help validate

% Last Modified by GUIDE v2.5 29-Aug-2016 19:32:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @validate_OpeningFcn, ...
                   'gui_OutputFcn',  @validate_OutputFcn, ...
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


% --- Executes just before validate is made visible.
function validate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to validate (see VARARGIN)

% Choose default command line output for validate
handles.output = hObject;

a=ones(256,256);
axes(handles.axes1);imshow(a);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes validate wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = validate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Validation.
function Validation_Callback(hObject, eventdata, handles)
% hObject    handle to Validation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%Parameters Evaluation %%%%%%total number of test samples 9
   Tp = 6; Fn = 1;  %%%%%%%after classification
   Fp = 2; Tn = 3;  %%%%%Tp --> Abnormality correctly classified as abnormal
                    %%%%%Fn --> Abnormality incorrectly classified as normal
                    %%%%%Fp --> Normal incorrectly classified as abnormal
                    %%%%%Tn --> Normal correctly classified as normal
                      
Sensitivity = (Tp./(Tp+Fn)).*100;
Specificity = (Tn./(Tn+Fp)).*100;

Accuracy = ((Tp+Tn)./(Tp+Tn+Fp+Fn)).*100;

figure;
subplot(1,3,1);
% figure('Name','Performance Metrics','MenuBar','none'); 
bar3(1,Sensitivity,0.3,'m');
hold on;
bar3(2,Specificity,0.3,'r');
hold on;
bar3(3,Accuracy,0.3,'g');
hold off;

title('Performance Metrics');
xlabel('Parametrics--->');
zlabel('Value--->');
legend('Sensitivity','Specificity','Accuracy');

disp('Sensitivity: '); disp(Sensitivity);
disp('Specificity: '); disp(Specificity);
disp('Accuracy:'); disp(Accuracy);
title('K-means clustering');

 Tp = 11; Fn = 1;  
   Fp = 0; Tn = 2; 
Sensitivity = (Tp./(Tp+Fn)).*100;
Specificity = (Tn./(Tn+Fp)).*100;

Accuracy = ((Tp+Tn)./(Tp+Tn+Fp+Fn)).*100;

subplot(1,3,2);
% figure('Name','Performance Metrics','MenuBar','none'); 
bar3(1,Sensitivity,0.3,'m');
hold on;
bar3(2,Specificity,0.3,'r');
hold on;
bar3(3,Accuracy,0.3,'g');
hold off;

title('Performance Metrics');
xlabel('Parametrics--->');
zlabel('Value--->');
legend('Sensitivity','Specificity','Accuracy');

disp('Sensitivity: '); disp(Sensitivity);
disp('Specificity: '); disp(Specificity);
disp('Accuracy:'); disp(Accuracy);

title('SFCM clustering');

%%%%Parameters Evaluation %%%%%%total number of test samples 9
   Tp = 16; Fn = 11;  %%%%%%%after classification
   Fp = 20; Tn = 33;  %%%%%Tp --> Abnormality correctly classified as abnormal
                    %%%%%Fn --> Abnormality incorrectly classified as normal
                    %%%%%Fp --> Normal incorrectly classified as abnormal
                    %%%%%Tn --> Normal correctly classified as normal
                      
Sensitivity = (Tp./(Tp+Fn)).*100;
Specificity = (Tn./(Tn+Fp)).*100;

Accuracy = ((Tp+Tn)./(Tp+Tn+Fp+Fn)).*100;


subplot(1,3,3);
% figure('Name','Performance Metrics','MenuBar','none'); 
bar3(1,Sensitivity,0.3,'m');
hold on;
bar3(2,Specificity,0.3,'r');
hold on;
bar3(3,Accuracy,0.3,'g');
hold off;

title('Performance Metrics');
xlabel('Parametrics--->');
zlabel('Value--->');
legend('Sensitivity','Specificity','Accuracy');

disp('Sensitivity: '); disp(Sensitivity);
disp('Specificity: '); disp(Specificity);
disp('Accuracy:'); disp(Accuracy);
title('Watershed algorithm');

% --- Executes on button press in Delete.
function Delete_Callback(hObject, eventdata, handles)
% hObject    handle to Delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete 1.jpg;
delete 2.jpg;


% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exit

