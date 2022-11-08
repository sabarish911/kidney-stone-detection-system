function varargout = gui_final_DWT(varargin)
% GUI_FINAL_DWT M-file for gui_final_DWT.fig
%      GUI_FINAL_DWT, by itself, creates a new GUI_FINAL_DWT or raises the existing
%      singleton*.
%
%      H = GUI_FINAL_DWT returns the handle to a new GUI_FINAL_DWT or the handle to
%      the existing singleton*.
%
%      GUI_FINAL_DWT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FINAL_DWT.M with the given input arguments.
%
%      GUI_FINAL_DWT('Property','Value',...) creates a new GUI_FINAL_DWT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_final_DWT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_final_DWT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_final_DWT

% Last Modified by GUIDE v2.5 12-Mar-2018 11:42:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_final_DWT_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_final_DWT_OutputFcn, ...
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


% --- Executes just before gui_final_DWT is made visible.
function gui_final_DWT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_final_DWT (see VARARGIN)

% Choose default command line output for gui_final_DWT
handles.output = hObject;

a=ones([256 256]);
axes(handles.axes1);imshow(a);
axes(handles.axes2);imshow(a);
axes(handles.axes5);imshow(a);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_final_DWT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_final_DWT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in inp_img.
function inp_img_Callback(hObject, eventdata, handles)
% hObject    handle to inp_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cd dataset
   [file,path] = uigetfile('*.jpg;*.bmp;*.gif;*.png', 'Pick an Image File');
   im = imread(file); 
cd ..
axes(handles.axes1);
imshow(im,[]);
imwrite(im,'fina.jpg');
handles.im = im;
save file file
% Update handles structure
guidata(hObject, handles);
% helpdlg('Test Image Selected');
% --- Executes on button press in Preprocessing.
function Preprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to Preprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
im = handles.im;

inp=imresize(im,[512 512]);
 
   if size(inp,3)>1
     inp = rgb2gray(inp);
   end
%    cd ..
   axes(handles.axes2);
   imshow(inp);
   title('Test Image');

handles.inp = inp;

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in dwt.
function dwt_Callback(hObject, eventdata, handles)
% hObject    handle to dwt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


inp = handles.inp;

[LL LH HL HH] = dwt2(inp,'db1');

aa = [LL LH;HL HH];

% % % % 2nd level decomp
[LL1 LH1 HL1 HH1] = dwt2(LL,'db1');

% aa1 = [LL1 LH1;HL1 HH1];

% % % 3rd level Decomp

[LL2 LH2 HL2 HH2] = dwt2(LL1,'db1');

% % % 4th level Decomp

[LL3 LH3 HL3 HH3] = dwt2(LL2,'db1');


v1=mean2(LL3)
v2=mean2(LH3)
v3=mean2(HL3)
v4=mean2(HH3)

aa1 = [LL3 LH3;HL3 HH3];

aa2 = [aa1 LH2;HL2 HH2];

aa3 = [aa2 LH1;HL1 HH1];
 
aa4  = [aa3 LH;HL HH];

axes(handles.axes5);
imshow(aa4,[]);
title('level-4 Wavelet Decomposition');

% % % Select the wavelet coefficients LH3 and HL3
% % % Haralick features for LH3

LH3 = uint8(LH3);
Min_val = min(min(LH3));
Max_val = max(max(LH3));
level = round(Max_val - Min_val);
GLCM = graycomatrix(LH3,'GrayLimits',[Min_val Max_val],'NumLevels',level);
stat_feature = graycoprops(GLCM);
Energy_fet1 = stat_feature.Energy;
Contr_fet1 = stat_feature.Contrast;
Corrla_fet1 = stat_feature.Correlation;
Homogen_fet1 = stat_feature.Homogeneity;

% % % % % Entropy
        R = sum(sum(GLCM));
        Norm_GLCM_region = GLCM/R;
        
        Ent_int = 0;
        for k = 1:length(GLCM)^2
            if Norm_GLCM_region(k)~=0
                Ent_int = Ent_int + Norm_GLCM_region(k)*log2(Norm_GLCM_region(k));
            end
        end
        Entropy_fet1 = -Ent_int;

%%%%%Haralick Features For HL3        
HL3 = uint8(HL3);
Min_val = min(min(HL3));
Max_val = max(max(HL3));
level = round(Max_val - Min_val);
GLCM = graycomatrix(HL3,'GrayLimits',[Min_val Max_val],'NumLevels',level);
stat_feature = graycoprops(GLCM);
Energy_fet2 = stat_feature.Energy;
Contr_fet2 = stat_feature.Contrast;
Corrla_fet2= stat_feature.Correlation;
Homogen_fet2 = stat_feature.Homogeneity;
% % % % % Entropy
        R = sum(sum(GLCM));
        Norm_GLCM_region = GLCM/R;
        
        Ent_int = 0;
        for k = 1:length(GLCM)^2
            if Norm_GLCM_region(k)~=0
                Ent_int = Ent_int + Norm_GLCM_region(k)*log2(Norm_GLCM_region(k));
            end
        end
% % % % % % Ent_int = entropy(GLCM);
        Entropy_fet2 = -Ent_int;

%%%%% Feature Sets

F1 = [Energy_fet1 Contr_fet1 Corrla_fet1 Homogen_fet1 Entropy_fet1];
F2 = [Energy_fet2 Contr_fet2 Corrla_fet2 Homogen_fet2 Entropy_fet2];

qfeat = [F1 F2]';
save qfeat qfeat;

disp('Query Features: ');
disp(qfeat);


handles.aa4=aa4;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in database.
function database_Callback(hObject, eventdata, handles)
% hObject    handle to database (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nnlearn;


% --- Executes on button press in nn_classifier.
function nn_classifier_Callback(hObject, eventdata, handles)
% hObject    handle to nn_classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inp = handles.inp;
% a2=handles.im;
% m2=handles.aa4;
% a2=imresize(a2,[64 64]);
% m2=imresize(m2,[64 64]);
load qfeat;
load dfeatures;
load netp;

load qfeat;
load netp;

qfeatx=max(max(qfeat));
%%%%%%classification

cout = sim(netp,qfeat);
cout = vec2ind(cout);

if isequal(cout,1)
warndlg('NORMAL')
elseif isequal(cout,2)
warndlg('stone 1');
% output = segment(inp);
elseif isequal(cout,3)
warndlg('stone 2');
elseif isequal(cout,4)
warndlg('stone 3');
elseif isequal(cout,5)
warndlg('stone 4');
else 
   helpdlg('Db updation required')
end  
handles.cout=cout;
 guidata(hObject,handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cout=handles.cout;
inp=handles.inp;

if cout>1 
cd datafin
a=inp;
k=0;
for  m =1:128:size(a,1)
    for n=1:128:size(a,2)
        k=k+1;
        imnew(:,:,k)=a(m:m+127,n:n+127);
        blknam=strcat(num2str(k),'.jpg');
        imwrite(imnew(:,:,k),strcat(num2str(k),'.jpg'));
    end
end
cd ..

for i=1:16
    blknam=strcat(num2str(i),'.jpg');
    cd datafin
    immg2(:,:,i)=imread(blknam);
    cd ..
I=immg2(:,:,i);
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
% figure(11)
% imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
L = watershed(gradmag);
Lrgb = label2rgb(L);
% figure(12), imshow(Lrgb), title('Watershed transform of gradient magnitude')
se = strel('disk', 20);
Io = imopen(I, se);
% figure
% imshow(Io), title('Opening (Io)')
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
% figure
% imshow(Iobr), title('Opening-by-reconstruction')
Ioc = imclose(Io, se);
% figure(13)
% imshow(Ioc), title('Opening-closing (Ioc)')
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
% figure(14)
% imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
fgm = imregionalmax(Iobrcbr);
% figure(15)
% imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
I2 = I;
I2(fgm) = 255;
% figure(16)
% imshow(I2), title('Regional maxima superimposed on original image (I2)')
se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
% figure(17)
% imshow(I3)
% title('Modified regional maxima superimposed on original image')
 bw = im2bw(Iobrcbr);
% imshow(bw), title('Thresholded opening-closing by reconstruction')
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
 figure(19)
% imshow(bgm), title('Watershed ridge lines')
gradmag2 = imimposemin(gradmag, bgm | fgm4);
I4 = I;
I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
% figure
% imshow(I4)
% title('Markers and object boundaries superimposed on fused image')
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
end

region=connectedregions(rgb2gray(Lrgb),16);
val=im2bw(I,.8);
I=region;
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
% figure(11)
% imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
L = watershed(gradmag);
Lrgb = label2rgb(L);
% figure(12), imshow(Lrgb), title('Watershed transform of gradient magnitude')
se = strel('disk', 20);
Io = imopen(I, se);
% figure
% imshow(Io), title('Opening (Io)')
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
% figure
% imshow(Iobr), title('Opening-by-reconstruction')
Ioc = imclose(Io, se);
% figure(13)
% imshow(Ioc), title('Opening-closing (Ioc)')
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
% figure(14)
% imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
fgm = imregionalmax(Iobrcbr);
% figure(15)
% imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
I2 = I;
I2(fgm) = 255;
% figure(16)
% imshow(I2), title('Regional maxima superimposed on original image (I2)')
se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
% figure(17)
% imshow(I3)
% title('Modified regional maxima superimposed on original image')
bw = im2bw(Iobrcbr);
% figure(18)
% imshow(bw), title('Thresholded opening-closing by reconstruction')
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
% figure(19)
% imshow(bgm), title('Watershed ridge lines')
gradmag2 = imimposemin(gradmag, bgm | fgm4);
I4 = I;
I4(imerode(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
% figure
% imshow(I4)
% title('Markers and object boundaries superimposed on fused image')
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure(20)
% imshow(Lrgb)
grayimg=rgb2gray(Lrgb);
% title('Colored watershed label matrix')
% figure
% imshow(I)
% himage = imshow(Lrgb);
% himage.AlphaData = 0.3;
% title('superimposed transparently on fused image')
output=im2bw(grayimg,.8);
boundary = bwboundaries(output);
axes(handles.axes1);

%Boundary Label the Filtered Image
    [L num]=bwlabel(output);
    STATS=regionprops(L,'all');
    cc=[];
    removed=0;
    hold on
    
%Remove the noisy regions 
   for i=1:num
    dd=STATS(i).Area;
    if (dd < 150)
    L(L==i)=0;
    removed = removed + 1;
    num=num-1;
    else
    end
   end
    figure
imshow(removed)
hold on;
for ii=1:1:length(boundary)
    btemp = boundary{ii};
    plot(btemp(:,2),btemp(:,1),'r','LineWidth',2);
end
hold off;
title('Segmented Image');
countt=0;
 for iii=1:size(output,1)
     for jjj=1:size(output,2)
         if output(iii,jjj)==0
             countt=countt+1;
         end
     end
 end
 disp('Area of K-Means segmented region is:');
 disp(countt);
perimeter=sqrt(length(boundary))*.264;
disp(perimeter);
end



% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cout=handles.cout;
inp=handles.inp;
save cout cout
a=inp;
figure;
if cout>1 
cd datafin
% if isequal(cout,2) 
k=0;
for  m =1:128:size(a,1)
    for n=1:128:size(a,2)
        k=k+1;
        imnew(:,:,k)=a(m:m+127,n:n+127);
        blknam=strcat(num2str(k),'.jpg');
        imwrite(imnew(:,:,k),strcat(num2str(k),'.jpg'));
    end
end
cd ..

for i=1:16
    subplot(4,4,i); 
    blknam=strcat(num2str(i),'.jpg');
    cd datafin
    immg2(:,:,i)=imread(blknam);
    cd ..
    imshow(immg2(:,:,i));
end
% end
 output = ksegment(immg2);
 axes(handles.axes1);
 imshow(output,[]);
 boundary = bwboundaries(output);
 for ii=1:1:length(boundary)
    btemp = boundary{ii};
    plot(btemp(:,2),btemp(:,1),'r','LineWidth',2);
end
hold on;
 countt=0;
 for iii=1:size(output,1)
     for jjj=1:size(output,2)
         if output(iii,jjj)==255
             countt=countt+1;
         end
     end
 end
area=countt;
disp('IDENTIFIED KMEANS AFFECTED REGION')
disp(area);
disp('percentage of affection')
percentagefin=(area/512*512)*100;
disp(percentagefin);

perimeter=length(boundary)*2*3.14;
disp('Perimeter of region:')
disp(perimeter);
end

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cout=handles.cout;
inp=handles.inp;
a=inp;
if cout>1 
cd datafin
% if isequal(cout,2) 
k=0;
for  m =1:128:size(a,1)
    for n=1:128:size(a,2)
        k=k+1;
        imnew(:,:,k)=a(m:m+127,n:n+127);
        blknam=strcat(num2str(k),'.jpg');
        imwrite(imnew(:,:,k),strcat(num2str(k),'.jpg'));
    end
end
cd ..

for i=1:16
    subplot(4,4,i); 
    blknam=strcat(num2str(i),'.jpg');
    cd datafin
    immg2(:,:,i)=imread(blknam);
    cd ..
    imshow(immg2(:,:,i));
end
% end
 output = segment(immg2);
 axes(handles.axes1);
 imshow(output,[]);
 boundary = bwboundaries(output);
 for ii=1:1:length(boundary)
    btemp = boundary{ii};
    plot(btemp(:,2),btemp(:,1),'r','LineWidth',2);
end
hold on;
 countt=0;
 for iii=1:size(output,1)
     for jjj=1:size(output,2)
         if output(iii,jjj)==255
             countt=countt+1;
         end
     end
 end
area=countt;
disp('IDENTIFIED SFCM AFFECTED REGION')
disp(area);
disp('percentage of affection')
percentagefin=(area/512*512)*100;
disp(percentagefin);

perimeter=length(boundary)*2*3.14;
disp('Perimeter of region:')
disp(perimeter);
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
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