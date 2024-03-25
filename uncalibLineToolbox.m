%     << Automatic Line-Image Extraction Toolbox for Uncalibrated Central
%        Systems with Revolution Symmetry (release v0.5 alpha) >> 
%     ====================================================================
%     Copyright (C) 2014  Jesus Bermudez-Cameo
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%     Acknowledgement:
%     This code has been implemented by Jesus Bermudez-Cameo and supervised 
%     by G.Lopez-Nicolas and J.J. Guerrero
% 
%     One may refer to the following paper:
% 
%     “Line extraction in uncalibrated central images with revolution symmetry”,
%     J. Bermudez-Cameo, G. Lopez-Nicolas and J. J. Guerrero,  
%     24th British Machine Vision Conference, BMVC, pp 1-11; Bristol, UK, Sept.  2013, 
% 
%     This work has been supported by the University of Zaragoza , the 
%     Spanish project VINEA DPI2012-31781, DGA-FSE(T04 and FEDER funds.
%     Jesus Bermudez-Cameo was supported by the FPU program AP2010-3849.

function varargout = uncalibLineToolbox(varargin)
%UNCALIBLINETOOLBOX M-file for uncalibLineToolbox.fig
%      UNCALIBLINETOOLBOX, by itself, creates a new UNCALIBLINETOOLBOX or raises the existing
%      singleton*.
%
%      H = UNCALIBLINETOOLBOX returns the handle to a new UNCALIBLINETOOLBOX or the handle to
%      the existing singleton*.
%
%      UNCALIBLINETOOLBOX('Property','Value',...) creates a new UNCALIBLINETOOLBOX using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to uncalibLineToolbox_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      UNCALIBLINETOOLBOX('CALLBACK') and UNCALIBLINETOOLBOX('CALLBACK',hObject,...) call the
%      local function named CALLBACK in UNCALIBLINETOOLBOX.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help uncalibLineToolbox

% Last Modified by GUIDE v2.5 22-Oct-2014 18:01:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @uncalibLineToolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @uncalibLineToolbox_OutputFcn, ...
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


% --- Executes just before uncalibLineToolbox is made visible.
function uncalibLineToolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for uncalibLineToolbox
handles.output = hObject;






set(handles.cannySlider, 'Min', 0);
set(handles.cannySlider, 'Max', 10);

handles.busyFlag = 0;
handles.plotFlag = 0;
handles.externFigure = 1; %Do not change this flag

handles.nAttempts = 100; %% It depends on ...
handles.nAttempts2P = 50; %% It depends on ...
handles.nAttemptsSingleBound = 50000;

handles.ratioShowLineImages = 1;

handles.preProcessingParams.cannyIndex = 2.5;
handles.preProcessingParams.acumDeltaThetaRatio = 0.5; % Ratio of Acumulated Histogram
handles.preProcessingParams.sizeThresh = [0.5;100;3000]; 
handles.preProcessingParams.subBoundThresh = 0.2;
handles.preProcessingParams.filterOmega = 2*pi*150;


handles.lineExtractionParams.normalizedShortPixelThreshold = 0.001;
handles.lineExtractionParams.normalizedrThreshold = 0.002;


handles.calibratedExtract.params.normalizedrThreshold = 0.002;
handles.calibratedExtract.params.aThreshold = 45/180*pi;
handles.calibratedExtract.params.thresholdInWhile = 0.3;
handles.calibratedExtract.params.normalizedSizeThreshold = 0.1242;
handles.calibratedExtract.params.maxLIperBound = 5;
handles.calibratedExtract.params.nAttempts = 100;


set(handles.cannySlider, 'Value', handles.preProcessingParams.cannyIndex);
set(handles.cannyIndexEdit,'string',sprintf('%f',handles.preProcessingParams.cannyIndex));

set(handles.proportionalSlider, 'Value', handles.preProcessingParams.sizeThresh(1));
set(handles.proportionalEdit,'string',sprintf('%f',handles.preProcessingParams.sizeThresh(1)));

set(handles.inferiorSlider, 'Value', handles.preProcessingParams.sizeThresh(2));
set(handles.inferiorEdit,'string',sprintf('%f',handles.preProcessingParams.sizeThresh(2)));

set(handles.superiorSlider, 'Value', handles.preProcessingParams.sizeThresh(3));
set(handles.superiorEdit,'string',sprintf('%f',handles.preProcessingParams.sizeThresh(3)));

set(handles.subBoundSlider, 'Value', handles.preProcessingParams.subBoundThresh);
set(handles.subBoundEdit,'string',sprintf('%f',handles.preProcessingParams.subBoundThresh));

set(handles.filterOmegaSlider, 'Value', handles.preProcessingParams.filterOmega);
set(handles.filterOmegaEdit,'string',sprintf('%f',handles.preProcessingParams.filterOmega));

set(handles.reStractSlider, 'Value', handles.ratioShowLineImages);
set(handles.ratioOfExBounds,'string',sprintf('%f',handles.ratioShowLineImages));





% Update handles structure
guidata(hObject, handles);

path(path,'edgesDetection\');
path(path,'preprocessing\');
path(path,'independentExtraction\');
path(path,'frameworkGeneral\');
path(path,'parametricLineImages\');
path(path,'projectionModel\');
path(path,'calibratedExtraction\');
% path(path,'optimPrincipalPoint\');
path(path,'hyperDedicated\');
path(path,'genericModelExtraction\');
path(path,'guideTools\');
path(path,'manualTools\');
path(path,'tools\');

path(path,'optimizations\');

% UIWAIT makes uncalibLineToolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = uncalibLineToolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function uipanel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function u0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to u0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of u0Edit as text
%        str2double(get(hObject,'String')) returns contents of u0Edit as a double
u_0 = str2num(get(handles.u0Edit,'string'));
if numel(u_0)~=1
    handles.u_0 = u_0;
end
figure(handles.figure1);



% --- Executes during object creation, after setting all properties.
function u0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to u0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v0Edit_Callback(hObject, eventdata, handles)
% hObject    handle to v0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v0Edit as text
%        str2double(get(hObject,'String')) returns contents of v0Edit as a double
v_0 = str2num(get(handles.v0Edit,'string'));
if numel(v_0)~=1
    handles.v_0 = v_0;
end
figure(handles.figure1);



% --- Executes during object creation, after setting all properties.
function v0Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v0Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in loadImageButton.
function loadImageButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.busyFlag
    ls('*.jpg');
    imgName = input('Input the name of the image:  ','s');

    if exist(imgName,'file')
        imgIn=imread(imgName);
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
        imshow(imgIn);    
        hold on;
        handles.imgIn=imgIn;    
        handles.imgGray = rgb2gray(imgIn);        
        [handles.rows,handles.columns,~] = size(imgIn);
        handles.mask = zeros(handles.rows,handles.columns);
        handles.u_0 = handles.columns/2;
        handles.v_0 = handles.rows/2;
        set(handles.ppMouseButton,'Enable','on');
        set(handles.cannyButton,'Enable','on');
        set(handles.preprocessingButton,'Enable','on');    
        
        
        
        
        
        guidata(hObject,handles);
    else
        disp('Unable to open the file');
    end
    
    iDot = find(imgName=='.');
    maskName = strcat(imgName(1:iDot(end)),'pbm');
    if exist(maskName,'file')
        fprintf('There is a file containing the mask for this image.\n');        
        reply = input('Do you want to load it? Y/N [Y]: ', 's');
        if isempty(reply)
            reply = 'Y';
        end
        
        if reply == 'Y'
            mask=imread(maskName);
            handles.mask=mask;
            
            [iCircle,jCircle] = find(edge(mask,'sobel')==true);                        
            [u_0,v_0] = getCircleFromNPoints(jCircle,iCircle,handles.rows);            
            handles.u_0 = u_0;
            handles.v_0 = v_0;                       
            set(handles.u0Edit,'string',sprintf('%f',u_0));
            set(handles.v0Edit,'string',sprintf('%f',v_0));  
            
            guidata(hObject,handles);                        
        end
        
    end
    
    
end
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function loadImageButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in ppMouseButton.
function ppMouseButton_Callback(hObject, eventdata, handles)
% hObject    handle to ppMouseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.busyFlag

    if isfield(handles,'imgIn')
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
        [m,n,~] = size(handles.imgIn);
        plot(n/2,m/2,'g*');
        [u_0,v_0] = ginput(1);
        handles.u_0 = u_0;
        handles.v_0 = v_0;
        imshow(handles.imgIn);
        line([0 n],[v_0,v_0],'Color','b');
        line([u_0 u_0],[0,m],'Color','b');
        set(handles.u0Edit,'string',sprintf('%f',u_0));
        set(handles.v0Edit,'string',sprintf('%f',v_0));  
        
        set(handles.genMaskButton,'Enable','on');   
        
        guidata(hObject,handles);
    else
        disp('There is not input image');
    end
end
figure(handles.figure1);

% --- Executes on button press in loadMaskButton.
function loadMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.busyFlag
    ls('*.pbm');
    imgName = input('Input the name of the mask:  ','s');

    if exist(imgName,'file')
        mask=imread(imgName);
        handles.mask=mask;
        
        [iCircle,jCircle] = find(edge(mask,'sobel')==true);                        
        [u_0,v_0] = getCircleFromNPoints(jCircle,iCircle,handles.rows);            
        handles.u_0 = u_0;
        handles.v_0 = v_0;                       
        set(handles.u0Edit,'string',sprintf('%f',u_0));
        set(handles.v0Edit,'string',sprintf('%f',v_0));          
        
        guidata(hObject,handles);
    else
        disp('Unable to open the mask file');
    end
end
figure(handles.figure1);



% --- Executes on button press in saveMaskButton.
function saveMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.busyFlag
    
    maskName = input('Input the name of the mask file without extension:  ','s');    
    maskName = sprintf('%s.pbm',maskName);
    
    if isfield(handles,'mask')        
        imwrite(handles.mask,maskName,'pbm');
        guidata(hObject,handles);    
    end
end
figure(handles.figure1);


% --- Executes on button press in genMaskButton.
function genMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to genMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.busyFlag
    if isfield(handles,'imgIn')
        
        numberOfCircles = input('Number of circles:  ');    

        switch numberOfCircles
            case 1
                if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
                [x_1,y_1] = ginput(1);
                r = norm([x_1-handles.u_0;y_1-handles.v_0]);
                theta = -pi:pi/360:pi;
                x = handles.u_0 + r*cos(theta);
                y = handles.v_0 + r*sin(theta);
                plot(x,y,'g-');

                [jImg,iImg] = meshgrid(1:1:size(handles.imgIn,2), 1:1:size(handles.imgIn,1));           
                jImg = reshape(jImg,numel(jImg),1); iImg = reshape(iImg,numel(iImg),1);        
                rImg = sqrt((jImg-handles.u_0).^2+(iImg-handles.v_0).^2);

                handles.mask = reshape(rImg>r,size(handles.imgIn,1),size(handles.imgIn,2));

                figure(2);
                imshow(handles.mask);        

                set(handles.saveMaskButton,'Enable','on');           
                guidata(hObject,handles);                
            case 2
                if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
                [x_1,y_1] = ginput(1);
                r1 = norm([x_1-handles.u_0;y_1-handles.v_0]);
                theta = -pi:pi/360:pi;
                x1 = handles.u_0 + r1*cos(theta);
                y1 = handles.v_0 + r1*sin(theta);
                plot(x1,y1,'g-');
                
                [x_2,y_2] = ginput(1);
                r2 = norm([x_2-handles.u_0;y_2-handles.v_0]);
                x2 = handles.u_0 + r2*cos(theta);
                y2 = handles.v_0 + r2*sin(theta);
                plot(x2,y2,'g-');
                

                [jImg,iImg] = meshgrid(1:1:size(handles.imgIn,2), 1:1:size(handles.imgIn,1));           
                jImg = reshape(jImg,numel(jImg),1); iImg = reshape(iImg,numel(iImg),1);        
                rImg = sqrt((jImg-handles.u_0).^2+(iImg-handles.v_0).^2);

                handles.mask = reshape((rImg>r1)|(rImg<r2),size(handles.imgIn,1),size(handles.imgIn,2));

                figure(2);
                imshow(handles.mask);        

                set(handles.saveMaskButton,'Enable','on');           
                guidata(hObject,handles);                
            otherwise
        end                        
    end   
end
figure(handles.figure1);


% --- Executes on button press in circleDetButton.
function circleDetButton_Callback(hObject, eventdata, handles)
% hObject    handle to circleDetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.busyFlag
    if isfield(handles,'imgGray'), 
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
        text(50,50,'Busy','Color','r') ;
        drawnow;
        figure(handles.figure1);
        handles.busyFlag = 1;
        [edges, ~, ~] = canny(handles.imgGray,1,[],[1 1]*handles.preProcessingParams.cannyIndex);
        [iImg,jImg] =find(edges==1);
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;                
        imshow(handles.imgIn);    
        hold on;
        plot(jImg,iImg,'r.','MarkerSize',2);
        handles.busyFlag = 0;
  
        xSel = zeros(3,3);
        for k = 1:3
            [x_1,y_1] = ginput(1);
            iii=getCloser([x_1;y_1],[jImg,iImg]');
            xSel(:,k) = [jImg(iii(1));iImg(iii(1));1];
            plot(xSel(1,k),xSel(2,k),'g*');
        end
                
        [u_0,v_0] = getCircleFromNPoints(xSel(1,1:3)',xSel(2,1:3)',handles.rows);            

        
        handles.u_0 = u_0;
        handles.v_0 = v_0;
        imshow(handles.imgIn);
        line([0 handles.columns],[v_0,v_0],'Color','b');
        line([u_0 u_0],[0,handles.rows],'Color','b');
        set(handles.u0Edit,'string',sprintf('%f',u_0));
        set(handles.v0Edit,'string',sprintf('%f',v_0));  
        
        set(handles.genMaskButton,'Enable','on');   
        
        guidata(hObject,handles);
        
        
    end
end
figure(handles.figure1);


% --- Executes on key press with focus on hyperRadButton and none of its controls.
function hyperRadButton_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to hyperRadButton (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function paraRadButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paraRadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called















function cannyIndexEdit_Callback(hObject, eventdata, handles)
% hObject    handle to cannyIndexEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cannyIndexEdit as text
%        str2double(get(hObject,'String')) returns contents of cannyIndexEdit as a double

cannyIndex = str2double(get(hObject,'String'));
handles.preProcessingParams.cannyIndex = cannyIndex;
guidata(hObject,handles);
figure(handles.figure1);



% --- Executes during object creation, after setting all properties.
function cannyIndexEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cannyIndexEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









% --- Executes on slider movement.
function cannySlider_Callback(hObject, eventdata, handles)
% hObject    handle to cannySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

cannyIndex = get(hObject,'Value');
handles.preProcessingParams.cannyIndex = cannyIndex;
set(handles.cannyIndexEdit,'string',sprintf('%f',cannyIndex));
guidata(hObject,handles);
figure(handles.figure1);



% --- Executes during object creation, after setting all properties.
function cannySlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cannySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in cannyButton.
function cannyButton_Callback(hObject, eventdata, handles)
% hObject    handle to cannyButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.busyFlag
    if isfield(handles,'imgGray'), 
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
        text(50,50,'Busy','Color','r') ;
        drawnow;
        figure(handles.figure1);
        handles.busyFlag = 1;
        [edges, ~, ~] = canny(handles.imgGray,1,[],[1 1]*handles.preProcessingParams.cannyIndex);
        [iImg,jImg] =find((edges&~handles.mask)==1);
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;                
        imshow(handles.imgIn);    
        hold on;
        plot(jImg,iImg,'r.','MarkerSize',2);
        handles.busyFlag = 0;
  
        set(handles.circleDetButton,'Enable','on');          
        guidata(hObject,handles);                           
    end
end
figure(handles.figure1);


% --- Executes on button press in preprocessingButton.
function preprocessingButton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessingButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.busyFlag

    if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
    imshow(handles.imgIn);
    hold on;
    
    if isfield(handles,'imgIn')

        text(50,50,'Busy','Color','r') ;
        drawnow;
        figure(handles.figure1);
        
        handles.busyFlag = 1;
        disp('Preprocessing ... Please wait'); 
        tic;
        [handles.boundaries,handles.subBounds,handles.sizeSubBounds,handles.kEnd,handles.boundary] = imagePreprocessingBM(handles.imgIn,handles.mask,handles.preProcessingParams,handles.u_0,handles.v_0);        
        disp('Preprocessing done ...');    
        toc;
        handles.busyFlag = 0;

        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
        imshow(handles.imgIn);

        handles.boundsIndexMap = zeros(handles.rows,handles.columns);        
        for k = 1:numel(handles.subBounds)
            xImg = handles.subBounds(k).xImg;
            handles.boundsIndexMap(coord2Index(xImg(2,:),xImg(1,:),handles.rows)) = k;
            xText = mean(xImg(1,:));
            yText = mean(xImg(2,:));
            text(xText,yText,sprintf('%i',k),'Color','r');           
            plot(xImg(1,:),xImg(2,:),'r.','MarkerSize',2);                                                
        end       
        
        set(handles.launchButton,'Enable','on');   

        set(handles.manual3PointsButton,'Enable','on'); %%NEW NEW
        set(handles.singleBoundButton,'Enable','on'); %%NEW NEW        
        
        guidata(hObject,handles);   
    else
       disp(' There is not input image or mask');    
    end
end
figure(handles.figure1);






% --- Executes on button press in launchButton.
function launchButton_Callback(hObject, eventdata, handles)
% hObject    handle to launchButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global secondOrderFlag;

if ~handles.busyFlag


    if isfield(handles,'boundaries')

        modelClass = getModelClassFromHandle(handles);
        algorithm = getAlgorithmFromHandle(handles);
        
        if algorithm == 3
            secondOrderFlag = true;
        else
            secondOrderFlag = false;
        end
        
        if(modelClass == 0)
            disp('Automatic device detection not yet implemented');        
            return;
        end;
        
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
        imshow(handles.imgIn);
        hold on;
        text(50,50,'Busy','Color','r') ;
        drawnow;
        figure(handles.figure1);
        
        handles.busyFlag = 1;
        disp('Extracting lines ... Please Wait');        
        
                        
        if algorithm == 1
            [segmentsAcum,n,calib] = uncalibratedExtractionFromBoundaries3PBasedBM(handles.boundaries,handles.subBounds,handles.sizeSubBounds,handles.kEnd,handles.nAttempts,modelClass,handles.lineExtractionParams,handles.rows,handles.u_0,handles.v_0,handles.plotFlag);                                                          
        else
            [segmentsAcum,n,calib] = uncalibratedExtractionFromBoundaries2PBasedBM(handles.boundaries,handles.subBounds,handles.sizeSubBounds,handles.kEnd,handles.nAttempts2P,modelClass,handles.lineExtractionParams,handles.rows,handles.u_0,handles.v_0,handles.plotFlag);
%             disp('Gradient algorithm not yet implemented');
%             return;               
        end                
        disp('Extracting lines done');                        
        
        if get(handles.optimCheckBox,'Value') && isreal(calib) && all(~isnan(calib)) && all(~isinf(calib)) && numel(segmentsAcum)>0                    
            if modelClass == 2
                disp('Optimization for hypercatadioptric systems not yet released');
%                Op_optim = optimFocalAndLIsHyper(segmentsAcum,n,calib(1),calib(2));
%                calib = Op_optim(1:2);                                                        
            else
                Op_optim = optimFocalAndLIs(segmentsAcum,n,calib(1),modelClass);
                calib = Op_optim(1);                         
            end
        end
        
        
        handles.busyFlag = 0;              

        
        set(handles.saveResultsButton,'Enable','on');
        set(handles.reExtractButton,'Enable','on');

        handles.calib = calib;
        handles.segments = segmentsAcum;
        handles.n = n;
               
        guidata(hObject,handles);  
        
        reExtractButton_Callback(hObject, eventdata, handles);
    else
       disp(' There is not input image or mask');    
    end
end
figure(handles.figure1);


% --- Executes on button press in reExtractButton.
function reExtractButton_Callback(hObject, eventdata, handles)
% hObject    handle to reExtractButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~handles.busyFlag
    
    if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
    imshow(handles.imgIn);
    hold on;
    text(50,50,'Busy','Color','r') ;
    drawnow;
    figure(handles.figure1);

    handles.busyFlag = 1;    
    
    
    modelClass = getModelClassFromHandle(handles);
    u_0 = handles.u_0;
    v_0 = handles.v_0;
    calib = handles.calib;
    rVL = getRVLFromCalib(calib,modelClass);

    %% re-Extraction using calibration.    
    rThreshold = handles.calibratedExtract.params.normalizedrThreshold*handles.rows; 
    thresholdInWhile = handles.calibratedExtract.params.thresholdInWhile;
    maxLIperBound = handles.calibratedExtract.params.maxLIperBound;
    nAttempts = handles.calibratedExtract.params.nAttempts;
    
    T_center = [[1 0 -u_0];[0 1 -v_0];[0 0 1]];
    
    [acceptedBinMatrix,nCalibrated,acceptedBinMatrixSizes] = calibratedLineImageExtraction(handles.subBounds,modelClass,calib,rThreshold,nAttempts,maxLIperBound,thresholdInWhile,handles.u_0,handles.v_0);    
    
    handles.calibratedExtract.acceptedBinMatrix = acceptedBinMatrix;
    handles.calibratedExtract.nCalibrated = nCalibrated;
    handles.calibratedExtract.acceptedBinMatrixSizes = acceptedBinMatrixSizes;
    
    
    n = [];
    sizeSegments = [];
    segments = cell(1,1);
    ik = 1;
    for k = 1:numel(handles.subBounds)        
        for kk = 1:acceptedBinMatrixSizes(1,k)
            segment = handles.subBounds(k).xHat(:,acceptedBinMatrix{k}(kk,:));
            localn = nCalibrated{k}(:,kk);          
            
            segments{ik} = segment;
            n = [n,localn];
            sizeSegments = [sizeSegments,size(segments{ik},2)];
            ik = ik+1;
        end                
    end
    
    segments = segments(sizeSegments>0);
    n = n(:,sizeSegments>0);
    sizeSegments = sizeSegments(sizeSegments>0);
    
    handles.calibratedExtract.segments = segments;
    handles.calibratedExtract.sizeSegments = sizeSegments;
    handles.calibratedExtract.n = n;
                            
    deltaThetaSize = zeros(1,numel(segments));
    for k = 1:numel(segments)
        deltaTheta=getDeltaTheta(segments{k}(1,:),segments{k}(2,:));      
        deltaThetaSize(k)=deltaTheta*sizeSegments(k);         
    end
    
    handles.calibratedExtract.deltaThetaSize = deltaThetaSize;
   
    set(handles.showReExButton,'Enable','on');   

    if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
    imshow(handles.imgIn);
    hold on;
    
    plotLineImages(segments,n,calib,modelClass,u_0,v_0);

    handles.busyFlag = 0;
    guidata(hObject,handles);

end
figure(handles.figure1);



% --- Executes on button press in showReExButton.
function showReExButton_Callback(hObject, eventdata, handles)
% hObject    handle to showReExButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.busyFlag    
              
    
    if get(handles.kappaRadButton,'Value'), 
        [~,sortIndex] = sort(handles.calibratedExtract.deltaThetaSize,'descend');        
    else
        [~,sortIndex] = sort(handles.calibratedExtract.sizeSegments,'descend');
    end;
            
    sortIndex = sortIndex(1:ceil(numel(sortIndex)*handles.ratioShowLineImages));       
    segments = handles.calibratedExtract.segments(sortIndex);
    n = handles.calibratedExtract.n(:,sortIndex);
    
    
    modelClass = getModelClassFromHandle(handles);
    u_0 = handles.u_0;
    v_0 = handles.v_0;
    calib = handles.calib;
    rVL = getRVLFromCalib(calib,modelClass);    
    
    if handles.externFigure, figure(3); else axes(handles.mainFigure); end;
    imshow(handles.imgIn);
    hold on;
    plotLineImages(segments,n,calib,modelClass,u_0,v_0);




end
figure(handles.figure1);







% --- Executes on button press in saveResultsButton.
function saveResultsButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveResultsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~handles.busyFlag
    
    resultsName = input('Input the name of the results file without extension:  ','s');    
    resultsName = sprintf('%s.mat',resultsName);    
    
    subBounds = handles.subBounds;    
    calibratedExtract = handles.calibratedExtract;    
    modelClass = getModelClassFromHandle(handles);
    u_0 = handles.u_0;
    v_0 = handles.v_0;
    calib = handles.calib;
    rVL = getRVLFromCalib(calib,modelClass);
    modelClass = getModelClassFromHandle(handles);
    

    save(resultsName,'subBounds','calibratedExtract','u_0','v_0','calib','rVL','modelClass');    
end 
figure(handles.figure1);



function proportionalEdit_Callback(hObject, eventdata, handles)
% hObject    handle to proportionalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of proportionalEdit as text
%        str2double(get(hObject,'String')) returns contents of proportionalEdit as a double
handles.preProcessingParams.sizeThresh(1) = str2double(get(hObject,'String'));
guidata(hObject,handles);
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function proportionalEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to proportionalEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function proportionalSlider_Callback(hObject, eventdata, handles)
% hObject    handle to proportionalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.preProcessingParams.sizeThresh(1) = get(hObject,'Value');
set(handles.proportionalEdit,'string',sprintf('%f',handles.preProcessingParams.sizeThresh(1)));
guidata(hObject,handles);
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function proportionalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to proportionalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function inferiorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to inferiorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of inferiorEdit as text
%        str2double(get(hObject,'String')) returns contents of inferiorEdit as a double

handles.preProcessingParams.sizeThresh(2) = round(str2double(get(hObject,'String')));
guidata(hObject,handles);
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function inferiorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inferiorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function inferiorSlider_Callback(hObject, eventdata, handles)
% hObject    handle to inferiorSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.preProcessingParams.sizeThresh(2) = round(get(hObject,'Value'));
set(handles.inferiorEdit,'string',sprintf('%i',uint32(handles.preProcessingParams.sizeThresh(2))));
guidata(hObject,handles);
figure(handles.figure1);


% --- Executes during object creation, after setting all properties.
function inferiorSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inferiorSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function superiorSlider_Callback(hObject, eventdata, handles)
% hObject    handle to superiorSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.preProcessingParams.sizeThresh(3) = round(get(hObject,'Value'));
set(handles.superiorEdit,'string',sprintf('%i',uint32(handles.preProcessingParams.sizeThresh(3))));
guidata(hObject,handles);
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function superiorSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to superiorSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function superiorEdit_Callback(hObject, eventdata, handles)
% hObject    handle to superiorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of superiorEdit as text
%        str2double(get(hObject,'String')) returns contents of superiorEdit as a double

handles.preProcessingParams.sizeThresh(3) = round(str2double(get(hObject,'String')));
guidata(hObject,handles);
figure(handles.figure1);


% --- Executes during object creation, after setting all properties.
function superiorEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to superiorEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function subBoundSlider_Callback(hObject, eventdata, handles)
% hObject    handle to subBoundSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.preProcessingParams.subBoundThresh = get(hObject,'Value');
set(handles.subBoundEdit,'string',sprintf('%f',handles.preProcessingParams.subBoundThresh));
guidata(hObject,handles);
figure(handles.figure1);



% --- Executes during object creation, after setting all properties.
function subBoundSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subBoundSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function subBoundEdit_Callback(hObject, eventdata, handles)
% hObject    handle to subBoundEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subBoundEdit as text
%        str2double(get(hObject,'String')) returns contents of subBoundEdit as a double
handles.preProcessingParams.subBoundThresh = str2double(get(hObject,'String'));
guidata(hObject,handles);
figure(handles.figure1);


% --- Executes during object creation, after setting all properties.
function subBoundEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subBoundEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function filterOmegaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to filterOmegaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.preProcessingParams.filterOmega = get(hObject,'Value');
set(handles.filterOmegaEdit,'string',sprintf('%f',handles.preProcessingParams.filterOmega));
guidata(hObject,handles);
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function filterOmegaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterOmegaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function filterOmegaEdit_Callback(hObject, eventdata, handles)
% hObject    handle to filterOmegaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filterOmegaEdit as text
%        str2double(get(hObject,'String')) returns contents of filterOmegaEdit as a double
handles.preProcessingParams.filterOmega = str2double(get(hObject,'String'));
guidata(hObject,handles);
figure(handles.figure1);


% --- Executes during object creation, after setting all properties.
function filterOmegaEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filterOmegaEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









% --- Executes on slider movement.
function reStractSlider_Callback(hObject, eventdata, handles)
% hObject    handle to reStractSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ratioShowLineImages = get(hObject,'Value');
set(handles.ratioOfExBounds,'string',sprintf('%f',handles.ratioShowLineImages));
guidata(hObject,handles);
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function reStractSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reStractSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ratioOfExBounds_Callback(hObject, eventdata, handles)
% hObject    handle to ratioOfExBounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ratioOfExBounds as text
%        str2double(get(hObject,'String')) returns contents of ratioOfExBounds as a double
handles.ratioShowLineImages = str2double(get(hObject,'String'));
figure(handles.figure1);

% --- Executes during object creation, after setting all properties.
function ratioOfExBounds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratioOfExBounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in onlySizeRadButton.
function onlySizeRadButton_Callback(hObject, eventdata, handles)
% hObject    handle to onlySizeRadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onlySizeRadButton


% --- Executes on button press in kappaRadButton.
function kappaRadButton_Callback(hObject, eventdata, handles)
% hObject    handle to kappaRadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kappaRadButton


% --- Executes on button press in manual3PointsButton.
function manual3PointsButton_Callback(hObject, eventdata, handles)
% hObject    handle to manual3PointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.busyFlag
    modelClass = getModelClassFromHandle(handles);
    if modelClass == 0;
        figure(handles.figure1);
        return;        
    end
    
    
    if isfield(handles,'imgIn')
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;        
        
        imshow(handles.imgIn);
        hold on;
        
        
        handles.busyFlag = 1;
        
        xImgW = [];
        gradxFW = [];
        gradyFW = [];
        xHatW = [];
        for kBound = 1:length(handles.subBounds)            
            gradxFW = [gradxFW,handles.subBounds(kBound).gradxF];
            gradyFW = [gradyFW,handles.subBounds(kBound).gradyF];
            xHatW = [xHatW,handles.subBounds(kBound).xHat];
            xImgW = [xImgW,handles.subBounds(kBound).xImg];
        end
        
        plot(xImgW(1,:),xImgW(2,:),'b.','MarkerSize',2);
        
        xSel = zeros(3,3);
        for k = 1:3
            [x_1,y_1] = ginput(1);
            iii=getCloser([x_1;y_1],xImgW(1:2,:));
            xSel(:,k) = xImgW(:,iii(1));
            
            plot(xSel(1,k),xSel(2,k),'r*');                        
        end
        
        T_center = [[1 0 -handles.u_0];[0 1 -handles.v_0];[0 0 1]];
        xHat = T_center*xSel;
        
        if modelClass==2                        
%             fH = gradientsToFocalHyperNoRVL([cos(gradThetaSel);sin(gradThetaSel)],xHat);      
%             [n, rVL, alpha, l1, l2, l3, r] = fitCurve(xHat,modelClass,fH);
%             calib = getCalibFromRVL(rVL,modelClass,fH);
        else
            [n, rVL, alpha, l1, l2, l3, r] = fitCurve(xHat,modelClass);
            calib = getCalibFromRVL(rVL,modelClass);                                                
        end

        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;                        
        imshow(handles.imgIn);
        hold on;   
        
        xHatParam=getParametricCurve(n,-pi:pi/360:pi,calib,modelClass);
        plot(xHatParam(1,:)+handles.u_0,xHatParam(2,:)+handles.v_0,'g-','LineWidth',5);  
        
        
        plot(xSel(1,:),xSel(2,:),'r*','MarkerSize',15,'LineWidth',4);            
        
        fprintf('u_0 = %f\n',handles.u_0);
        fprintf('v_0 = %f\n',handles.v_0);
        fprintf('rVL = %f\n',rVL);
        if modelClass == 2
            fprintf('fH = %f\n',calib(1));
            fprintf('xi = %f\n',cos(calib(2)));
        end

        handles.calib = calib;               
        guidata(hObject,handles);          
        
        
        handles.busyFlag = 0;
    end
    
end


% --- Executes on button press in singleBoundButton.
function singleBoundButton_Callback(hObject, eventdata, handles)
% hObject    handle to singleBoundButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global secondOrderFlag;

if ~handles.busyFlag
    modelClass = getModelClassFromHandle(handles);
    algorithm = getAlgorithmFromHandle(handles);

    if algorithm == 3
        secondOrderFlag = true;
    else
        secondOrderFlag = false;
    end

    if modelClass == 0;
        figure(handles.figure1);
        return;        
    end
    
    
    if isfield(handles,'imgIn')
        if handles.externFigure, figure(3); else axes(handles.mainFigure); end;        
        
        imshow(handles.imgIn);
        hold on;
        
        
        handles.busyFlag = 1;

        xImg = [];
        gradxF = [];
        gradyF = [];
        xHat = [];
        for kBound = 1:length(handles.subBounds)            
            gradxF = [gradxF,handles.subBounds(kBound).gradxF];
            gradyF = [gradyF,handles.subBounds(kBound).gradyF];
            xHat = [xHat,handles.subBounds(kBound).xHat];
            xImg = [xImg,handles.subBounds(kBound).xImg];
        end
        
        plot(xImg(1,:),xImg(2,:),'b.','MarkerSize',2);        
       
        nAttempts = handles.nAttemptsSingleBound;

        disp('Extracting lines ... Please Wait'); 
        shortPixelThreshold = handles.lineExtractionParams.normalizedShortPixelThreshold*handles.rows;


        if algorithm == 1
            [acceptedIndex, sizeAccepted, r_vlAcum] = uncalibratedLineExtraction3Points(xHat,xHat,nAttempts,shortPixelThreshold,modelClass,handles.rows);
        else
            gradThetaSel = atan2(gradyF,gradxF);    
            [acceptedIndex, sizeAccepted, r_vlAcum] = uncalibratedLineExtractionDeviceDistance(xHat,xHat,gradThetaSel,nAttempts,shortPixelThreshold,getSign(modelClass),modelClass,handles.rows); %% El ultimo campo permitiria el uso de fH si se conociese
            acceptedIndex = acceptedIndex(r_vlAcum>0);
            sizeAccepted = sizeAccepted(r_vlAcum>0);
            r_vlAcum = r_vlAcum(r_vlAcum>0);                                    
        end
      
        disp('Extracting lines done');                        
        

        binSel = imag(r_vlAcum)==0 & r_vlAcum>0;
        sizeAccepted = sizeAccepted(binSel);
        r_vlAcum = r_vlAcum(binSel);
        acceptedIndex = acceptedIndex(binSel);
                               
        [maxv,maxi] = max(sizeAccepted);
        
        fprintf('rVL = %f\n',r_vlAcum(maxi(1)));
        calib = getCalibFromRVL(r_vlAcum(maxi(1)),modelClass);
        handles.calib = calib;               
        guidata(hObject,handles);          
        handles.busyFlag = 0;
        set(handles.saveResultsButton,'Enable','on');
        set(handles.reExtractButton,'Enable','on');        
        reExtractButton_Callback(hObject, eventdata, handles);

    end    
    
end


% --- Executes on button press in manual2PointsButton.
function manual2PointsButton_Callback(hObject, eventdata, handles)
% hObject    handle to manual2PointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in optimCheckBox.
function optimCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to optimCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optimCheckBox
