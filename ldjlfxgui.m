function varargout = ldjlfxgui(varargin)
% LDJLFXGUI MATLAB code for ldjlfxgui.fig
%      LDJLFXGUI, by itself, creates a new LDJLFXGUI or raises the existing
%      singleton*.
%
%      H = LDJLFXGUI returns the handle to a new LDJLFXGUI or the handle to
%      the existing singleton*.
%
%      LDJLFXGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LDJLFXGUI.M with the given input arguments.
%
%      LDJLFXGUI('Property','Value',...) creates a new LDJLFXGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ldjlfxgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ldjlfxgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ldjlfxgui

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ldjlfxgui_OpeningFcn, ...
                   'gui_OutputFcn',  @ldjlfxgui_OutputFcn, ...
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


% --- Executes just before ldjlfxgui is made visible.
function ldjlfxgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ldjlfxgui (see VARARGIN)

% Choose default command line output for ldjlfxgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ldjlfxgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ldjlfxgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N=1000;
tmin=0;
tmax=10;                           %PRI取值范围设定
K=201;   %PRI箱个数，K越大，PRI变换图越精确
N1=464;
N2=328;
N3=207;
p1=1;
p2=sqrt(2);
p3=sqrt(5);
b=(tmax-tmin)/K;                    %每个箱子的宽度
for k=1:K,
    d1(k)=0;
    d2(k)=0;
    d3(k)=0;  
    x(k)=(k-1/2)*b+tmin;           %每个箱子的中心坐标
    for n1=1:N1-1,
        for m1=0:n1-1,
                tao1=(n1-m1)*p1;
                          %tao=tn-tm
            if (tao1>=x(k)-b/2)&(tao1<=x(k)+b/2)
               d1(k)=d1(k)+exp(2*pi*i*(tmin+n1*p1)/tao1);
            else
                d1(k)=d1(k);
            end
        end
    end
    for n2=1:N2-1,
        for m2=0:n2-1,
     
            tao2=(n2-m2)*p2;
                   
            if (tao2>=x(k)-b/2)&(tao2<=x(k)+b/2)
               d2(k)=d2(k)+exp(2*pi*i*(tmin+n2*p2)/tao2);                                   
            else
                d2(k)=d2(k);
            end
        end
    end
     for n3=1:N3-1,
        for m3=0:n3-1,
             tao3=(n3-m3)*p3;
            if (tao3>=x(k)-b/2)&(tao3<=x(k)+b/2)
               d3(k)=d3(k)+exp(2*pi*i*(tmin+n3*p3)/tao3);
            else
                d3(k)=d3(k);
            end
        end
     end
    y1(k)=abs(d1(k));                    
    y2(k)=abs(d2(k));
    y3(k)=abs(d3(k));         %Dk取绝对值
      
end
 axes(handles.axes1);
 plot(x,y1,x,y2,x,y3);                 %画图
 title('雷达PRI脉冲序列(峰值为分选出的雷达信号)')
 xlabel('脉冲间隔')
 ylabel('PRI变换值')
 legend('PRI为1的脉冲序列','PRI为√2的脉冲序列','PRI为√5的脉冲序列')


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla;
mgd=str2double(get(handles.edit1,'string'));
load xyz;
[clusters,clusterInds,clusterBounds] = clusterData(xyz,[0.25 0.4 mgd]);
colors = jet(numel(clusters));
inds = randperm(numel(clusters));
colors = colors(inds,:);
axes(handles.axes1);
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),'b.')
xlabel('x');ylabel('y');zlabel('z');
title('雷达网格聚类分选(更精细的分选)');
hold on
for ii = 1:numel(clusters)
plot3(clusters{ii}(:,1),clusters{ii}(:,2),clusters{ii}(:,3),...
    'linestyle','none','linewidth',2,'marker','o','color',colors(ii,:));
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla;
load kmeans;
X = meas(:,1:2); 
[n,p] = size(X);
rng(3); 
k = 3;
Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};
SCtext = {'true','false'};
nSC = numel(SharedCovariance);
d = 500;
x1 = linspace(min(X(:,1)) - 2,max(X(:,1)) + 2,d);
x2 = linspace(min(X(:,2)) - 2,max(X(:,2)) + 2,d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000); 

cluster0 = {[ones(n-8,1); [2; 2; 2; 2]; [3; 3; 3; 3]];...
            randsample(1:k,n,true); randsample(1:k,n,true); 'plus'};
converged = nan(4,1);
axes(handles.axes1);
j = 2;

gmfit = fitgmdist(X,k);
clusterX = cluster(gmfit,X);
mahalDist = mahal(gmfit,X0);
h1 = gscatter(X(:,1),X(:,2),clusterX);

title('雷达信号K均值聚类分选')

sum(converged);
