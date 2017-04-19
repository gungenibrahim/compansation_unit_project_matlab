function varargout = busadd(varargin)
% BUSADD MATLAB code for busadd.fig
%      BUSADD, by itself, creates a new BUSADD or raises the existing
%      singleton*.
%
%      H = BUSADD returns the handle to a new BUSADD or the handle to
%      the existing singleton*.
%
%      BUSADD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BUSADD.M with the given input arguments.
%
%      BUSADD('Property','Value',...) creates a new BUSADD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before busadd_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to busadd_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help busadd

% Last Modified by GUIDE v2.5 19-May-2016 15:55:32

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @busadd_OpeningFcn, ...
                   'gui_OutputFcn',  @busadd_OutputFcn, ...
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


% --- Executes just before busadd is made visible.
function busadd_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to busadd (see VARARGIN)
global bus;
global busunit;
bus = zeros;
busunit=0;
global activehour;
global reactivehour;
activehour = 0;
reactivehour = 0;
global activepower;
global reactivepower;
activepower = zeros;
reactivepower = zeros;
global loadconnectedbus;
loadconnectedbus = zeros;
global loadcounter;
loadcounter = 0;
% Choose default command line output for busadd
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes busadd wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = busadd_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function SendingBus_Callback(hObject, eventdata, handles)
% hObject    handle to SendingBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SendingBus as text
%        str2double(get(hObject,'String')) returns contents of SendingBus as a double


% --- Executes during object creation, after setting all properties.
function SendingBus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SendingBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ReceivingBus_Callback(hObject, eventdata, handles)
% hObject    handle to ReceivingBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ReceivingBus as text
%        str2double(get(hObject,'String')) returns contents of ReceivingBus as a double


% --- Executes during object creation, after setting all properties.
function ReceivingBus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReceivingBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ParallelCables_Callback(hObject, eventdata, handles)
% hObject    handle to ParallelCables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ParallelCables as text
%        str2double(get(hObject,'String')) returns contents of ParallelCables as a double


% --- Executes during object creation, after setting all properties.
function ParallelCables_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParallelCables (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CableType.
function CableType_Callback(hObject, eventdata, handles)
% hObject    handle to CableType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CableType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CableType


% --- Executes during object creation, after setting all properties.
function CableType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CableType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Distance_Callback(hObject, eventdata, handles)
% hObject    handle to Distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Distance as text
%        str2double(get(hObject,'String')) returns contents of Distance as a double


% --- Executes during object creation, after setting all properties.
function Distance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Temperature_Callback(hObject, eventdata, handles)
% hObject    handle to Temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Temperature as text
%        str2double(get(hObject,'String')) returns contents of Temperature as a double


% --- Executes during object creation, after setting all properties.
function Temperature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function SendingBus_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to SendingBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in AddBus.
function AddBus_Callback(hObject, eventdata, handles)
% hObject    handle to AddBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global bus;
global busunit;
busunit = busunit +1;
 bus(busunit,1) = str2double(get(handles.SendingBus, 'string'));
 bus(busunit,2) = str2double(get(handles.ReceivingBus, 'string'));
 bus(busunit,3) = str2double(get(handles.ParallelCables, 'string'));
 bus(busunit,4)= str2double(get(handles.CableType,'string'));
 bus(busunit,5)= str2double(get(handles.Distance, 'string'));
 bus(busunit,6)= str2double(get(handles.Temperature, 'string'));
 





function LoadBus_Callback(hObject, eventdata, handles)
% hObject    handle to LoadBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LoadBus as text
%        str2double(get(hObject,'String')) returns contents of LoadBus as a double


% --- Executes during object creation, after setting all properties.
function LoadBus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadBus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ActivePower_Callback(hObject, eventdata, handles)
% hObject    handle to ActivePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ActivePower as text
%        str2double(get(hObject,'String')) returns contents of ActivePower as a double


% --- Executes during object creation, after setting all properties.
function ActivePower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ActivePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ReactivePower_Callback(hObject, eventdata, handles)
% hObject    handle to ReactivePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ReactivePower as text
%        str2double(get(hObject,'String')) returns contents of ReactivePower as a double


% --- Executes during object creation, after setting all properties.
function ReactivePower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReactivePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ActiveAdd.
function ActiveAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ActiveAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global activepower;
global loadcounter;
global activehour;

if activehour<23
    activepower(activehour+1,loadcounter) = str2double(get(handles.ActivePower, 'string'));
    set(handles.ActiveHour, 'string', activehour+1);
else
    activepower(activehour+1,loadcounter) = str2double(get(handles.ActivePower, 'string'));
    activehour=0;
    set(handles.ActiveHour, 'string', activehour);
end

activehour = activehour +1;

% --- Executes on button press in ReactiveAdd.
function ReactiveAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ReactiveAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global reactivepower;
global loadcounter;
global reactivehour;

if reactivehour<23
    reactivepower(reactivehour+1,loadcounter) = str2double(get(handles.ReactivePower, 'string'));
    set(handles.ReactiveHour, 'string', reactivehour+1);

else
    reactivepower(reactivehour+1,loadcounter) = str2double(get(handles.ReactivePower, 'string'));
    reactivehour=0;
    set(handles.ReactiveHour, 'string', reactivehour);
end

reactivehour = reactivehour +1;

% --- Executes on button press in DefineLoad.
function DefineLoad_Callback(hObject, eventdata, handles)
% hObject    handle to DefineLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global loadconnectedbus;
global loadcounter;
loadcounter = loadcounter +1;
loadconnectedbus(loadcounter,1) = str2double(get(handles.LoadBus, 'string'));



% --- Executes on button press in Calculate.
function Calculate_Callback(hObject, eventdata, handles)
% hObject    handle to Calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global loadconnectedbus;
global bus;
global activepower;
global reactivepower;
global hour;
global Drop1;
global Drop2;
global Factor1;
global Factor2;
global Loading1;
global Loading2;
global busunit;
busunit = busunit -1;
global Sub1;
global Sub2;
global XSUB;
global XLOAD;
hour = 1:24;
[Drop1,Drop2,Factor1,Factor2,Loading1,Loading2,Sub1,Sub2,XSUB,XLOAD] = LetMeCalculate(bus,loadconnectedbus,activepower,reactivepower);



function [Voltagedrop_case1,Voltagedrop_case2,pfsubcase1,pfsubcase2_comp,itotal_busload_case1,itotal_busload,Isub_case1,Isub_case2,Xsub,xload]=LetMeCalculate(buses,loadsbus,PkW,QkVAr)
Vt = 400; %line to line voltage
freq = 50;
f=50;
Xcable1 = 0.04;
Xcable2 = 0.07;
Acable1 = 240;
Acable2 = 70;
Xcopper = 56;
%%busall
% column 1---->sending bus 
% column 2----> receiving bus
% column 3-----> cable type
% column 4-----> parallel lines
% column 5----> distance
% column 6----> temprature


%% case1

len = size(PkW);

for i = 1:len(2)
    for j=1:24
    if PkW(j,i) == 0 && QkVAr(j,i) == 0
        pf(j,i) = 1;
    else
    pf(j,i) = PkW(j,i) / sqrt((PkW(j,i))^2+ (QkVAr(j,i))^2);
    end
    end
end

for i = 1:1:len(2)

    for j = 1:24
    if mean(QkVAr(j,i))>0
        Qmax(j,i) = sqrt(((PkW(j,i)/0.98)^2)-(PkW(j,i)^2));
    elseif mean(QkVAr(i))<0
        Qmax(j,i) = -sqrt(((PkW(j,i)/0.95)^2)-(PkW(j,i)^2));
    else 
        Qmax(j,i) = 0;
    end
    
    end
end

for i = 1:1:len(2)
    
    for j = 1:24
    
    QCmin(j,i) = QkVAr(j,i)- Qmax(j,i);
    
    end
end

Qtot = zeros(1,3);


% pek gerekli bi parta benzemiyor.

for i = 1:1:len(2)
   a=0;
    for j = 1:24
    
    if QkVAr(j,i)>0
        if pf(j,i)<0.98
            Qtot(i) = Qtot(i) + QCmin(j,i);
            a = a+1;
            A(a,i) = QCmin(j,i);
        end
    elseif QkVAr(j,i)<0
        if pf(j,i)<0.95
            Qtot(i) = Qtot(i) + QCmin(j,i);
            a = a+1;
            A(a,i) = QCmin(j,i);
        end 
    end    
    end
end
maximum = zeros(1,len(2));

minimum = zeros(1,len(2));

for j=1:len(2)
    
    B = QCmin(:,j);
    if mean(B)>0
    if mean(B)~= 0
        maximum(1,j) = max(B(B>0));
        minimum(1,j) = min(B(B>0));
    else
        maximum(1,j) = 0;
        minimum(1,j) = 0;
    
    end
    else
        if mean(B)~= 0
            maximum(1,j) = max(B(B<0));
            minimum(1,j) = min(B(B<0));
        else
            maximum(1,j) = 0;
            minimum(1,j) = 0;
        end
    end
    
end




% for loads find switch capacitance or inductance X means both component

for j = 1:len(2)
    if mean(QkVAr(:,j))>=0 
        Xmaxload(j) = (maximum(j)+0.15*maximum(j))*10^3/(400^2*(2*pi*f));
        Xminload(j) = (minimum(j)+0.15*minimum(j))*10^3/(400^2*(2*pi*f));
    else
        Xmaxload(j) = -400^2/((2*pi*f)*(minimum(j)+0.15*minimum(j))*10^3);
        Xminload(j) = -400^2/((2*pi*f)*(maximum(j)+0.15*maximum(j))*10^3);
    end

end

% code in below find capacitor values
for j = 1:len(2)
    xloada = linspace(Xminload(j), Xmaxload(j),3);
    xload(:,j)= xloada;
    Qcompa =  linspace(1.15*minimum(j),1.15*maximum(j),3);
    Qcomp(:,j) = Qcompa;
end
%minimum compansete values increased by %20 to get better values
QkVAr1 = QkVAr; 

for j = 1:len(2)
    for i = 1:24

        %load1
            if QkVAr(i,j)>= Qcomp(1,j) && QkVAr(i,j)< Qcomp(2,j)
                QkVAr1(i,j) = QkVAr(i,j) - Qcomp(1,j);
                

            elseif QkVAr(i,j)>= Qcomp(2,j) && QkVAr(i,j)< Qcomp(3,j)
                QkVAr1(i,j) = QkVAr(i,j) - Qcomp(2,j);
                
            elseif QkVAr(i,j)>= Qcomp(3,j) 
                QkVAr1(i,j) = QkVAr(i,j) - Qcomp(3,j);

            end
    end
end


for i = 1:len(2)
    for j=1:24
    if PkW(j,i) == 0 && QkVAr1(j,i) == 0
        pf1(j,i) = 1;
    else
    pf1(j,i) = PkW(j,i) / sqrt((PkW(j,i))^2+ (QkVAr1(j,i))^2);
    end
    end
end

%% case2 


for j = 1:len(2)
    for i = 1:24
        S_case2_loads(i,j) = sqrt(PkW(i,j)^2+QkVAr(i,j)^2); % kVars
        S_case1_load_comp(i,j) = sqrt(PkW(i,j)^2+QkVAr1(i,j)^2);%kVars
        
        iload_active(i,j) = (PkW(i,j)/(sqrt(3)*Vt))*1000;  %Ampers
        iload_reactive(i,j) = QkVAr(i,j)/(sqrt(3)*Vt)*1000; %ampers
        
        
        iload_reactive_case1(i,j) = QkVAr1(i,j)/(sqrt(3)*Vt)*1000; %amper
            
        
        iload_total_uncomp(i,j)= S_case2_loads(i,j)/(sqrt(3)*Vt)*1000; %amps
        iload_total_comp(i,j)= S_case1_load_comp(i,j)/(sqrt(3)*Vt)*1000; %amps
    end 
end

%%bus
% column 1---->sending bus 
% column 2----> receiving bus
% column 3-----> cable type
% column 4-----> parallel lines
% column 5----> distance
% column 6----> temprature



% below code gives us each bus current active and reactive ones.
t = size(buses);
num_of_bus = max(buses(:,2));
ibusload_active = zeros(24,num_of_bus);
ibusload_reactive = zeros(24,num_of_bus);
ibload_reactive_case1= zeros(24,num_of_bus);
for f = 1: num_of_bus
    j = num_of_bus +1 - f;
    for i = 1:len(2)
        if loadsbus(i,1) == j
            for m = 1:24
            ibusload_active(m,j) = ibusload_active(m,j) + iload_active(m,i);
            ibusload_reactive(m,j) = ibusload_reactive(m,j) + iload_reactive(m,i);
            ibload_reactive_case1(m,j) = ibload_reactive_case1(m,j) +iload_reactive_case1(m,i);
            end
        end
    end
    for b = 1:num_of_bus
        if buses(b,1) == buses(j,2)
            w = buses(b,1);
            q = buses(b,2);
            
            for m = 1:24
                ibusload_active(m,w) =  ibusload_active(m,w) + ibusload_active(m,q) ;
                ibusload_reactive(m,w) = ibusload_reactive(m,w) + ibusload_reactive(m,q);
                ibload_reactive_case1(m,j) = ibload_reactive_case1(m,j) +ibload_reactive_case1(m,q);
            end
        end
    end
end

%below code gives each buses resistance and reactance 
row2=1;
for f = 1:num_of_bus
    row = (buses(:,2)==f);
    for i=1:num_of_bus
        if(row(i) == 1)
        row2 = i;
        end
    end
        if buses(row2,4) == 1
            Acable = Acable1;
            Xcable = Xcable1;
        elseif buses(row2,4) == 2
            Acable = Acable2;
            Xcable = Xcable2;
        end
        
        Rbus(f,1) = buses(row2,5)/(Xcopper*Acable*buses(row2,3)); %ohm
        %parallel line and cable type change
        Xbus(f,1) = buses(row2,5)*(10^(-3))*Xcable/buses(row2,3); %ohm
end
    %cable type will define
    
    Pbusline=zeros(24,num_of_bus);
    
for i = 1:num_of_bus
    for j = 1:24
        itotal_busload(j,i) = sqrt(ibusload_active(j,i)^2 + ibusload_reactive(j,i)^2);
        
        itotal_busload_case1(j,i) = sqrt(ibusload_active(j,i)^2 + ibload_reactive_case1(j,i)^2);
        
        
        Pbusline(j,i) = (itotal_busload(j,i)^2)*Rbus(i,1)/1000; %kwatt
        Qbusline(j,i) = (itotal_busload(j,i)^2)*Xbus(i,1)/1000; %kvar
        Qbusline_case1(j,i) =  (itotal_busload_case1(j,i)^2)*Xbus(i,1)/1000; %kvar
    end
end


for i = 1:24
    Ptotalline(i) = sum(Pbusline(i,:));
    Qtotalline(i) = sum(Qbusline(i,:));
    Ptotalload(i) = sum(PkW(i,:));
    Qtotalload(i) = sum (QkVAr(i,:));
    Qtotalload_case1(i) = sum(QkVAr1(i,:));
    Qtotalline_case1(i) = sum(Qbusline_case1(i,:));
    
    
    Psub(i) = Ptotalline(i)+Ptotalload(i);
    Qsub(i) = Qtotalline(i)+Qtotalload(i);
    Qsub_case1(i)= Qtotalline_case1(i)+Qtotalload_case1(i);
    
    Ssub_case1(i) = sqrt(Psub(i)^2+Qsub_case1(i)^2);
    Isub_case1(i) = 1000*Ssub_case1(i)/(sqrt(3)*Vt);
    Ssub_case2(i) = sqrt(Psub(i)^2+Qsub(i)^2);
    Isub_case2(i) = 1000*Ssub_case2(i)/(sqrt(3)*Vt);
    if Psub(i) == 0 && Qsub(i) == 0
        pfsubcase2(i) = 1;
    else
        pfsubcase2(i) = Psub(i) / sqrt((Psub(i))^2+ (Qsub(i))^2);
    end
    
    if Psub(i) == 0 && Qsub_case1(i) == 0
        pfsubcase1(i) = 1;
    else
        pfsubcase1(i) = Psub(i) / sqrt((Psub(i))^2+ (Qsub_case1(i))^2);
    end
    
end


for i = 1:24
    
    if mean(Qsub)>0
        Qmaxsub(i) = sqrt((Psub(i)/0.98).^2-Psub(i).^2);
    elseif mean(Qsub)<0
        Qmaxsub(i) = -sqrt((Psub(i)/0.95).^2-Psub(i).^2);
    else
        Qmaxsub(i) = 0;
    end
    
    Qcomp_min_sub(i) = Qsub(i) - Qmaxsub(i);
    
end


Qtotsub = 0;
for i = 1:24
    if Qsub(i)>0
        if pfsubcase2(i)<0.98
            Qtotsub = Qtotsub + Qcomp_min_sub(i);
            a = a+1;
            Asub(a) = Qcomp_min_sub(i);
        end
    elseif Qsub(i)<0
        if pfsubcase2(i)<0.95
            Qtotsub = Qtotsub + Qcomp_min_sub(i);
            a = a+1;
            Asub(a) = Qcomp_min_sub(i);
        end
    else
        Asub = zeros;
    end
end

Z = Qcomp_min_sub;

if mean(Qsub)>0    
    min_sub = min(Z(Z>0));
    max_sub = max(Z(Z>0));
else 
    min_sub = min(Z(Z<0));
    max_sub = max(Z(Z<0));
end

if max_sub>=0 && min_sub>=0
    Xsubmax   = (max_sub+0.15*max_sub)*10^3/(400^2*(2*pi*f));
    Xsubmin = (min_sub+0.15*min_sub)*10^3/(400^2*(2*pi*f));
elseif max_sub<=0 && min_sub<=0 
    Xsubmax = 400^2/((2*pi*f)*(min_sub+0.15*min_sub)*10^3);
    Xsubmin = 400^2/((2*pi*f)*(max_sub+0.15*max_sub)*10^3);
end

Xsub = linspace(Xsubmin,Xsubmax,5);
Qcom_sub = linspace(1.15*min_sub,1.15*max_sub,5);
Qsub_comp=zeros(24,1);
pfsubcase2_comp = pfsubcase2;
for i = 1:24
      
        if Qsub(i)>= Qcom_sub(1) && Qsub(i)< Qcom_sub(2)
            Qsub_comp(i) = Qsub(i) - Qcom_sub(1)- 0.1* Qcom_sub(1);
        
        elseif Qsub(i)>= Qcom_sub(2) && Qsub(i)< Qcom_sub(3)
            Qsub_comp(i) = Qsub(i) - Qcom_sub(2)- 0.1* Qcom_sub(2);
        
        elseif Qsub(i)>= Qcom_sub(3) && Qsub(i)< Qcom_sub(4)
            Qsub_comp(i) = Qsub(i) - Qcom_sub(3)- 0.1* Qcom_sub(3);  
        
        elseif Qsub(i)>= Qcom_sub(4) && Qsub(i)< Qcom_sub(5)
            Qsub_comp(i) = Qsub(i) - Qcom_sub(4)- 0.1* Qcom_sub(4);
        
        elseif Qsub(i)>= Qcom_sub(5) 
            Qsub_comp(i) = Qsub(i) - Qcom_sub(5)- 0.1* Qcom_sub(5);
      
        end
    
    
   
        pfsubcase2_comp(i) = Psub(i) / sqrt((Psub(i))^2+ (Qsub_comp(i))^2);
     
end

for j = 1:len(2)
    for i=1:24
        Qbus_fromloads_case1(i,j) = sqrt(3)*ibload_reactive_case1(i,j)*Vt; %var
        Pbus_fromloads(i,j) = sqrt(3)*ibusload_active(i,j)*Vt; %watt
        
        % for case1
        Voltagedrop_case1(i,j) = (Pbus_fromloads(i,j).*Rbus(j) + Qbus_fromloads_case1(i,j).*Xbus(j))/400; % (%V)


    %     for case 2

      Qbus_fromloads_case2(i,j) = sqrt(3)*ibusload_reactive(i,j)*Vt;
      
        
       
     Voltagedrop_case2(i,j) = (Pbus_fromloads(i,j).*Rbus(j) + Qbus_fromloads_case2(i,j).*Xbus(j))/400; % (%V)
    end
end


% --- Executes on button press in BusProfile.
function BusProfile_Callback(hObject, eventdata, handles)
% hObject    handle to BusProfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global bus;
    global busunit;
    bus = xlsread('BusProfile.xlsx',-1);
    a=size(bus);
    busunit=a(1);

% --- Executes on button press in LoadConnection.
function LoadConnection_Callback(hObject, eventdata, handles)
% hObject    handle to LoadConnection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global loadconnectedbus;
    global loadcounter;
    loadconnectedbus = xlsread('LoadConnection.xlsx',-1);
    loadcounter = max(loadconnectedbus(:,1));

% --- Executes on button press in LoadActivePower.
function LoadActivePower_Callback(hObject, eventdata, handles)
% hObject    handle to LoadActivePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global activepower;
    activepower = xlsread('LoadProfile.xlsx',-1);

% --- Executes on button press in LoadReactivePower.
function LoadReactivePower_Callback(hObject, eventdata, handles)
% hObject    handle to LoadReactivePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global reactivepower;
    reactivepower = xlsread('LoadProfile.xlsx',-1);



function Number_Callback(hObject, eventdata, handles)
% hObject    handle to Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Number as text
%        str2double(get(hObject,'String')) returns contents of Number as a double


% --- Executes during object creation, after setting all properties.
function Number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowButton.
function ShowButton_Callback(hObject, eventdata, handles)
% hObject    handle to ShowButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Drop1;
global Drop2;
global Factor1;
global Factor2;
global Loading1;
global Loading2;
global Sub1;
global Sub2;
global hour;

a = str2double(get(handles.Number, 'string'));
if(a==0)
    axes(handles.VoltageDropAxes);
plot(hour,0,hour,0);
legend('Case 1','Case 2');
xlabel('Hour');
ylabel('Regulation');
title('Regulation vs Hour Graph');
    axes(handles.PowerFactorAxes);
plot(hour,Factor1,hour,Factor2);
xlabel('Hour');
ylabel('Power Factor');
title('Power Factor vs Hour Graph');
legend('Case 1','Case 2');
    axes(handles.LineLoadingAxes);
plot(hour,Sub1(:,1),hour,Sub2(:,1));
xlabel('Hour');
ylabel('Line Loadings');
title('Line Loadings vs Hour Graph');
legend('Case 1','Case 2');
else
        axes(handles.VoltageDropAxes);
plot(hour,Drop1(:,a),hour,Drop2(:,a));
legend('Case 1','Case 2');
xlabel('Hour');
ylabel('Regulation');
title('Regulation vs Hour Graph');
    axes(handles.PowerFactorAxes);
plot(hour,Factor1,hour,Factor2);
xlabel('Hour');
ylabel('Power Factor');
title('Power Factor vs Hour Graph');
legend('Case 1','Case 2');
    axes(handles.LineLoadingAxes);
plot(hour,Loading1(:,a),hour,Loading2(:,a));
xlabel('Hour');
ylabel('Line Loadings');
title('Line Loadings vs Hour Graph');
legend('Case 1','Case 2');
end



function LoadNo_Callback(hObject, eventdata, handles)
% hObject    handle to LoadNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LoadNo as text
%        str2double(get(hObject,'String')) returns contents of LoadNo as a double


% --- Executes during object creation, after setting all properties.
function LoadNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ShowCapacitance.
function ShowCapacitance_Callback(hObject, eventdata, handles)
% hObject    handle to ShowCapacitance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global XSUB;
global XLOAD;
a = str2double(get(handles.LoadNo, 'string'));
if(a==0)
    set(handles.Capacitances, 'string', XSUB);
else
set(handles.Capacitances, 'string', XLOAD(:,a));
end
