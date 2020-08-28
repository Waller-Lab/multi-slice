function varargout = view3d(varargin)
%VIEW3D GUI for interactive viewing of 3D Volumes
%   VIEW3D is used view orthographic slices of 3D volumes
%    
%   Type in an expression that generates a 3D array
%   then press the Display button
%
%   3D expressions such as: rand(50,40,30) or
%   the name or a 3D array variable in the workspace
%
%   The 3 views (all except the lower right one)
%   display orthographic projections
%
%   Use the scroll bars to change the number of the slice viewed
%
%   Use the transpose, flipud, or fliplr to
%   transpose the view, flip it vertically, or horizontally
%
%   Use Update 3D to view the slices in a 3D view
%
%   Check auto to obtain an automatic update 
%   of the 3D view of the slices
%   (Note: this may affect performance)
%
%   Use cla to clear the 3D view this may improve performance
%
%   Change the span values to the volume''s physical dimensions
%   so the aspect ratio is displayed properly
%   (note: you can use relative values
%   for example use 1,3,2 instead of 0.5,1.5,1.0
%   then press Span
%
%   See also: SLICE, MONTAGE, ISOSURFACE
%
%   (c) Ghassan Hamarneh, ghamarneh@yahoo.com

%%% view3d([1 2 3])  --> only span vals

if ~isempty(varargin) & (all(size(varargin{1})==[3 1])  | all(size(varargin{1})==[1 3]))
    spanvar=varargin{1};
end

if nargin == 0 | exist('spanvar')% LAUNCH GUI    
    fig = openfig(mfilename,'reuse');        % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);    
    if nargout > 0
        varargout{1} = fig;
    end
    
    if exist('spanvar')
        set(handles.edit2,'string',num2str(spanvar(1))) 
        set(handles.edit3,'string',num2str(spanvar(2)))         
        set(handles.edit4,'string',num2str(spanvar(3))) 
    end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
end

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'view3d_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.view3d, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
myplot(handles,1);

% --------------------------------------------------------------------
function varargout = slider2_Callback(h, eventdata, handles, varargin)
myplot(handles,2);

% --------------------------------------------------------------------
function varargout = slider3_Callback(h, eventdata, handles, varargin)
myplot(handles,3);

% --------------------------------------------------------------------
function varargout = edit1_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)

a=get(handles.edit1,'String');
handles.vol=double(squeeze(evalin('base',a)));
if ndims(handles.vol)~=3,
    disp('not 3d')
    return
end

[handles.sx,handles.sy,handles.sz]=size(handles.vol);
set(handles.slider1,'min',1);
set(handles.slider2,'min',1);
set(handles.slider3,'min',1);
set(handles.slider1,'max',handles.sx);
set(handles.slider2,'max',handles.sy);
set(handles.slider3,'max',handles.sz);
set(handles.slider1,'value',round(handles.sx/2)+1);
set(handles.slider2,'value',round(handles.sy/2)+1);
set(handles.slider3,'value',round(handles.sz/2)+1);
cla(handles.axes4);axis([1 handles.sx 1 handles.sy 1 handles.sz]); axis vis3d


%axes(handles.axes1);imagesc(squeeze(handles.vol(1,:,:)));axis image;
%axes(handles.axes2);imagesc(squeeze(handles.vol(:,1,:)));axis image;
%axes(handles.axes3);imagesc(squeeze(handles.vol(:,:,1)));axis image;

set(gcf,'DoubleBuffer','on');

myplot(handles,[1 2 3])

%%% produced error in matlab 7.0
%if ~isfield(handles,'clrmnu')
%    handles.clrmnu=0;
%end
%if ~handles.clrmnu; 
%    colormenu; 
%    handles.clrmnu=1;
%end

guidata(h,handles);

% --------------------------------------------------------------------
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)
myplot(handles,1);
% --------------------------------------------------------------------
function varargout = checkbox2_Callback(h, eventdata, handles, varargin)
myplot(handles,2);
% --------------------------------------------------------------------
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)
myplot(handles,3);
% --------------------------------------------------------------------
function varargout = checkbox4_Callback(h, eventdata, handles, varargin)
myplot(handles,1);
% --------------------------------------------------------------------
function varargout = checkbox5_Callback(h, eventdata, handles, varargin)
myplot(handles,2);
% --------------------------------------------------------------------
function varargout = checkbox6_Callback(h, eventdata, handles, varargin)
myplot(handles,3);
% --------------------------------------------------------------------
function varargout = checkbox7_Callback(h, eventdata, handles, varargin)
myplot(handles,1);
% --------------------------------------------------------------------
function varargout = checkbox8_Callback(h, eventdata, handles, varargin)
myplot(handles,2);
% --------------------------------------------------------------------
function varargout = checkbox9_Callback(h, eventdata, handles, varargin)
myplot(handles,3);

% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
helpdlg({'3D Volume Orthoslice Viewer','(c) Ghassan Hamarneh 2002-2004'})

% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)
s1=round(get(handles.slider1,'value'));
s2=round(get(handles.slider2,'value'));
s3=round(get(handles.slider3,'value'));
axes(handles.axes4); hslc=slice(handles.vol,s1,s2,s3);%rotate3d on;
axis tight; set(hslc(1:3),'LineStyle','none');
xlabel 'x' ;ylabel 'y' ;zlabel 'z';
% --------------------------------------------------------------------
function varargout = checkbox10_Callback(h, eventdata, handles, varargin)
myplot(handles,[1 2 3]);

% --------------------------------------------------------------------
function myplot(handles,n)

sp1=1/str2num(get(handles.edit2,'string'));
sp2=1/str2num(get(handles.edit3,'string'));
sp3=1/str2num(get(handles.edit4,'string'));

vmx=max(handles.vol(:));
vmn=min(handles.vol(:));

s1=round(get(handles.slider1,'value'));
if any(n==1) % x -- y z
    I=squeeze(handles.vol(s1,:,:));
    if get(handles.checkbox1,'value'), I=I'; end
    if get(handles.checkbox4,'value'), I=flipud(I); end
    if get(handles.checkbox7,'value'), I=fliplr(I); end
    %    axes(handles.axes1);imagesc(I);
    axes(handles.axes1); 
    if get(handles.checkbox11,'value');imagesc(I,[vmn,vmx]);else imagesc(I);end
    colormap(gray)
    pbaspect('manual');pbaspect(handles.axes1,[sp2,sp3,sp1]);
    %    c=get(handles.axes1,'Children'); set(c(1),'cdata',I);
    set(handles.text1,'string',['x=',num2str(s1)])
    title('y-z')
end
s2=round(get(handles.slider2,'value'));
if any(n==2)% y -- x z
    I=squeeze(handles.vol(:,s2,:));
    if get(handles.checkbox2,'value'), I=I'; end
    if get(handles.checkbox5,'value'), I=flipud(I); end
    if get(handles.checkbox8,'value'), I=fliplr(I); end
    %    axes(handles.axes2);imagesc(I);
    axes(handles.axes2); 
    if get(handles.checkbox11,'value');imagesc(I,[vmn,vmx]);else imagesc(I);end
    colormap(gray)
    pbaspect('manual');pbaspect(handles.axes2,[sp1,sp3,sp2]);
    %    c=get(handles.axes2,'Children'); set(c(1),'cdata',I);
    set(handles.text2,'string',['y=',num2str(s2)])
    title('x-z')
end
s3=round(get(handles.slider3,'value'));
if any(n==3)% z -- x y
    I=squeeze(handles.vol(:,:,s3));
    if get(handles.checkbox3,'value'), I=I'; end
    if get(handles.checkbox6,'value'), I=flipud(I); end
    if get(handles.checkbox9,'value'), I=fliplr(I); end
    %    axes(handles.axes3); imagesc(I);
    axes(handles.axes3);
    if get(handles.checkbox11,'value');imagesc(I,[vmn,vmx]);else imagesc(I);end
    colormap(gray)
    pbaspect('manual');pbaspect(handles.axes3,[sp1,sp2,sp3]);
    %    c=get(handles.axes3,'Children'); set(c(1),'cdata',I);
    set(handles.text3,'string',['z=',num2str(s3)])
    title('x-y')
end
if get(handles.checkbox10,'value'), 
    pushbutton3_Callback([],[],handles,[]);
end

% --------------------------------------------------------------------
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
cla(handles.axes4);

% --------------------------------------------------------------------
function varargout = edit2_Callback(h, eventdata, handles, varargin)
myplot(handles,[1 2 3]);

% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)
myplot(handles,[1 2 3]);

% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)
myplot(handles,[1 2 3]);

% --------------------------------------------------------------------
function pushbutton5_Callback(hObject, eventdata, handles)
helpdlg({'- Type in an expression that generates a 3D array                     ',...
        '   then press the Load  button                                         ',...
        '- 3D expressions such as: rand(50,40,30) or                            ',...
        '    the name or a 3D array variable in the workspace                   ',...
        '- The 3 views (all except the lower right one)                         ',...
        '    display orthographic projections                                   ',...
        '- Use the scroll bars to change the number of the slice viewed         '...
        '- Use the transpose, flipud, or fliplr to                              '...
        '    transpose the view, flip it vertically, or horizontally            ',...
        '- Use update 3d to view the slices in a 3D view                        ',...
        '- Check auto to obtain an automatic update of the 3D view of the slices',...
        '    (Note: this may affect performance)                                ',...
        '- Use cla to clear the 3D view                                         ',...
        '    this may improve performance                                       ',...
        '- Change the span values to the volume''s physical dimensions          ',...
        '    so the aspect ratio is displayed properly                          ',...
        '    (note: you can use relative values                                 ',...
        '    for example use 1,3,2 instead of 0.5,1.5,1.0)                      ',...
        '    then press Apply                                                   '})

function pushbutton6_Callback(hObject, eventdata, handles)
if strcmp(questdlg('Exit View3D?','View3D','Yes','No','No'),'Yes')
    close(handles.view3d);
end

function pushbutton7_Callback(hObject, eventdata, handles)
myplot(handles,[1 2 3]);


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11


