function sliderDisplayImVC2(data, figMod, circ)
% This function creates plots for data, a multidimensional array (3D or
% more).
% Sliders are created to navigate through data.
% 
% For example, let's assume data is 3D, containing images (in the first two 
% dimensions) for different values of defocus (third dimension). 
% sliderDisplay will display these images for a certain defocus value that
% can be changed with the slider.
% 
% For first time use, sliderDisplay() with no parameter will automatically
% generate data to plot.
% figMod is a string of figure commands (such as caxis, colorbar, etc.)
% 
% To adapt sliderDisplay to other situations, there are only two sections
% to modify: "Values to change" in the beginning and "Here is what you want
% to display" in the very end.

%% Data generation (only for testing)

if nargin == 0
    data = zeros(100,100,7,6);
    for i = 1:size(data,3)
        for j = 1:size(data,4)
            data(:,:,i,j) = rand(100);
        end
    end
    circ=zeros(size(data,3),6);
    figMod = '';
elseif nargin == 1
    circ=zeros(size(data,3),6);
    figMod = '';
elseif nargin == 2
    circ=zeros(size(data,3),6);
end

%% Values to change

if numel(size(circ))<3
    sliderNumber = size(size(circ),2) - 1;
else
    sliderNumber = size(size(circ),2) - 2;
end
sliderDimension1 = 1;  % Dimension corresponding to Slider 1 (3 for the example above)
sliderDimension2 = 4;  % Same for Slider 2, won't be used if nSlider < 2
sliderDimension3 = 5;  % Same for Slider 3, won't be used if nSlider < 3

screenWidth = 1920;
sliderLeftPosition = 0;  % Change to round(5*screenWidth/12) if it covers the image
sliderLength = round(screenWidth/6);

% For display
sliderName1 = 'Time';  % Description of sliderDimension1 (defocus for the example above)
sliderName2 = '';
sliderName3 = '';
sliderStep1 = 1;  % Step for each dimension (defocus distance between two images for the example above)
sliderStep2 = 1;
sliderStep3 = 1;

%% Other initializations

sliderTotal1 = size(circ,sliderDimension1);
sliderTotal2 = size(circ,sliderDimension2);
sliderTotal3 = size(circ,sliderDimension3);
if strcmp(sliderName1,'')
    sliderName1 = 'Slider 1 position';
end
if strcmp(sliderName2,'')
    sliderName2 = 'Slider 2 position';
end
if strcmp(sliderName3,'')
    sliderName3 = 'Slider 3 position';
end

%% Display initialization

currentPosition1 = 1;
currentPosition2 = 1;
currentPosition3 = 1;
lastPosition1 = currentPosition1;
lastPosition2 = currentPosition2;
lastPosition3 = currentPosition3;

figure();

%% Slider and text creation

if sliderNumber >= 1
    handleSlider1 = uicontrol('Style', 'slider',...
        'Min',1,'Max',sliderTotal1,'Value',currentPosition1,...
        'sliderStep', [min(1/(sliderTotal1-1),0.1) 0.2], ...
        'Position', [sliderLeftPosition 25 sliderLength 25]);
    addlistener(handleSlider1,'ContinuousValueChange',@callback1);
    
    % Text
    handleText1 = uicontrol('Style','text',...
        'Position',[sliderLeftPosition 0 sliderLength 25],...
        'String',[sliderName1 ': ' num2str(currentPosition1*sliderStep1)]);
end
if sliderNumber >= 2
    handleSlider2 = uicontrol('Style', 'slider',...
        'Min',1,'Max',sliderTotal2,'Value',currentPosition2,...
        'sliderStep', [min(1/(sliderTotal2-1),0.1) 0.2], ...
        'Position', [sliderLeftPosition 75 sliderLength 25]);
    addlistener(handleSlider2,'ContinuousValueChange',@callback2);
    
    % Text
    handleText2 = uicontrol('Style','text',...
        'Position',[sliderLeftPosition 50 sliderLength 25],...
        'String',[sliderName2 ': ' num2str(currentPosition2*sliderStep2)]);
end
if sliderNumber >= 3
    handleSlider3 = uicontrol('Style', 'slider',...
        'Min',1,'Max',sliderTotal3,'Value',currentPosition3,...
        'sliderStep', [min(1/(sliderTotal3-1),0.1) 0.2], ...
        'Position', [sliderLeftPosition 125 sliderLength 25]);
    addlistener(handleSlider3,'ContinuousValueChange',@callback3);
    
    % Text
    handleText3 = uicontrol('Style','text',...
        'Position',[sliderLeftPosition 100 sliderLength 25],...
        'String',[sliderName3 ': ' num2str(currentPosition3*sliderStep3)]);
end

%% Display

display();

%% Functions called when slider is moved

    function callback1(varargin)
        currentPosition1 = round(get(handleSlider1,'Value'));
        set(handleSlider1,'Value',currentPosition1);
        
        if currentPosition1 ~= lastPosition1
            currentPosition2 = lastPosition2;
            currentPosition3 = lastPosition3;
            
            display();
            
            lastPosition1 = currentPosition1;
        end
    end
    function callback2(varargin)
        currentPosition2 = round(get(handleSlider2,'Value'));
        set(handleSlider2,'Value',currentPosition2);
        
        if currentPosition2 ~= lastPosition2
            currentPosition1 = lastPosition1;
            currentPosition3 = lastPosition3;
            
            display();
            
            lastPosition2 = currentPosition2;
        end
    end
    function callback3(varargin)
        currentPosition3 = round(get(handleSlider3,'Value'));
        set(handleSlider3,'Value',currentPosition3);
        
        if currentPosition3 ~= lastPosition3
            currentPosition1 = lastPosition1;
            currentPosition2 = lastPosition2;
            
            display();
            
            lastPosition3 = currentPosition3;
        end
    end

%% Display

    function display()
        % Change text
        if sliderNumber >= 1
            delete(handleText1);
            handleText1 = uicontrol('Style','text',...
                'Position',[sliderLeftPosition 0 sliderLength 25],...
                'String',[sliderName1 ': ' num2str(currentPosition1*sliderStep1)]);
        end
        if sliderNumber >= 2
            delete(handleText2);
            handleText2 = uicontrol('Style','text',...
                'Position',[sliderLeftPosition 50 sliderLength 25],...
                'String',[sliderName2 ': ' num2str(currentPosition2*sliderStep2)]);
        end
        if sliderNumber >= 3
            delete(handleText3);
            handleText3 = uicontrol('Style','text',...
                'Position',[sliderLeftPosition 100 sliderLength 25],...
                'String',[sliderName3 ': ' num2str(currentPosition3*sliderStep3)]);
        end
        
        cla();  % Clear axes
        
        %% Here is what you want to display
        % currentPosition1 represents the current value for
        % sliderDimension1. For the example above, to display images at
        % a given defocus distance, use:
        %     imshow(data(:,:,currentPosition1));
        
        imagesc(squeeze(data(:,:,currentPosition1)));
        %imagesc(squeeze(data(:,:,currentPosition1,currentPosition2,currentPosition3)));
        col={'r','g','m','b','c','w','y','k','r','g','m','b','c'};
%         axis(451+100*[-1 1 -1 1]); 
        axis image; 
        hold on;
        for kk=1:size(circ,3)
            viscircles(circ(currentPosition1,1:2,kk,currentPosition2,currentPosition3),circ(currentPosition1,3,kk,currentPosition2,currentPosition3),'EdgeColor',col{kk},'LineWidth',0.25);
        end
        
        for ii=1:length(figMod)
            eval(figMod{ii});
        end
        
        
    end
end