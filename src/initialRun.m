%% Initialize
set(0,'DefaultFigureWindowStyle','docked'); 

% Set a consistent, readable font size (12 is usually the "sweet spot")
fontSize = 12; 
set(0, 'DefaultTextFontSize', fontSize);
set(0, 'DefaultAxesFontSize', fontSize);

% Line and Marker settings
set(0, 'DefaultLineLineWidth', 1.5); % 2 can be a bit thick for complex plots
set(0, 'DefaultLineMarkerSize', 6);

% Grid and Axis look
set(0, 'DefaultAxesTickLength', [0.02, 0.02]);
set(0, 'DefaultAxesLineWidth', 1);