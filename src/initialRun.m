%% Initialize
set(0,'DefaultFigureWindowStyle','docked'); 

%% --- PORTABLE SETUP (GEFs REPO) ---
% Getting the location of this script
[currentDir, ~, ~] = fileparts(mfilename('fullpath'));
% Getting MATLAB wrappers path
wrapperPath = fullfile(currentDir, 'GEFs_core', 'matlab_wrappers');
if exist(wrapperPath, 'dir')
    addpath(wrapperPath);
else
    error('Cannot find wrappers at: %s. Check your folder structure!', wrapperPath);
end
% Linking to the Python Virtual Environment (.venv)
% It is located 2 levels up from this script (GEFs/src/GEFs_core)
projectRoot = fileparts(fileparts(currentDir));
venvPath = fullfile(projectRoot, '.venv', 'Scripts', 'python.exe');

if exist(venvPath, 'file')
    if isempty(pyenv().Executable) || ~contains(pyenv().Executable, '.venv')
        pyenv('Version', venvPath, 'ExecutionMode', 'OutOfProcess');
    end
end
% 4. Adding the GEFs_core folder to Python search path
pythonSrc = fullfile(currentDir, 'GEFs_core');
if count(py.sys.path, pythonSrc) == 0
    insert(py.sys.path, int32(0), pythonSrc);
end

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
