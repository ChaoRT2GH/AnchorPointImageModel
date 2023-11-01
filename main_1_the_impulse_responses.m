%%**************************************************
% This is the mainfile of the APIM room imuplse response generation
% approach, which is an extension of the traditional image model method. 
% In this project, we consider the source is direcitonal, and the sensor
% are omni-directional. 
% ** one can simply run the file to get an example of the APIM approach
%%**************************************************
clc;clear all; close all;
%% setting the simulation coditions
% param. sound speed, sampling rate, and length of the impulse response
c = 340;     % sound speed, in m/s
fs = 16e3;    % sampling rate, in Hz
Lh = 1024*2;   % length of IR,

% param. size parameters of the room, in meter
vec_L = [4;4;3];

% param. reflection coefficients
mat_beta = [0.96,0.8;...
    0.96,0.9;
    0.5,0.5];

% param. position of souce and microphone, in meter
vec_rs = [3;3;1];     % position of the souce
vec_as = [3.1;3.1;1]; % position of the anchor piont
vec_rm = [1.5;1.5;1]; % position of the microphone, the microphone here is
                      % omni-directional;therefore, no-need to provide an anchor point

%% the impulse response with the traditional image model method
Qmax = -1;  % degen. to the traditional image model method
[hVec] = cal_IRismFdirectional(vec_L, vec_rs, vec_as,...
    Qmax, vec_rm, mat_beta, fs, Lh,0);

%% the impulse responses with APIM approach under different Qmax's
vec_Qmax = [-1:1:5].';
mat_h = zeros(Lh,length(vec_Qmax));
for ii = 1:length(vec_Qmax)
    Qmax = vec_Qmax(ii) 
    [hVec] = cal_IRismFdirectional(vec_L, vec_rs, vec_as,...
        Qmax, vec_rm, mat_beta, fs, Lh,0); 
    mat_h(:,ii) = hVec;
end
 
idx= (1:1000).';%/fs*1e3; 
%% show the results
figure; 
%************************************************************** 
subplot(3,1,1) ; 
plot(idx, mat_h(idx,1),'linestyle','-', 'color',...
    [0.7,0.7,0.7], 'linewidth',2); hold on;
plot(idx, mat_h(idx,2),'linestyle','-', 'color', ...
    [0.,0.,1], 'linewidth',0.5,...
    'marker','none', 'markersize',3); hold on; 
grid on;  
xlabel('Time index' );
ylabel('Impulse response'); 
%********************************************************
subplot(3,1,2) 
plot(idx, mat_h(idx,1),'linestyle','-', 'color',...
    [0.7,0.7,0.7], 'linewidth',2); hold on;
plot(idx, mat_h(idx,3),'linestyle','-', 'color', ...
    [0.,0.,1], 'linewidth',0.5,...
    'marker','none', 'markersize',3); hold on; 
grid on;  
xlabel('Time index');
ylabel('Impulse response');  
%********************************************************
subplot(3,1,3) 
plot(idx, mat_h(idx,1),'linestyle','-', 'color',...
    [0.7,0.7,0.7], 'linewidth',2); hold on;
plot(idx, mat_h(idx,4),'linestyle','-', 'color', ...
    [0.,0.,1], 'linewidth',0.5,...
    'marker','none', 'markersize',3); hold on; 
grid on; 
xlabel('Time index');
ylabel('Impulse response'); 






