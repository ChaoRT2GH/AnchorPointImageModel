function [hVec] = cal_IRismF(vec_L, vec_rs, vec_rm, mat_beta, fs, Lh,isHIGHPASS)
% cal. the room impulse response according to the image source method, where 
% the fractional time delay problem is handled (the "F" in the function title
% stands for "fractional").  
% Inputs:
%   vec_L = [Lx;Ly;Lz] carries the size of the room, in m
%   vec_rs = [xs;ys;zs] carries the position of the source, in m
%   vec_rm = [xm;ym;zm] carries the poistion of the microphone, in m
%   mat_beta = [betax0, betax1;
%               betay0, betay1;
%               betaz0, betaz1]; carries the reflection coefficients of six
%               walls.
%   fs, the sampling frequencies, in Hz
%   Lh, the length the impulse response 
%   isHIGHPASS, ~=0, using hipass filter 
% Outputs:
%   hVec, the calculated impulse response, a vector of length Lh x 1.
%
% ************************************************************************@
% Auther: Chao Pan, Ph.D candidate, NWPU
% Copyright 2016 (c)
% ************************************************************************@ 


% param. sound speed, sampling rate, and length of the impulse response
c = 340;     % sound speed, in m/s

% param. for fractional time delay filter
N2w = 32;    %  2*N2w+1 is the filter length
vec_wind = hanning(N2w*2+1); % a hanning window to reduce Gibs phenomenon 
D = 0;       % fractional delay
vec_hD = vec_wind.*sinc((-N2w:1:N2w).'-D); % FIR filter to do the fractional delay

% param. size parameters of the room, in meter
Lx = vec_L(1);
Ly = vec_L(2);
Lz = vec_L(3);

% param. reflection coefficients
betax0 = mat_beta(1,1);
betax1 = mat_beta(1,2);
betay0 = mat_beta(2,1);
betay1 = mat_beta(2,2);
betaz0 = mat_beta(3,1);
betaz1 = mat_beta(3,2);

% functions for calculating the beta and distance
fun_betas = @(qx,qy,qz,px,py,pz)...
    betax0^(abs(qx-px))*betax1^(abs(qx))*...
    betay0^(abs(qy-py))*betay1^(abs(qy))*...
    betaz0^(abs(qz-pz))*betaz1^(abs(qz));
fun_dt = @(qx,qy,qz,px,py,pz)...
    norm(...
    (-1).^([px;py;pz]).*vec_rs + 2*[qx;qy;qz].*[Lx;Ly;Lz] - vec_rm,...
    2);

% cal. the impulse response
hVec = zeros(Lh+2*N2w, 1);  % +2*N2w is for fractional time delay,N2w for front and end
Qx = round(Lh*c/(2*Lx*fs)+1);
Qy = round(Lh*c/(2*Ly*fs)+1);
Qz = round(Lh*c/(2*Lz*fs)+1);

for qx = -Qx:Qx
    for qy = -Qy:Qy
        for qz = -Qz:Qz
            % px,py,pz = 0/1,0/1,0/1
            for px=0:1
                for py = 0:1
                    for pz = 0:1
                        betas = fun_betas(qx,qy,qz,px,py,pz); % the beta
                        dt = fun_dt(qx,qy,qz,px,py,pz);       % distance, in meter 
                        ih = round(dt*fs/c);               % the index of the reflect
                        D = dt*fs/c - ih;   % the fractional delay
                        
                        ih = ih + N2w;   % plus the global delay of the IR
                        vec_hD = vec_wind.*sinc((-N2w:1:N2w).'-D); % FIR filter to do the fractional delay

                        if ih < Lh+N2w
                            hVec(ih-N2w:ih+N2w) = hVec(ih-N2w:ih+N2w) ...
                                + (betas / (4*pi*dt))*vec_hD;           %
                        end
                    end
                end
            end
            %
        end
    end
end
 
%  % add a lowpass filter to the filter out the negative influence in the
%  high frequecy band
% d = fdesign.lowpass('N,Fc,Ap,Ast',16,0.9,0.5,40);
% designmethods(d)
% hd = design(d,'equiripple');
hlowpass=[-0.0056646; 0.031756; -0.033528; 0.050612; -0.065185; 0.079288;...
    -0.090623;0.098038;0.89939;0.098038;-0.090623;0.079288;-0.065185;...
    0.050612;-0.033528;0.031756;-0.0056646];
% fvtool(hlowpass);
hVec = filter(hlowpass, 1, hVec);
hVec = hVec(N2w+1:N2w+Lh);

% high pass the impulse response
%[vec_b,vec_a]=butter(3,0.06,'high');
%fvtool(vec_b,vec_a)
vec_b = [0.8280   -2.4839    2.4839   -0.8280]';
vec_a = [1.0000   -2.6236    2.3147   -0.6855]';

% highpass the IR if necessary
if isHIGHPASS ~= 0
    hVec = filter(vec_b, vec_a, hVec);
end




