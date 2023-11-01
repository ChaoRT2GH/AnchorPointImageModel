function [hVec] = cal_IRismFdirectional(vec_L, vec_rs, vec_as, Qmax, vec_rm, mat_beta, fs, Lh,isHIGHPASS)
% cal. the room impulse response according to the image source method, where 
% the fractional time delay problem is handled (the "F" in the function title
% stands for "fractional").  
% Inputs:
%   vec_L = [Lx;Ly;Lz] carries the size of the room, in m
%   vec_rs = [xs;ys;zs] carries the position of the source, in m
%   vec_as = [xa; ya; za], the anchor point of the source, rs-as is the
%           orientation of the source;
%   Qmax, the maximum order of the directional sources
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
% Auther: Dr. Chao Pan, NWPU
% Copyright 2022 (c) 
% Ref. Chao Pan, Lei Zhang, Yilong Lu, Jilu Jin, Lin Qiu, Jingdong Chen,
% Jacob Benesty. An anchor-point based image-model for room impulse
% response simulation with directional source radiation and sensor
% directivity patterns. https://arxiv.org/abs/2308.10543
% ************************************************************************@ 


% param. sound speed, sampling rate, and length of the impulse response
c = 340;     % sound speed, in m/s

% param. for fractional time delay filter
% N2w = 32;    %  2*N2w+1 is the filter length
% vec_wind = hanning(N2w*2+1); % a hanning window to reduce Gibs phenomenon 
% zeta = 0;       % fractional delay
% vec_hD = vec_wind.*sinc((-N2w:1:N2w).'-zeta); % FIR filter to do the fractional delay
% the parameters for the fractictional delay
% Qmax = 3;
D = 32;
vec_ell = (0:1:2*D).';

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
hVec = zeros(Lh+2*D, 1);  % +2*D is for fractional time delay,D for front and end
Qx = floor(Lh*c/(2*Lx*fs)+1);
Qy = floor(Lh*c/(2*Ly*fs)+1);
Qz = floor(Lh*c/(2*Lz*fs)+1);
% the fkHz to cal. the radiation pattern 
De = 64;
vec_fkHz = (0:De).'/(2*De+1)*fs/1000; 
vec_rho = log(1 + 0.6743*vec_fkHz + 0.3776*vec_fkHz.^2 - 0.0540*vec_fkHz.^3+...
       0.020*vec_fkHz.^4);
   

for qx = -Qx:Qx
    for qy = -Qy:Qy
        for qz = -Qz:Qz
            % px,py,pz = 0/1,0/1,0/1
            for px=0:1
                for py = 0:1
                    for pz = 0:1
                        betas = fun_betas(qx,qy,qz,px,py,pz); % the beta
                        dt = fun_dt(qx,qy,qz,px,py,pz);       % distance, in meter 
                        tau = round(dt*fs/c);                 % the index of the reflect
                        zeta = dt*fs/c - tau;   % the fractional delay
                        vec_win = 0.54 - 0.46*cos((pi*(vec_ell-zeta))/D); % the hamming window
                        
                        if abs(qx)<=Qmax && abs(qy)<=Qmax && abs(qz)<=Qmax
                            % 1) cal. the angle of the orientation
                            %cosagl
                            varphi = (-1).^([px;py;pz]).*vec_rs +...
                                2*[qx;qy;qz].*[Lx;Ly;Lz] - vec_rm;
                            vartheta = (-1).^([px;py;pz]).*(vec_rs - vec_as);
                            cosagl = -(varphi'*vartheta)/(norm(vartheta, 2)*...
                                norm(varphi, 2));
                            % 2) cal. the response
                            vec_g = (0.5*(1+cosagl)).^vec_rho;
                            vec_a = (1+vec_fkHz).^(-2)*(0.5*(1-cosagl)).^8;
                            vec_B = vec_a.*(1-vec_g) + vec_g;
                            vec_PH = exp(-1i*2*pi*(0:De).'/(2*De+1)*(zeta+De));
                            vec_B = vec_B.*vec_PH;
                            vec_B2 = [vec_B; flipud(conj(vec_B(2:De+1)))];
                            vec_en = real(ifft(vec_B2, 2*De+1));
                            vec_bn = vec_en(De-D+1:De+D+1);
                        else % no source radiation pattern
                            vec_bn = sinc(vec_ell-D-zeta);
                        end 

                        if tau + D < Lh && tau-D>0
                            hVec(tau-D:tau+D) = hVec(tau-D:tau+D) ...
                                + (betas / (4*pi*dt))*vec_win.*vec_bn;           %
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
hVec = hVec(D+1:D+Lh);

% high pass the impulse response
%[vec_b,vec_a]=butter(3,0.06,'high');
%fvtool(vec_b,vec_a)
vec_b = [0.8280   -2.4839    2.4839   -0.8280]';
vec_a = [1.0000   -2.6236    2.3147   -0.6855]';

% highpass the IR if necessary
if isHIGHPASS ~= 0
    hVec = filter(vec_b, vec_a, hVec);
end




