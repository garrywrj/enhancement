function GSC()
%==================================================
% The script is for Speech Enhancement.
% A GSC Beamformer add a Postfiler in this script
% Author              :Xiukiun Wu
% Date                :30-06-2017
% Related paper       :
% ===================================================
% try
clear all
% ---------------the paremeter set---------------------
global sc ;
len = 1024;
%fs  = 32e3;
fs  = 16e3;
El  = 0.03;                 % ElemnetSpace                  
Sp  = 340;                  % the speed of the voice in air
%Ne  = 4;                    % the dimesion of microphones
Ne  = 2;                    % the dimesion of microphones
ti  = 1;
ki  = 1;
nFFT     = 2 * len;
len_buf  = len / 2;
%angSteer = [0;0];
angSteer = [85;0];
% --------------------------------------
%y_buf = zeros(len,4);
y_buf = zeros(len,2);
d_buf = zeros(len_buf,1);
ds_buf = zeros(len_buf,1);
win   = hanning(len);  
%--------------Later set------------------
% pcm1 = './实验室1录音数据/主声源90度干扰30度_主声源提前/4mic原始数据/chann0_32k.pcm';
% pcm2 = './实验室1录音数据/主声源90度干扰30度_主声源提前/4mic原始数据/chann1_32k.pcm';
% pcm3 = './实验室1录音数据/主声源90度干扰30度_主声源提前/4mic原始数据/chann2_32k.pcm';
% pcm4 = './实验室1录音数据/主声源90度干扰30度_主声源提前/4mic原始数据/chann3_32k.pcm';
% fileId1 = fopen(pcm1,'r');
% fileId2 = fopen(pcm2,'r');
% fileId3 = fopen(pcm3,'r');
% fileId4 = fopen(pcm4,'r');
% x1 = (fread(fileId1,inf,'int16'))/(2^15);
% x2 = (fread(fileId2,inf,'int16'))/(2^15);
% x3 = (fread(fileId3,inf,'int16'))/(2^15);
% x4 = (fread(fileId4,inf,'int16'))/(2^15);
yf = audioread('M:\论文\4.27_1\10dB\主_130_噪_85_MicroArrayData.wav');

% yf(:,1) = x1;
% yf(:,2) = x2;
% yf(:,3) = x3;
% yf(:,4) = x4;
%-------------------Linear microphone array set-------------------
 himc = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 20e3]);
 sULA = phased.ULA('NumElements',Ne,'ElementSpacing',El,'Element',himc);
% rGSC = phased.GSCBeamformer('SensorArray',sULA,'SampleRate',fs,'PropagationSpeed'...
%                             ,Sp,'Direction',angSteer,'FilterLength',10);
sc.delay   = phased.ElementDelay('SensorArray',sULA,'PropagationSpeed',Sp);                        
filereader = dsp.SignalSource('Signal',yf,'SamplesPerFrame',len);  
sc.fs      = fs;
sc.Ne      = Ne;
sc.L       = 10;
sc.ang     = angSteer;
sc.Wanc1    = zeros(sc.L,Ne-1);
sc.Wanc2    = zeros(sc.L,Ne-1);
Wds   = 1/sc.Ne * ones(sc.Ne,1);
DataBuffer = zeros(sc.L-1,sc.Ne);
while ~filereader.isDone
	audin = step(filereader);   
    if ki ==1
        y_buf(1:len_buf,:) = audin(len_buf+1:len,:);
     else
     y_buf(len_buf + 1:len,:) = audin(1:len_buf,:);
     
     y2     = bsxfun(@times,win,audin);
     y1     = bsxfun(@times,win,y_buf);
     d1     = rDS([DataBuffer;y1]);
     d2     = rDS([DataBuffer;y2]);
%    df1    = fft(d1,nFFT);
%    df2    = fft(d2,nFFT); 
     d_b1   = BM(d1);  %block matrix output
     d_b2   = BM(d2);  %block matrix output
     d_out1 = ANC1(d1,d_b1,y1);
     d_out2 = ANC2(d2,d_b2,y2);
     
     y_buf(1:len_buf,:)            = audin(len_buf + 1 : len,:);
     GSCout(ki:ki+len_buf-1,1)     = d_buf + d_out1(1:len_buf);
     GSCout(ki+len_buf:ki+len-1,1) = d_out2(1:len_buf,1) + d_out1(len_buf+1:len,1);
     DSout(ki:ki+len_buf-1,1)     = ds_buf + d1(1+9:len_buf+9,:)*Wds;
     DSout(ki+len_buf:ki+len-1,1) = d2(1+9:len_buf+9,:)*Wds + d1(len_buf+1+9:len+9,:)*Wds;
     d_buf  = d_out2(len_buf+1:len);
     ds_buf = d2(len_buf+1+9:len+9,:)*Wds;
     end
     ki = ki + len;

end
out_GSC = sprintf('./wav1/%s.wav','GSCout');
out_DS = sprintf('./wav1/%s.wav','DSout');
filewriter1 = dsp.AudioFileWriter(out_GSC,'SampleRate',fs);
filewriter2 = dsp.AudioFileWriter(out_DS,'SampleRate',fs);
step(filewriter1,GSCout);
step(filewriter2,DSout);
release(filewriter1);
release(filewriter2);
end

%--------------Define Time Delay--------------
function ds = rDS(x)  
    global sc;
    td = step(sc.delay,sc.ang);
    ds  = delayseq(x,-td,sc.fs);
end
% --------------Define the block--------------
function ybm = BM(x)
	global sc;
	M  = sc.Ne;
    M2 = nextpow2(M);
    if ~isequal(2^M2,M)
       Ws = toeplitz([1 zeros(1,M-2)],[1 -1 zeros(1,M-2)]);
    else
       OM = hadamard(2^M2);                 %  OrthoMatrix = 
       Ws = OM(2:end,:);
    end        
	x_block = x;
	ybm       = x_block * Ws';               % x_null_block    
end
%-------------Define the ANC--------------------
% function ds = FANC(x)
%      global sc;
%     L         = sc.L;
%     alpha     = 1/L;
% 	yanc      = zeros(size(x1,1),1);
% 	Wds       = 1/sc.Ne * ones(sc.Ne,1);
%     row       = size(x2,1);
%     WF        = zeros(L,1);
%     yc        = zeros(row,1); 
% 	yout      = zeros(size(x,1),1); 
%     Wanc      = sc.Wanc2;
%     yds       = x1 * Wds;                                        %DS output  
%     WF(end) = 1;    
%     for k = 1:size(x,1)
%         
%         
%         
%         
%     end 
%   
% end
function yout = ANC1(x1,x2,x)
    global sc ;
    L         = sc.L;
    alpha     = 1/L;
	yanc      = zeros(size(x1,1),1);
	Wds       = 1/sc.Ne * ones(sc.Ne,1);
    row       = size(x2,1);
    WF        = zeros(L,1);
    yc        = zeros(row,1);
    yout      = zeros(size(x,1),1); 
    yds       = x1 * Wds;                                    %DS output  
    Wanc      = sc.Wanc1; 
    WF(end) = 1;    
    for k = 1:size(x,1)
        xnull     = x2(k:k + L - 1,:);                           %anc input in L
        yc(k,1)   = sum(WF .* yds(k:k + L - 1,1),1);             %anc output in L
        P_xnull   = sum(sum(abs(xnull).^2,1),2);
        epa       = yc(k,1) - sum(sum(Wanc.*xnull),2)/3; 
        Wanc      = Wanc + alpha /P_xnull * conj(epa) * xnull;
        yanc(k,1) = sum(sum(Wanc.*xnull,1),2);
        yout(k,1) = yc(k,1) - yanc(k,1);
    end 
    sc.Wanc1 = Wanc;
end
function yout = ANC2(x1,x2,x)
    global sc;
    L         = sc.L;
    alpha     = 1/L;
	yanc      = zeros(size(x1,1),1);
	Wds       = 1/sc.Ne * ones(sc.Ne,1);
    row       = size(x2,1);
    WF        = zeros(L,1);
    yc        = zeros(row,1); 
	yout      = zeros(size(x,1),1); 
    Wanc      = sc.Wanc2;
    yds       = x1 * Wds;                                        %DS output  
    WF(end) = 1;    
    for k = 1:size(x,1)
        xnull     = x2(k:k + L - 1,:);                           %anc input in L
        yc(k,1)   = sum(WF .* yds(k:k + L - 1,1),1);             %anc output in L
        P_xnull   = sum(sum(abs(xnull).^2,1),2);
        epa       = yc(k,1) - sum(sum(Wanc.*xnull),2); 
        Wanc      = Wanc + alpha /P_xnull * conj(epa) * xnull;
        yanc(k,1) = sum(sum(Wanc.*xnull,1),2);
        yout(k,1) = yc(k,1) - yanc(k,1);
    end 
    sc.Wanc2 = Wanc;
end






