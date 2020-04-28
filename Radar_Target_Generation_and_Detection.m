clc
close all
clear 
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = 77e9;           % Frequency of operation = 77GHz
R_max = 200;        % Max Range = 200m
R_res = 1;          % Range Resolution = 1 m
v_max = 100;        % Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8;            %speed of light = 3e8
%% User Defined Range and Velocity of target
% define the target's initial position and velocity. Note : Velocity remains contant
R = 110;        % initial range of the target (Cann't exceed the max. range of the radar)
v = -20;        % speed of the target

%% FMCW Waveform Generation

%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

% the range resolution defines the chirp Bandwidth
B = c / (2*R_res);

% the max. range defines the chirp time
Tchirp = 5.5 * (2*R_max)/c;

% the slope of the chirp
slope = B / Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
% linspace(x1,x2,n) generates n points  in the interval [x1,x2]. The spacing between the points is (x2-x1)/(n-1).


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));  % range update
td=zeros(1,length(t));   % time delay


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = R + v * t(i);
    td(i) = (2*r_t(i)) / c;
    
    %For each time sample we need update the transmitted and received signal. 
    Tx(i) = cos(2*pi*(fc*t(i) + ((slope * (t(i))^2) / 2)));
    
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + ((slope*(t(i)-td(i))^2) / 2)));
    
    %Now by mixing(suptracting) the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and Receiver Signal
    %This process in turn works as frequency subtraction.
    Mix(i) = Tx(i) * Rx(i); 
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
% Nr is the size of Range FFT samples (rows = 1024)(The number of samples on each chirp)
% Nd is the size of Doppler FFT samples (columns = 128) (The number of chirps in one sequence)
Mix = reshape(Mix,[Nr,Nd]);  % 1024 x 128

%run the FFT on the beat signal along the range bins dimension (Nr) and normalize.
signal_fft = fft(Mix);       % 1024 x 128
%Y = fft(X) computes the discrete Fourier transform (DFT) of X using a fast Fourier transform (FFT) algorithm.
%1. If X is a vector, then fft(X) returns the Fourier transform of the vector.
%2. If X is a matrix, then fft(X) treats the columns of X as vectors 
%and returns the Fourier transform of each column.

% Normalization
signal_fft = signal_fft ./ Nr;   % 1024 x 128
% Take the absolute value of FFT output
signal_fft = abs(signal_fft);    % 1024 x 128

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:Nr/2);   % 1 x 512  the first column (the first chirp "all chirps the same")

%plotting the range
figure ('Name','Range from First FFT')

% plot FFT output
plot(signal_fft) 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
signal_fft2 = fft2(Mix,Nr,Nd);   % 1024 x 128

% Taking just one side of signal from Range dimension.
signal_fft2 = signal_fft2(1:Nr/2,1:Nd);   % 512 x 128
signal_fft2 = fftshift (signal_fft2);
RDM = abs(signal_fft2);
% To make the signal stength is defined in logarithmic (db)
%RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure('Name','the output of 2DFFT')
surf(doppler_axis,range_axis,RDM);
colorbar;
%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 8;

%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% offset the threshold by SNR value in dB
offset = 8;

%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR


No_training_cells = ((2*Tr+2*Gr+1)*(2*Td+2*Gd+1))-((2*Gr+1)*(2*Gd+1));

for i = Tr+Gr+1:(Nr/2)-(Gr+Tr)
    for j = Td+Gd+1:Nd-(Gd+Td)
        
       %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        
        % Calculate noise SUM in the area around CUT (sum of all training cells)
        for p = i-(Tr+Gr) : i+(Tr+Gr)
            for q = j-(Td+Gd) : j+(Td+Gd)
                if (abs(i-p) > Gr || abs(j-q) > Gd)
                    noise_level = noise_level + RDM(p,q);
                    % if the signal stength is defined in logarithmic (db)
                    % use this equation insteed the above
                    %noise_level = noise_level + db2pow(RDM(p,q));   
                    
                    %y = db2pow(ydb) returns the power measurements, y, 
                    %that correspond to the decibel (dB) values specified in ydb. 
                    %The relationship between power and decibels is ydb = 10 log10(y).
                end
            end
        end
        
        % Calculate threshould from noise average then add the offset
        
        threshold = noise_level/No_training_cells;
        threshold = threshold * offset;
        
        % if the signal stength is defined in logarithmic (db)
        %threshold = pow2db(noise_level/No_training_cells);
        %threshold = threshold + offset;
        
        CUT = RDM(i,j);
        
        if (CUT < threshold)
            RDM(i,j) = 0;
        else
            RDM(i,j) = 1;
        end
        
    end
end


% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 

RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows
RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns 

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure('Name','CA-CFAR Filtered RDM')
surf(doppler_axis,range_axis,RDM);
colorbar;