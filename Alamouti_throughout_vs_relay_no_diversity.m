%% Initialising data
clear all
clc
close all

size = 10000;
data1(1,:) = round(rand(1,size));   % generates s0 - signal 1 to be sent
data1(2,:) = round(rand(1,size)); %
data1(3,:) = round(rand(1,size));
data1(4,:) = round(rand(1,size));
data1(~data1)= -1;    %turns vector to ones and negative ones

data2(1,:) = round(rand(1,size));   % generates s1 - signal 2 to be sent
data2(2,:) = round(rand(1,size)); % 
data2(3,:) = round(rand(1,size));
data2(4,:) = round(rand(1,size));
data2(~data2)= -1;    %turns vector to ones and negative ones


pointA1 = complex(3,3);
pointA2 = complex(3,1);
pointA3 = complex(1,3);
pointA4 = complex(1,1);

pointB1 = complex(3,-1);
pointB2 = complex(3,-3);
pointB3 = complex(1,-1);
pointB4 = complex(1,-3);

pointC1 = complex(-1,-1);
pointC2 = complex(-1,-3);
pointC3 = complex(-3,-1);
pointC4 = complex(-3,-3);

pointD1 = complex(-1,3);
pointD2 = complex(-1,1);
pointD3 = complex(-3,3);
pointD4 = complex(-3,1);

%% loop creating 2-d signal position based on data
% QPSK has 4 bits, bits 1 & 2 will rperesent which quadrant
% bits 3 & 4 will represent position within quadrant - values= -3, -1, 1, 3

for loop = 1:size         
    %real part of s0 - real part deals with bits 1 & 3
if ((data1(1,loop) > 0) && (data1(3,loop) >0) )
    Tx1(1,loop) = 3;
end
if ((data1(1,loop) > 0) && (data1(3,loop) <=0) )
    Tx1(1,loop) = 1;
end
if ((data1(1,loop) <= 0) && (data1(3,loop) >0) )
    Tx1(1,loop) = -1;
end
if ((data1(1,loop) <= 0) && (data1(3,loop) <=0) )
    Tx1(1,loop) = -3;
end

    %complex part pf s0 - complex part deals with bits 2 & 4
if ((data1(2,loop) > 0) && (data1(4,loop) >0) )
    Tx1(2,loop) = 3;
end
if ((data1(2,loop) > 0) && (data1(4,loop) <=0) )
    Tx1(2,loop) = 1;
end
if ((data1(2,loop) <= 0) && (data1(4,loop) >0) )
    Tx1(2,loop) = -1;
end
if ((data1(2,loop) <= 0) && (data1(4,loop) <=0) )
    Tx1(2,loop) = -3;
end
%

    %real part of s1
if ((data2(1,loop) > 0) && (data2(3,loop) >0) )
    Tx2(1,loop) = 3;
end
if ((data2(1,loop) > 0) && (data2(3,loop) <=0) )
    Tx2(1,loop) = 1;
end
if ((data2(1,loop) <= 0) && (data2(3,loop) >0) )
    Tx2(1,loop) = -1;
end
if ((data2(1,loop) <= 0) && (data2(3,loop) <=0) )
    Tx2(1,loop) = -3;
end

    %complex part pf s1
if ((data2(2,loop) > 0) && (data2(4,loop) >0) )
    Tx2(2,loop) = 3;
end
if ((data2(2,loop) > 0) && (data2(4,loop) <=0) )
    Tx2(2,loop) = 1;
end
if ((data2(2,loop) <= 0) && (data2(4,loop) >0) )
    Tx2(2,loop) = -1;
end
if ((data2(2,loop) <= 0) && (data2(4,loop) <=0) )
    Tx2(2,loop) = -3;
end
end

%% Transmission Loop to Destination & Relay + Coding & Decoding & erroer checking

s0(1,:) = complex(Tx1(1,:), Tx1(2,:));    %Transmitter is complex of dataX and dataY
s1(1,:) = complex(Tx2(1,:), Tx2(2,:));    %data is a combination of the real and imaginary

Es = 10;  % in this case for points; (1,1), (1,3), (3,1), (3,3) - average symbol energy = 10
pyth = 1/(sqrt(2));     %the real and imaginary parts carry equal weight of the signal 1/root(2)
relay_fading = (10/6)^2;  %in this case the relay networks is set up with a relay node dist = 6, destination node = 10 from source

z=1;            %loop counter
max_db = 20;    %max Es/No (SNR) limit for graph / test - can set here
step = 1;       %step size

for db = 0:step:max_db      %variation of Es/No (SNR)
    db
    No = Es/(10^(db/10));   %calculation of noise using known SNR
    var = No/2;             %variance = sig^2 = No/2
    std = sqrt(var);        %standar deviation

    %noise functions are Gaussian functions (randn) - different noise for
    %each channel
    n0(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    n1(1,:) = std.*randn(1,size) + std.*1i*randn(1,size);
    
    n2(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    n3(1,:) = std.*randn(1,size) + std.*1i*randn(1,size);

    %Fading functions are Rayleigh functions (2 Gaussians = complex) -
    %different fading for each channel
    h0 = (pyth).*randn(1,size)+ (pyth).*1i*randn(1,size); %complex fading
    h1 = (pyth).*randn(1,size)+ (pyth).*1i*randn(1,size);
    
    h2 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size); %complex fading
    h3 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size);
    
    %Alamoutis equations for transmitting data
    r0(1,:) = h0(1,:).*s0(1,:) + h1(1,:).*s1(1,:) + n0(1,:);  %alamouti eqns
    r1(1,:) = -h0(1,:).*conj(s1(1,:)) + h1(1,:).*conj(s0(1,:)) + n1(1,:);
    
    r2(1,:) = h2(1,:).*s0(1,:) + h3(1,:).*s1(1,:) + n2(1,:);  %alamouti eqns
    r3(1,:) = -h2(1,:).*conj(s1(1,:)) + h3(1,:).*conj(s0(1,:)) + n3(1,:);
    
    %Alamoutis Eqns for resolving data at Rx
    s0_tilda(1,:) = conj(h0(1,:)).*r0(1,:) + h1(1,:).*conj(r1(1,:)); %alamouti eqns
    s1_tilda(1,:) = conj(h1(1,:)).*r0(1,:) - h0(1,:).*conj(r1(1,:));
    
    %Relay resolving data
    s0_tilda_Relay(1,:) = conj(h2(1,:)).*r2(1,:) + h3(1,:).*conj(r3(1,:)); %alamouti eqns
    s1_tilda_Relay(1,:) = conj(h3(1,:)).*r2(1,:) - h2(1,:).*conj(r3(1,:));
    
    % magnitude of fading - fading assumed to be known at the receiver.
    h_mag(1,:) = (abs(h0(1,:)).^2) + (abs(h1(1,:)).^2); %(alpha0^2 + alpha1^2)
    h_mag_Relay(1,:) = (abs(h2(1,:)).^2) + (abs(h3(1,:)).^2); %(alpha0^2 + alpha1^2)
    
    % can remove fading at the receiver as it is assumed known
    s0_rx(1,:) = s0_tilda(1,:)./h_mag(1,:);
    s1_rx(1,:) = s1_tilda(1,:)./h_mag(1,:);
    
    s0_rx_Relay(1,:) = s0_tilda_Relay(1,:)./h_mag_Relay(1,:); %Relay removing fading
    s1_rx_Relay(1,:) = s1_tilda_Relay(1,:)./h_mag_Relay(1,:);    
   

    %no diversity    
%     faded(1,:) = s0(1,:) .* h0(1,:); 
%     
%     Rx(1,:) = faded(1,:) + n0(1,:);
%     Rx(1,:) = Rx(1,:) ./ h0(1,:);
    
    %end of no diversity
    
    
    % 2 sets of noise and fading for the transmission of both s0 & s1 
    % These are channels from relay to destination
    n4(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    h4 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size);
    
    n5(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    h5 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size);
    
     h_mag_yo(1,:) = (abs(h4(1,:)).^2) + (abs(h5(1,:)).^2); %(alpha0^2 + alpha1^2)
    %initialising counts
      count = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    count5 = 0;
    count6 = 0;
    count13 = 0;
    count14 = 0;
    count15= 0;
    count16 = 0;
    count17 = 0;
    count18 = 0;
    count19 = 0;
    count20 = 0;
    
    
    %loop to turn estimate data values & check errors for given SNR
    for i = 1:size  
        
        % 
        % Data on link from source straight to destination - s0
        %
        if (real(s0_rx(1,i)) > 0) %bits 1 & 3 on Re axis
            output(1,i) = 1;        %check faded real
            if (real(s0_rx(1,i)) > 2)
                output(3,i) = 1; 
            else
                output(3,i) = -1; 
            end
        else
            output(1,i) = -1;
            if (real(s0_rx(1,i)) > -2)
                output(3,i) = 1; 
            else
                output(3,i) = -1; 
            end
        end
         if (imag(s0_rx(1,i)) > 0)     %check faded imaginary
            output(2,i) = 1;
            if (imag(s0_rx(1,i)) > 2)
                output(4,i) = 1;
            else
                output(4,i) = -1;
            end
        else
            output(2,i) = -1;
            if (imag(s0_rx(1,i)) > -2)
                output(4,i) = 1;
            else
                output(4,i) = -1;
            end
         end
        
       
        %error checking by comparing to intiial data
         
        if (data1(1,i) == output(1,i))      
            count = count + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output(2,i))
            count = count + 1;              %count will get b-e-r
        end
        if (data1(3,i) == output(3,i))
            count = count + 1;              %count will get b-e-r
        end
        if (data1(4,i) == output(4,i))
            count = count + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output(1,i)) | (data1(2,i) ~= output(2,i)) | (data1(3,i) ~= output(3,i)) | (data1(4,i) ~= output(4,i)) 
            count2 = count2 + 1;            %count 2 will get s-e-r
        end
          
        % 
        % Data on link from source straight to destination - s1
        %
         if (real(s1_rx(1,i)) > 0)
             output2(1,i) = 1;        %check faded real
            if (real(s1_rx(1,i)) > 2)
                output2(3,i) = 1; 
            else
                output2(3,i) = -1; 
            end
        else
            output2(1,i) = -1;
            if (real(s1_rx(1,i)) > -2)
                output2(3,i) = 1; 
            else
                output2(3,i) = -1; 
            end
        end
         if (imag(s1_rx(1,i)) > 0)     %check faded imaginary
            output2(2,i) = 1;
            if (imag(s1_rx(1,i)) > 2)
               output2(4,i) = 1; 
            else
               output2(4,i) = -1; 
            end
        else
            output2(2,i) = -1;
             if (imag(s1_rx(1,i)) > -2)
               output2(4,i) = 1; 
            else
               output2(4,i) = -1; 
            end
         end
        
         % error checking by comparing to initial data2
         %
        if (data2(1,i) == output2(1,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == output2(2,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == output2(3,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == output2(4,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= output2(1,i)) | (data2(2,i) ~= output2(2,i)) | (data2(3,i) ~= output2(3,i)) | (data2(4,i) ~= output2(4,i))
            count4 = count4 + 1;            %count 4 will get s-e-r
        end
        
        %no diversity section - from source to destination
%         if (real(Rx(1,i)) > 0)
%             output3(1,i) = 1;        %check faded real
%         else
%             output3(1,i) = -1;
%         end
%          if (imag(Rx(1,i)) > 0)     %check faded imaginary
%             output3(2,i) = 1;
%         else
%             output3(2,i) = -1;
%         end
%         if (data1(1,i) == output3(1,i))
%             count5 = count5 + 1;              %count will get b-e-r
%         end
%         if (data1(2,i) == output3(2,i))
%             count5 = count5 + 1;              %count will get b-e-r
%         end
%         
%         if (data1(1,i) ~= output3(1,i)) | (data1(2,i) ~= output3(2,i))
%             count6 = count6 + 1;            %count 2 will get s-e-r
%         end
%           

        
        %
        %
        % Relay node decode and re-transmit.
        %
        %
        
        
        % 
        % Data on link from source to relay - s0
        % Re_Tx1 will be used to store the 2-d point estimates
        %
         if (real(s0_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
                 
            if (real(s0_rx_Relay(1,i)) > 2)
                Re_Tx1(1,i) = 3;
            else
                Re_Tx1(1,i) = 1; 
            end
        else
            
            if (real(s0_rx_Relay(1,i)) > -2)
                Re_Tx1(1,i) = -1;
            else
                Re_Tx1(1,i) = -3;
            end
        end
         if (imag(s0_rx_Relay(1,i)) > 0)     %check faded imaginary
            
            if (imag(s0_rx_Relay(1,i)) > 2)
                Re_Tx1(2,i) = 3;
            else
                Re_Tx1(2,i) = 1;
            end
        else
            
            if (imag(s0_rx_Relay(1,i)) > -2)
                Re_Tx1(2,i) = -1;
            else
                Re_Tx1(2,i) = -3;
            end
         end
         
         
        % 
        % Data on link from source to relay - s1
        % Re_Tx2 will be used to store the 2-d point estimates
        %
         
         if (real(s1_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
                 
            if (real(s1_rx_Relay(1,i)) > 2)
                Re_Tx2(1,i) = 3;
            else
                Re_Tx2(1,i) = 1; 
            end
        else
            
            if (real(s1_rx_Relay(1,i)) > -2)
                Re_Tx2(1,i) = -1;
            else
                Re_Tx2(1,i) = -3;
            end
        end
         if (imag(s1_rx_Relay(1,i)) > 0)     %check faded imaginary
            
            if (imag(s1_rx_Relay(1,i)) > 2)
                Re_Tx2(2,i) = 3;
            else
                Re_Tx2(2,i) = 1;
            end
        else
            
            if (imag(s1_rx_Relay(1,i)) > -2)
                Re_Tx2(2,i) = -1;
            else
                Re_Tx2(2,i) = -3;
            end
         end
         
         
         %
         % MIDPOINT CALCS 
         % At Relay node - check accuracy of data received so far.
         % Signals decoded & estimated & estimates checked for errors vs.
         % original signal
         %
          if (real(s0_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
            output7(1,i) = 1;        %check faded real
            if (real(s0_rx_Relay(1,i)) > 2)
                output7(3,i) = 1; 
            else
                output7(3,i) = -1; 
            end
        else
            output7(1,i) = -1;
            if (real(s0_rx_Relay(1,i)) > -2)
                output7(3,i) = 1; 
            else
                output7(3,i) = -1; 
            end
        end
         if (imag(s0_rx_Relay(1,i)) > 0)     %bits 2 & 4 on Im axis
            output7(2,i) = 1;
            if (imag(s0_rx_Relay(1,i)) > 2)
                output7(4,i) = 1;
            else
                output7(4,i) = -1;
            end
        else
            output7(2,i) = -1;
            if (imag(s0_rx_Relay(1,i)) > -2)
                output7(4,i) = 1;
            else
                output7(4,i) = -1;
            end
         end
        
       
        % error checking - count 13 will track errors for here
         
        if (data1(1,i) == output7(1,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output7(2,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        if (data1(3,i) == output7(3,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        if (data1(4,i) == output7(4,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output7(1,i)) | (data1(2,i) ~= output7(2,i)) | (data1(3,i) ~= output7(3,i)) | (data1(4,i) ~= output7(4,i)) 
            count14 = count14 + 1;            %count 2 will get s-e-r
        end
          
        
        % Still at Relay node Rx, checkin s1 now.
        
        if (real(s1_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
            output8(1,i) = 1;        %check faded real
            if (real(s1_rx_Relay(1,i)) > 2)
                output8(3,i) = 1; 
            else
                output8(3,i) = -1; 
            end
        else
            output8(1,i) = -1;
            if (real(s1_rx_Relay(1,i)) > -2)
                output8(3,i) = 1; 
            else
                output8(3,i) = -1; 
            end
        end
         if (imag(s1_rx_Relay(1,i)) > 0)     %check faded imaginary
            output8(2,i) = 1;
            if (imag(s1_rx_Relay(1,i)) > 2)
                output8(4,i) = 1;
            else
                output8(4,i) = -1;
            end
        else
            output8(2,i) = -1;
            if (imag(s1_rx_Relay(1,i)) > -2)
                output8(4,i) = 1;
            else
                output8(4,i) = -1;
            end
         end
        
       %count 15 records errors for s1 rx at relay node.       
         
        if (data2(1,i) == output8(1,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == output8(2,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == output8(3,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == output8(4,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= output8(1,i)) | (data2(2,i) ~= output8(2,i)) | (data2(3,i) ~= output8(3,i)) | (data2(4,i) ~= output8(4,i)) 
            count16 = count16 + 1;            %count 2 will get s-e-r
        end
        
    end
    
    s4(1,:) = complex(Re_Tx1(1,:), Re_Tx1(2,:));    % Relay transmitter is complex of dataX and dataY
    s5(1,:) = complex(Re_Tx2(1,:), Re_Tx2(2,:));    % data estimates ready for transmission now ro destination
    
    
    faded_Re1(1,:) = s4(1,:) .* h4(1,:); 
    faded_Re2(1,:) = s5(1,:) .* h5(1,:); 
        
    %noise added, and fading undone.
    Rx1(1,:) = faded_Re1(1,:) + n4(1,:);
    Rx1(1,:) = Rx1(1,:) ./ h4(1,:);
    
    %noise added, and fading undone.
    Rx2(1,:) = faded_Re2(1,:) + n5(1,:);
    Rx2(1,:) = Rx2(1,:) ./ h5(1,:);  
    
    %end of no diversity
    
    count7 = 0;
    count8 = 0;
    count9 = 0;
    count10 = 0;
    count11 = 0;
    count12 = 0;
    %      HERE
    
     for i = 1:size
          
        if (real(Rx1(1,i)) > 0) %bits 1 & 3 on Re axis
            output4(1,i) = 1;        %check faded real
            if (real(Rx1(1,i)) > 2)
                output4(3,i) = 1; 
            else
                output4(3,i) = -1; 
            end
        else
            output4(1,i) = -1;
            if (real(Rx1(1,i)) > -2)
                output4(3,i) = 1; 
            else
                output4(3,i) = -1; 
            end
        end
         if (imag(Rx1(1,i)) > 0)     %check faded imaginary
            output4(2,i) = 1;
            if (imag(Rx1(1,i)) > 2)
                output4(4,i) = 1;
            else
                output4(4,i) = -1;
            end
        else
            output4(2,i) = -1;
            if (imag(Rx1(1,i)) > -2)
                output4(4,i) = 1;
            else
                output4(4,i) = -1;
            end
         end
        
       
        
         
        if (data1(1,i) == output4(1,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output4(2,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        if (data1(3,i) == output4(3,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        if (data1(4,i) == output4(4,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output4(1,i)) | (data1(2,i) ~= output4(2,i)) | (data1(3,i) ~= output4(3,i)) | (data1(4,i) ~= output4(4,i)) 
            count8 = count8 + 1;            %count 2 will get s-e-r
        end
          
        
        if (real(Rx2(1,i)) > 0) %bits 1 & 3 on Re axis
            output5(1,i) = 1;        %check faded real
            if (real(Rx2(1,i)) > 2)
                output5(3,i) = 1; 
            else
                output5(3,i) = -1; 
            end
        else
            output5(1,i) = -1;
            if (real(Rx2(1,i)) > -2)
                output5(3,i) = 1; 
            else
                output5(3,i) = -1; 
            end
        end
         if (imag(Rx2(1,i)) > 0)     %check faded imaginary
            output5(2,i) = 1;
            if (imag(Rx2(1,i)) > 2)
                output5(4,i) = 1;
            else
                output5(4,i) = -1;
            end
        else
            output5(2,i) = -1;
            if (imag(Rx2(1,i)) > -2)
                output5(4,i) = 1;
            else
                output5(4,i) = -1;
            end
         end
        
       
        
         
        if (data2(1,i) == output5(1,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == output5(2,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == output5(3,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == output5(4,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= output5(1,i)) | (data2(2,i) ~= output5(2,i)) | (data2(3,i) ~= output5(3,i)) | (data2(4,i) ~= output5(4,i)) 
            count10 = count10 + 1;            %count 2 will get s-e-r
        end
        
        
        %COMBINING 2 SIGNALS USING EUCLIDIAN DISTANCE
        
        %output = s0_rx (s0 - source to dest)
        %output4 = Rx1 (s0 - relay to dest)
        %if their estimates are not equal - go Euclidian distance and set a
        %final output
        
        source1(1,:) = s0_rx(1,:);
        relay1(1,:) = Rx1(1,:);
        
         source2(1,:) = s1_rx(1,:);
        relay2(1,:) = Rx2(1,:);
        
       
        
        if (output4(1,i) ~= output(1,i)) | (output4(2,i) ~= output(2,i)) | (output4(3,i) ~= output(3,i)) | (output4(4,i) ~= output(4,i)) 
            
            
            %WEIGHTED RELAY AT DESTINATION COMBINER
            dist_s0(1,1) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA1)^2;
            dist_s0(1,2) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA2)^2;
            dist_s0(1,3) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA3)^2;
            dist_s0(1,4) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA4)^2;
            dist_s0(1,5) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB1)^2;
            dist_s0(1,6) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB2)^2;
            dist_s0(1,7) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB3)^2;
            dist_s0(1,8) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB4)^2;
            dist_s0(1,9) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC1)^2;
            dist_s0(1,10) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC2)^2;
            dist_s0(1,11) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC3)^2;
            dist_s0(1,12) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC4)^2;
            dist_s0(1,13) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD1)^2;
            dist_s0(1,14) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD2)^2;
            dist_s0(1,15) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD3)^2;
            dist_s0(1,16) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD4)^2;
            
            
            %NO WEIGHTINGS AT COMBINER
%             dist_s0(1,1) = abs(source1(1,i) - pointA1)^2 + abs(relay1(1,i) - pointA1)^2;
%             dist_s0(1,2) = abs(source1(1,i) - pointA2)^2 + abs(relay1(1,i) - pointA2)^2;
%             dist_s0(1,3) = abs(source1(1,i) - pointA3)^2 + abs(relay1(1,i) - pointA3)^2;
%             dist_s0(1,4) = abs(source1(1,i) - pointA4)^2 + abs(relay1(1,i) - pointA4)^2;
%             dist_s0(1,5) = abs(source1(1,i) - pointB1)^2 + abs(relay1(1,i) - pointB1)^2;
%             dist_s0(1,6) = abs(source1(1,i) - pointB2)^2 + abs(relay1(1,i) - pointB2)^2;
%             dist_s0(1,7) = abs(source1(1,i) - pointB3)^2 + abs(relay1(1,i) - pointB3)^2;
%             dist_s0(1,8) = abs(source1(1,i) - pointB4)^2 + abs(relay1(1,i) - pointB4)^2;
%             dist_s0(1,9) = abs(source1(1,i) - pointC1)^2 + abs(relay1(1,i) - pointC1)^2;
%             dist_s0(1,10) = abs(source1(1,i) - pointC2)^2 + abs(relay1(1,i) - pointC2)^2;
%             dist_s0(1,11) = abs(source1(1,i) - pointC3)^2 + abs(relay1(1,i) - pointC3)^2;
%             dist_s0(1,12) = abs(source1(1,i) - pointC4)^2 + abs(relay1(1,i) - pointC4)^2;
%             dist_s0(1,13) = abs(source1(1,i) - pointD1)^2 + abs(relay1(1,i) - pointD1)^2;
%             dist_s0(1,14) = abs(source1(1,i) - pointD2)^2 + abs(relay1(1,i) - pointD2)^2;
%             dist_s0(1,15) = abs(source1(1,i) - pointD3)^2 + abs(relay1(1,i) - pointD3)^2;
%             dist_s0(1,16) = abs(source1(1,i) - pointD4)^2 + abs(relay1(1,i) - pointD4)^2;
            
            [mini, posn] = min(dist_s0);
            if posn == 1
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 1; %A1
            end
            if posn == 2
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 0; %A2
            end
            if posn == 3
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 1; %A3
            end
            if posn == 4
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 0; %A4
            end
            if posn == 5
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 1;
            end
            if posn == 6
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 0;
            end
            if posn == 7
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 1;
            end
            if posn == 8
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 0;
            end
            if posn == 9
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 1;
            end
            if posn == 10
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 0;
            end
            if posn == 11
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 1;
            end
            if posn == 12
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 0;
            end
            if posn == 13
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 1;
            end
            if posn == 14
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 0;
            end
            if posn == 15
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 1;
            end
            if posn == 16
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 0;
            end
            Combo(~Combo)= -1;
        else 
            Combo(1,i) = output4(1,i);
            Combo(2,i) = output4(2,i);
            Combo(3,i) = output4(3,i);
            Combo(4,i) = output4(4,i);
        end
        
        
        
        %output2 = s1_rx (s1 - source to dest)
        %output5 = Rx2 (s1 - relay to dest)
        %if their estimates are not equal - go Euclidian distance and set a
        %final output
       
        
        if (output5(1,i) ~= output2(1,i)) | (output5(2,i) ~= output2(2,i)) | (output5(3,i) ~= output2(3,i)) | (output5(4,i) ~= output2(4,i)) 
            
            %WEIGHTED RELAY AT DESTINATION COMBINER
            dist_s1(1,1) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA1)^2;
            dist_s1(1,2) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA2)^2;
            dist_s1(1,3) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA3)^2;
            dist_s1(1,4) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA4)^2;
            dist_s1(1,5) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB1)^2;
            dist_s1(1,6) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB2)^2;
            dist_s1(1,7) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB3)^2;
            dist_s1(1,8) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB4)^2;
            dist_s1(1,9) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC1)^2;
            dist_s1(1,10) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC2)^2;
            dist_s1(1,11) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC3)^2;
            dist_s1(1,12) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC4)^2;
            dist_s1(1,13) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD1)^2;
            dist_s1(1,14) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD2)^2;
            dist_s1(1,15) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD3)^2;
            dist_s1(1,16) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD4)^2;
 
            
            %NO WEIGHTINGS ON COMBINER
%             dist_s1(1,1) = abs(source2(1,i) - pointA1)^2 + abs(relay2(1,i) - pointA1)^2;
%             dist_s1(1,2) = abs(source2(1,i) - pointA2)^2 + abs(relay2(1,i) - pointA2)^2;
%             dist_s1(1,3) = abs(source2(1,i) - pointA3)^2 + abs(relay2(1,i) - pointA3)^2;
%             dist_s1(1,4) = abs(source2(1,i) - pointA4)^2 + abs(relay2(1,i) - pointA4)^2;
%             dist_s1(1,5) = abs(source2(1,i) - pointB1)^2 + abs(relay2(1,i) - pointB1)^2;
%             dist_s1(1,6) = abs(source2(1,i) - pointB2)^2 + abs(relay2(1,i) - pointB2)^2;
%             dist_s1(1,7) = abs(source2(1,i) - pointB3)^2 + abs(relay2(1,i) - pointB3)^2;
%             dist_s1(1,8) = abs(source2(1,i) - pointB4)^2 + abs(relay2(1,i) - pointB4)^2;
%             dist_s1(1,9) = abs(source2(1,i) - pointC1)^2 + abs(relay2(1,i) - pointC1)^2;
%             dist_s1(1,10) = abs(source2(1,i) - pointC2)^2 + abs(relay2(1,i) - pointC2)^2;
%             dist_s1(1,11) = abs(source2(1,i) - pointC3)^2 + abs(relay2(1,i) - pointC3)^2;
%             dist_s1(1,12) = abs(source2(1,i) - pointC4)^2 + abs(relay2(1,i) - pointC4)^2;
%             dist_s1(1,13) = abs(source2(1,i) - pointD1)^2 + abs(relay2(1,i) - pointD1)^2;
%             dist_s1(1,14) = abs(source2(1,i) - pointD2)^2 + abs(relay2(1,i) - pointD2)^2;
%             dist_s1(1,15) = abs(source2(1,i) - pointD3)^2 + abs(relay2(1,i) - pointD3)^2;
%             dist_s1(1,16) = abs(source2(1,i) - pointD4)^2 + abs(relay2(1,i) - pointD4)^2;
            
            
            [mini2, posn2] = min(dist_s1);
            if posn2 == 1
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 1; %A1
            end
            if posn2 == 2
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 0; %A2
            end
            if posn2 == 3
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 1; %A3
            end
            if posn2 == 4
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 0; %A4
            end
            if posn2 == 5
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 1;
            end
            if posn2 == 6
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 0;
            end
            if posn2 == 7
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 1;
            end
            if posn2 == 8
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 0;
            end
            if posn2 == 9
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 1;
            end
            if posn2 == 10
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 0;
            end
            if posn2 == 11
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 1;
            end
            if posn2 == 12
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 0;
            end
            if posn2 == 13
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 1;
            end
            if posn2 == 14
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 0;
            end
            if posn2 == 15
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 1;
            end
            if posn2 == 16
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 0;
            end
            Combo2(~Combo2)= -1;
            
            
        else 
            Combo2(1,i) = output5(1,i);
            Combo2(2,i) = output5(2,i);
            Combo2(3,i) = output5(3,i);
            Combo2(4,i) = output5(4,i);
        end
        
        if (data1(1,i) == Combo(1,i))      
            count17 = count17 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == Combo(2,i))
            count17 = count17 + 1;              %count will get b-e-r
        end
        if (data1(3,i) == Combo(3,i))
            count17 = count17 + 1;              %count will get b-e-r
            
        end
        if (data1(4,i) == Combo(4,i))
            count17 = count17 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= Combo(1,i)) | (data1(2,i) ~= Combo(2,i)) | (data1(3,i) ~= Combo(3,i)) | (data1(4,i) ~= Combo(4,i)) 
            count18 = count18 + 1;            %count 2 will get s-e-r
        end
        
        
        
        
        if (data2(1,i) == Combo2(1,i))      
            count19 = count19 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == Combo2(2,i))
            count19 = count19 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == Combo2(3,i))
            count19 = count19 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == Combo2(4,i))
            count19 = count19 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= Combo2(1,i)) | (data2(2,i) ~= Combo2(2,i)) | (data2(3,i) ~= Combo2(3,i)) | (data2(4,i) ~= Combo2(4,i)) 
            count20 = count20 + 1;            %count 2 will get s-e-r
        end
       
    end
    
    ber_s0(1,z) = ((size*4) - count) / (size*4);  
    ser_s0(1,z) = (count2 / size);
    count_total_s0(1,z) = count;   %two count vectors in case of errors
    count_total_s0_2(1,z) = count2;
    
    
    %for non-fading - bit error-rate & symbol error-rate
     ber_s1(1,z) = ((size*4) - count3) / (size*4);
     ser_s1(1,z) = (count4 / size);
     count_total_s1(1,z) = count3;
     count_total_s1_2(1,z) = count4;
     
    %for no diversity
     ber_no_diversity(1,z) = ((size*2) - count5) / (size*2);
     ser_no_diversity(1,z) = (count6 / size);
     count_total_no_diversity(1,z) = count5;
     count_total_no_diversity_2(1,z) = count6;
       
    
    ber_relay_s0_mp(1,z) = ((size*4) - count13) / (size*4);
    ber_relay_s1_mp(1,z) = ((size*4) - count15) / (size*4);
    count_total_relay_ber_s0_mp(1,z) = count13;
    count_total_relay_ser_s0_mp(1,z) = count14;
    count_total_relay_ber_s1_mp(1,z) = count15;
    count_total_relay_ser_s1_mp(1,z) = count16;
    
    ber_relay_s0(1,z) = ((size*4) - count7) / (size*4);
    ber_relay_s1(1,z) = ((size*4) - count9) / (size*4);
    count_total_relay_ber_s0(1,z) = count7;
    count_total_relay_ser_s0(1,z) = count8;
    count_total_relay_ber_s1(1,z) = count9;
    count_total_relay_ser_s1(1,z) = count10;
    
    
    %for fading - bit error-rate & symbol error-rate
      
     
    z = z +1;
      
end
    
%% Plotting

ber_alamouti(1,:) = (ber_s0(1,:) + ber_s1(1,:))/2;
ber_relay(1,:) = (ber_relay_s0(1,:) + ber_relay_s1(1,:))/2;

db = [0:step:max_db];
%axis([0 15 10e-10 10e-1]);
% semilogy (db, ber_alamouti, 'r');


%semilogy (db, ber_s1,'r');
semilogy (db,ber_relay,'g');
title ('Alamouti Code Error Rate: 16-QAM Relaying');
xlabel ('SNR of source-destination channel (dB)');
ylabel ('Error rate');

hold on;
% legend('Alamouti QPSK - direct', 'No diversity QPSK - Relay');

%% 2 GRAPHS YO
clear all
clc
% close all

size = 200000;
data1(1,:) = round(rand(1,size));   % generates s0
data1(2,:) = round(rand(1,size)); %
data1(3,:) = round(rand(1,size));
data1(4,:) = round(rand(1,size));
data1(~data1)= -1;    %turns vector to ones and negative ones

data2(1,:) = round(rand(1,size));   % generates s1
data2(2,:) = round(rand(1,size)); % 
data2(3,:) = round(rand(1,size));
data2(4,:) = round(rand(1,size));
data2(~data2)= -1;    %turns vector to ones and negative ones


pointA1 = complex(3,3);
pointA2 = complex(3,1);
pointA3 = complex(1,3);
pointA4 = complex(1,1);

pointB1 = complex(3,-1);
pointB2 = complex(3,-3);
pointB3 = complex(1,-1);
pointB4 = complex(1,-3);

pointC1 = complex(-1,-1);
pointC2 = complex(-1,-3);
pointC3 = complex(-3,-1);
pointC4 = complex(-3,-3);

pointD1 = complex(-1,3);
pointD2 = complex(-1,1);
pointD3 = complex(-3,3);
pointD4 = complex(-3,1);

Tx1 = zeros(2,size);
Tx2 = zeros(2,size);

for loop = 1:size
    %real part of s0
if ((data1(1,loop) > 0) && (data1(3,loop) >0) )
    Tx1(1,loop) = 3;
end
if ((data1(1,loop) > 0) && (data1(3,loop) <=0) )
    Tx1(1,loop) = 1;
end
if ((data1(1,loop) <= 0) && (data1(3,loop) >0) )
    Tx1(1,loop) = -1;
end
if ((data1(1,loop) <= 0) && (data1(3,loop) <=0) )
    Tx1(1,loop) = -3;
end

    %complex part pf s0
if ((data1(2,loop) > 0) && (data1(4,loop) >0) )
    Tx1(2,loop) = 3;
end
if ((data1(2,loop) > 0) && (data1(4,loop) <=0) )
    Tx1(2,loop) = 1;
end
if ((data1(2,loop) <= 0) && (data1(4,loop) >0) )
    Tx1(2,loop) = -1;
end
if ((data1(2,loop) <= 0) && (data1(4,loop) <=0) )
    Tx1(2,loop) = -3;
end
%

    %real part of s1
if ((data2(1,loop) > 0) && (data2(3,loop) >0) )
    Tx2(1,loop) = 3;
end
if ((data2(1,loop) > 0) && (data2(3,loop) <=0) )
    Tx2(1,loop) = 1;
end
if ((data2(1,loop) <= 0) && (data2(3,loop) >0) )
    Tx2(1,loop) = -1;
end
if ((data2(1,loop) <= 0) && (data2(3,loop) <=0) )
    Tx2(1,loop) = -3;
end

    %complex part pf s1
if ((data2(2,loop) > 0) && (data2(4,loop) >0) )
    Tx2(2,loop) = 3;
end
if ((data2(2,loop) > 0) && (data2(4,loop) <=0) )
    Tx2(2,loop) = 1;
end
if ((data2(2,loop) <= 0) && (data2(4,loop) >0) )
    Tx2(2,loop) = -1;
end
if ((data2(2,loop) <= 0) && (data2(4,loop) <=0) )
    Tx2(2,loop) = -3;
end
end

output = zeros(4,size);
output2 = zeros(4,size);
output3 = zeros(2,size);
output4 = zeros(4,size);
output5 = zeros(4,size);
output7 = zeros(4,size);
output8 = zeros(4,size);
Re_Tx1 = zeros(2,size);
Re_Tx2 = zeros(2,size);
Combo = zeros(4,size);
Combo2 = zeros(4,size);

 
s0(1,:) = complex(Tx1(1,:), Tx1(2,:));    %Transmitter is complex of dataX and dataY
s1(1,:) = complex(Tx2(1,:), Tx2(2,:)); 

Es = 10;  
pyth = 1/(sqrt(2));
relay_fading = (10/6)^2;

z=1;
max_db = 20;
step = 3;

for db = 0:step:max_db   %variation of Es/No
    No = Es/(10^(db/10));     %5dB    %incorrect implemenation?
    var = No/2;     %variance = sig^2 = No/2
    std = sqrt(var);
    db
    
    n0(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    n1(1,:) = std.*randn(1,size) + std.*1i*randn(1,size);
    
    n2(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    n3(1,:) = std.*randn(1,size) + std.*1i*randn(1,size);

    h0 = (pyth).*randn(1,size)+ (pyth).*1i*randn(1,size); %complex fading
    h1 = (pyth).*randn(1,size)+ (pyth).*1i*randn(1,size);
    
    h2 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size); %complex fading
    h3 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size);
    
    r0(1,:) = h0(1,:).*s0(1,:) + h1(1,:).*s1(1,:) + n0(1,:);  %alamouti eqns
    r1(1,:) = -h0(1,:).*conj(s1(1,:)) + h1(1,:).*conj(s0(1,:)) + n1(1,:);
    
    r2(1,:) = h2(1,:).*s0(1,:) + h3(1,:).*s1(1,:) + n2(1,:);  %alamouti eqns
    r3(1,:) = -h2(1,:).*conj(s1(1,:)) + h3(1,:).*conj(s0(1,:)) + n3(1,:);
    
    s0_tilda(1,:) = conj(h0(1,:)).*r0(1,:) + h1(1,:).*conj(r1(1,:)); %alamouti eqns
    s1_tilda(1,:) = conj(h1(1,:)).*r0(1,:) - h0(1,:).*conj(r1(1,:));
    
    s0_tilda_Relay(1,:) = conj(h2(1,:)).*r2(1,:) + h3(1,:).*conj(r3(1,:)); %alamouti eqns
    s1_tilda_Relay(1,:) = conj(h3(1,:)).*r2(1,:) - h2(1,:).*conj(r3(1,:));
    
    h_mag(1,:) = (abs(h0(1,:)).^2) + (abs(h1(1,:)).^2); %(alpha0^2 + alpha1^2)
    h_mag_Relay(1,:) = (abs(h2(1,:)).^2) + (abs(h3(1,:)).^2); %(alpha0^2 + alpha1^2)
    
    s0_rx(1,:) = s0_tilda(1,:)./h_mag(1,:);
    s1_rx(1,:) = s1_tilda(1,:)./h_mag(1,:);
    
    s0_rx_Relay(1,:) = s0_tilda_Relay(1,:)./h_mag_Relay(1,:);
    s1_rx_Relay(1,:) = s1_tilda_Relay(1,:)./h_mag_Relay(1,:);    

    %no diversity    
    faded(1,:) = s0(1,:) .* h0(1,:); 
    
    Rx(1,:) = faded(1,:) + n0(1,:);
    Rx(1,:) = Rx(1,:) ./ h0(1,:);
    
    %end of no diversity
    
    
    % CHANNEL INFO FOR RELAY
    n4(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    h4 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size);
    
    n5(1,:) = std.*randn(1,size) + std.*1i*randn(1,size); %complex noise
    h5 = relay_fading*(pyth).*randn(1,size)+ relay_fading*(pyth).*1i*randn(1,size);
       
    h_mag_yo(1,:) = (abs(h4(1,:)).^2) + (abs(h5(1,:)).^2); %(alpha0^2 + alpha1^2)
   
    %
    %
    
    count = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    count5 = 0;
    count6 = 0;
    count13 = 0;
    count14 = 0;
    count15= 0;
    count16 = 0;
    count17 = 0;
    count18 = 0;
    count19 = 0;
    count20 = 0;
    
    for i = 1:size
          
        if (real(s0_rx(1,i)) > 0) %bits 1 & 3 on Re axis
            output(1,i) = 1;        %check faded real
            if (real(s0_rx(1,i)) > 2)
                output(3,i) = 1; 
            else
                output(3,i) = -1; 
            end
        else
            output(1,i) = -1;
            if (real(s0_rx(1,i)) > -2)
                output(3,i) = 1; 
            else
                output(3,i) = -1; 
            end
        end
         if (imag(s0_rx(1,i)) > 0)     %check faded imaginary
            output(2,i) = 1;
            if (imag(s0_rx(1,i)) > 2)
                output(4,i) = 1;
            else
                output(4,i) = -1;
            end
        else
            output(2,i) = -1;
            if (imag(s0_rx(1,i)) > -2)
                output(4,i) = 1;
            else
                output(4,i) = -1;
            end
         end
        
       
        
         
        if (data1(1,i) == output(1,i))
            count = count + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output(2,i))
            count = count + 1;              %count will get b-e-r
        end
        if (data1(3,i) == output(3,i))
            count = count + 1;              %count will get b-e-r
        end
        if (data1(4,i) == output(4,i))
            count = count + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output(1,i)) | (data1(2,i) ~= output(2,i)) | (data1(3,i) ~= output(3,i)) | (data1(4,i) ~= output(4,i)) 
            count2 = count2 + 1;            %count 2 will get s-e-r
        end
          
          
          
         %checking non-faded stuff
         if (real(s1_rx(1,i)) > 0)
             output2(1,i) = 1;        %check faded real
            if (real(s1_rx(1,i)) > 2)
                output2(3,i) = 1; 
            else
                output2(3,i) = -1; 
            end
        else
            output2(1,i) = -1;
            if (real(s1_rx(1,i)) > -2)
                output2(3,i) = 1; 
            else
                output2(3,i) = -1; 
            end
        end
         if (imag(s1_rx(1,i)) > 0)     %check faded imaginary
            output2(2,i) = 1;
            if (imag(s1_rx(1,i)) > 2)
               output2(4,i) = 1; 
            else
               output2(4,i) = -1; 
            end
        else
            output2(2,i) = -1;
             if (imag(s1_rx(1,i)) > -2)
               output2(4,i) = 1; 
            else
               output2(4,i) = -1; 
            end
         end
        
         
         
        if (data2(1,i) == output2(1,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == output2(2,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == output2(3,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == output2(4,i))
            count3 = count3 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= output2(1,i)) | (data2(2,i) ~= output2(2,i)) | (data2(3,i) ~= output2(3,i)) | (data2(4,i) ~= output2(4,i))
            count4 = count4 + 1;            %count 4 will get s-e-r
        end
        
        %no diversity section
        if (real(Rx(1,i)) > 0)
            output3(1,i) = 1;        %check faded real
        else
            output3(1,i) = -1;
        end
         if (imag(Rx(1,i)) > 0)     %check faded imaginary
            output3(2,i) = 1;
        else
            output3(2,i) = -1;
        end
        if (data1(1,i) == output3(1,i))
            count5 = count5 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output3(2,i))
            count5 = count5 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output3(1,i)) | (data1(2,i) ~= output3(2,i))
            count6 = count6 + 1;            %count 2 will get s-e-r
        end
          
        
        %
        %
        % Relay node decode and re-transmit.
        %
        %
        
         if (real(s0_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
                 
            if (real(s0_rx_Relay(1,i)) > 2)
                Re_Tx1(1,i) = 3;
            else
                Re_Tx1(1,i) = 1; 
            end
        else
            
            if (real(s0_rx_Relay(1,i)) > -2)
                Re_Tx1(1,i) = -1;
            else
                Re_Tx1(1,i) = -3;
            end
        end
         if (imag(s0_rx_Relay(1,i)) > 0)     %check faded imaginary
            
            if (imag(s0_rx_Relay(1,i)) > 2)
                Re_Tx1(2,i) = 3;
            else
                Re_Tx1(2,i) = 1;
            end
        else
            
            if (imag(s0_rx_Relay(1,i)) > -2)
                Re_Tx1(2,i) = -1;
            else
                Re_Tx1(2,i) = -3;
            end
         end
         
         
          if (real(s1_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
                 
            if (real(s1_rx_Relay(1,i)) > 2)
                Re_Tx2(1,i) = 3;
            else
                Re_Tx2(1,i) = 1; 
            end
        else
            
            if (real(s1_rx_Relay(1,i)) > -2)
                Re_Tx2(1,i) = -1;
            else
                Re_Tx2(1,i) = -3;
            end
        end
         if (imag(s1_rx_Relay(1,i)) > 0)     %check faded imaginary
            
            if (imag(s1_rx_Relay(1,i)) > 2)
                Re_Tx2(2,i) = 3;
            else
                Re_Tx2(2,i) = 1;
            end
        else
            
            if (imag(s1_rx_Relay(1,i)) > -2)
                Re_Tx2(2,i) = -1;
            else
                Re_Tx2(2,i) = -3;
            end
         end
         
         
         %
         % MIDPOINT CALCS
         %
         %
          if (real(s0_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
            output7(1,i) = 1;        %check faded real
            if (real(s0_rx_Relay(1,i)) > 2)
                output7(3,i) = 1; 
            else
                output7(3,i) = -1; 
            end
        else
            output7(1,i) = -1;
            if (real(s0_rx_Relay(1,i)) > -2)
                output7(3,i) = 1; 
            else
                output7(3,i) = -1; 
            end
        end
         if (imag(s0_rx_Relay(1,i)) > 0)     %check faded imaginary
            output7(2,i) = 1;
            if (imag(s0_rx_Relay(1,i)) > 2)
                output7(4,i) = 1;
            else
                output7(4,i) = -1;
            end
        else
            output7(2,i) = -1;
            if (imag(s0_rx_Relay(1,i)) > -2)
                output7(4,i) = 1;
            else
                output7(4,i) = -1;
            end
         end
        
       
        
         
        if (data1(1,i) == output7(1,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output7(2,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        if (data1(3,i) == output7(3,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        if (data1(4,i) == output7(4,i))
            count13 = count13 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output7(1,i)) | (data1(2,i) ~= output7(2,i)) | (data1(3,i) ~= output7(3,i)) | (data1(4,i) ~= output7(4,i)) 
            count14 = count14 + 1;            %count 2 will get s-e-r
        end
          
        
        if (real(s1_rx_Relay(1,i)) > 0) %bits 1 & 3 on Re axis
            output8(1,i) = 1;        %check faded real
            if (real(s1_rx_Relay(1,i)) > 2)
                output8(3,i) = 1; 
            else
                output8(3,i) = -1; 
            end
        else
            output8(1,i) = -1;
            if (real(s1_rx_Relay(1,i)) > -2)
                output8(3,i) = 1; 
            else
                output8(3,i) = -1; 
            end
        end
         if (imag(s1_rx_Relay(1,i)) > 0)     %check faded imaginary
            output8(2,i) = 1;
            if (imag(s1_rx_Relay(1,i)) > 2)
                output8(4,i) = 1;
            else
                output8(4,i) = -1;
            end
        else
            output8(2,i) = -1;
            if (imag(s1_rx_Relay(1,i)) > -2)
                output8(4,i) = 1;
            else
                output8(4,i) = -1;
            end
         end
        
       
        
         
        if (data2(1,i) == output8(1,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == output8(2,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == output8(3,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == output8(4,i))
            count15 = count15 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= output8(1,i)) | (data2(2,i) ~= output8(2,i)) | (data2(3,i) ~= output8(3,i)) | (data2(4,i) ~= output8(4,i)) 
            count16 = count16 + 1;            %count 2 will get s-e-r
        end
        
    end
%     s4(1,:) = complex(Re_Tx1(1,:), Re_Tx1(2,:));    %Transmitter is complex of dataX and dataY
%     s5(1,:) = complex(Re_Tx2(1,:), Re_Tx2(2,:)); 

    s4(1,:) = complex(Tx1(1,:), Tx1(2,:));    %Transmitter is complex of dataX and dataY
    s5(1,:) = complex(Tx2(1,:), Tx2(2,:)); 
    
     r4(1,:) = h4(1,:).*s4(1,:) + h5(1,:).*s5(1,:) + n4(1,:);  %alamouti eqns
    r5(1,:) = -h4(1,:).*conj(s5(1,:)) + h5(1,:).*conj(s4(1,:)) + n5(1,:);
    
    s0_tilda_4(1,:) = conj(h4(1,:)).*r4(1,:) + h5(1,:).*conj(r5(1,:)); %alamouti eqns
    s1_tilda_5(1,:) = conj(h5(1,:)).*r4(1,:) - h4(1,:).*conj(r5(1,:));
    Rx1(1,:) = s0_tilda_4(1,:)./h_mag_yo(1,:);
    Rx2(1,:) = s1_tilda_5(1,:)./h_mag_yo(1,:);
    
   
    count7 = 0;
    count8 = 0;
    count9 = 0;
    count10 = 0;
    count11 = 0;
    count12 = 0;
    
    for i = 1:size
          
        if (real(Rx1(1,i)) > 0) %bits 1 & 3 on Re axis
            output4(1,i) = 1;        %check faded real
            if (real(Rx1(1,i)) > 2)
                output4(3,i) = 1; 
            else
                output4(3,i) = -1; 
            end
        else
            output4(1,i) = -1;
            if (real(Rx1(1,i)) > -2)
                output4(3,i) = 1; 
            else
                output4(3,i) = -1; 
            end
        end
         if (imag(Rx1(1,i)) > 0)     %check faded imaginary
            output4(2,i) = 1;
            if (imag(Rx1(1,i)) > 2)
                output4(4,i) = 1;
            else
                output4(4,i) = -1;
            end
        else
            output4(2,i) = -1;
            if (imag(Rx1(1,i)) > -2)
                output4(4,i) = 1;
            else
                output4(4,i) = -1;
            end
         end
        
       
        
         
        if (data1(1,i) == output4(1,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == output4(2,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        if (data1(3,i) == output4(3,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        if (data1(4,i) == output4(4,i))
            count7 = count7 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= output4(1,i)) | (data1(2,i) ~= output4(2,i)) | (data1(3,i) ~= output4(3,i)) | (data1(4,i) ~= output4(4,i)) 
            count8 = count8 + 1;            %count 2 will get s-e-r
        end
          
        
        if (real(Rx2(1,i)) > 0) %bits 1 & 3 on Re axis
            output5(1,i) = 1;        %check faded real
            if (real(Rx2(1,i)) > 2)
                output5(3,i) = 1; 
            else
                output5(3,i) = -1; 
            end
        else
            output5(1,i) = -1;
            if (real(Rx2(1,i)) > -2)
                output5(3,i) = 1; 
            else
                output5(3,i) = -1; 
            end
        end
         if (imag(Rx2(1,i)) > 0)     %check faded imaginary
            output5(2,i) = 1;
            if (imag(Rx2(1,i)) > 2)
                output5(4,i) = 1;
            else
                output5(4,i) = -1;
            end
        else
            output5(2,i) = -1;
            if (imag(Rx2(1,i)) > -2)
                output5(4,i) = 1;
            else
                output5(4,i) = -1;
            end
         end
        
       
        
         
        if (data2(1,i) == output5(1,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == output5(2,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == output5(3,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == output5(4,i))
            count9 = count9 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= output5(1,i)) | (data2(2,i) ~= output5(2,i)) | (data2(3,i) ~= output5(3,i)) | (data2(4,i) ~= output5(4,i)) 
            count10 = count10 + 1;            %count 2 will get s-e-r
        end
    end
     %for no diversity
         
     %COMBINING 2 SIGNALS USING EUCLIDIAN DISTANCE
        
        %output = s0_rx (s0 - source to dest)
        %output4 = Rx1 (s0 - relay to dest)
        %if their estimates are not equal - go Euclidian distance and set a
        %final output
        
        source1(1,:) = s0_rx(1,:);
        relay1(1,:) = Rx1(1,:);
        
         source2(1,:) = s1_rx(1,:);
        relay2(1,:) = Rx2(1,:);
        
       for i = 1:size
        
        if (output4(1,i) ~= output(1,i)) | (output4(2,i) ~= output(2,i)) | (output4(3,i) ~= output(3,i)) | (output4(4,i) ~= output(4,i)) 
            
            
            %WEIGHTED RELAY AT DESTINATION COMBINER
            dist_s0(1,1) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA1)^2;
            dist_s0(1,2) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA2)^2;
            dist_s0(1,3) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA3)^2;
            dist_s0(1,4) = (h_mag(1,i)^2)*abs(source1(1,i) - pointA4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointA4)^2;
            dist_s0(1,5) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB1)^2;
            dist_s0(1,6) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB2)^2;
            dist_s0(1,7) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB3)^2;
            dist_s0(1,8) = (h_mag(1,i)^2)*abs(source1(1,i) - pointB4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointB4)^2;
            dist_s0(1,9) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC1)^2;
            dist_s0(1,10) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC2)^2;
            dist_s0(1,11) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC3)^2;
            dist_s0(1,12) = (h_mag(1,i)^2)*abs(source1(1,i) - pointC4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointC4)^2;
            dist_s0(1,13) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD1)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD1)^2;
            dist_s0(1,14) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD2)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD2)^2;
            dist_s0(1,15) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD3)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD3)^2;
            dist_s0(1,16) = (h_mag(1,i)^2)*abs(source1(1,i) - pointD4)^2 + (h_mag_yo(1,i)^2)*abs(relay1(1,i) - pointD4)^2;
            
            
            %NO WEIGHTINGS AT COMBINER
%             dist_s0(1,1) = abs(source1(1,i) - pointA1)^2 + abs(relay1(1,i) - pointA1)^2;
%             dist_s0(1,2) = abs(source1(1,i) - pointA2)^2 + abs(relay1(1,i) - pointA2)^2;
%             dist_s0(1,3) = abs(source1(1,i) - pointA3)^2 + abs(relay1(1,i) - pointA3)^2;
%             dist_s0(1,4) = abs(source1(1,i) - pointA4)^2 + abs(relay1(1,i) - pointA4)^2;
%             dist_s0(1,5) = abs(source1(1,i) - pointB1)^2 + abs(relay1(1,i) - pointB1)^2;
%             dist_s0(1,6) = abs(source1(1,i) - pointB2)^2 + abs(relay1(1,i) - pointB2)^2;
%             dist_s0(1,7) = abs(source1(1,i) - pointB3)^2 + abs(relay1(1,i) - pointB3)^2;
%             dist_s0(1,8) = abs(source1(1,i) - pointB4)^2 + abs(relay1(1,i) - pointB4)^2;
%             dist_s0(1,9) = abs(source1(1,i) - pointC1)^2 + abs(relay1(1,i) - pointC1)^2;
%             dist_s0(1,10) = abs(source1(1,i) - pointC2)^2 + abs(relay1(1,i) - pointC2)^2;
%             dist_s0(1,11) = abs(source1(1,i) - pointC3)^2 + abs(relay1(1,i) - pointC3)^2;
%             dist_s0(1,12) = abs(source1(1,i) - pointC4)^2 + abs(relay1(1,i) - pointC4)^2;
%             dist_s0(1,13) = abs(source1(1,i) - pointD1)^2 + abs(relay1(1,i) - pointD1)^2;
%             dist_s0(1,14) = abs(source1(1,i) - pointD2)^2 + abs(relay1(1,i) - pointD2)^2;
%             dist_s0(1,15) = abs(source1(1,i) - pointD3)^2 + abs(relay1(1,i) - pointD3)^2;
%             dist_s0(1,16) = abs(source1(1,i) - pointD4)^2 + abs(relay1(1,i) - pointD4)^2;
            
            [mini, posn] = min(dist_s0);
            if posn == 1
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 1; %A1
            end
            if posn == 2
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 0; %A2
            end
            if posn == 3
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 1; %A3
            end
            if posn == 4
                Combo(1,i) = 1; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 0; %A4
            end
            if posn == 5
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 1;
            end
            if posn == 6
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 0;
            end
            if posn == 7
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 1;
            end
            if posn == 8
                Combo(1,i) = 1; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 0;
            end
            if posn == 9
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 1;
            end
            if posn == 10
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 1; Combo(4,i) = 0;
            end
            if posn == 11
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 1;
            end
            if posn == 12
                Combo(1,i) = 0; Combo(2,i) = 0; Combo(3,i) = 0; Combo(4,i) = 0;
            end
            if posn == 13
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 1;
            end
            if posn == 14
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 1; Combo(4,i) = 0;
            end
            if posn == 15
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 1;
            end
            if posn == 16
                Combo(1,i) = 0; Combo(2,i) = 1; Combo(3,i) = 0; Combo(4,i) = 0;
            end
            Combo(~Combo)= -1;
        else 
            Combo(1,i) = output4(1,i);
            Combo(2,i) = output4(2,i);
            Combo(3,i) = output4(3,i);
            Combo(4,i) = output4(4,i);
        end
        
        
        
        %output2 = s1_rx (s1 - source to dest)
        %output5 = Rx2 (s1 - relay to dest)
        %if their estimates are not equal - go Euclidian distance and set a
        %final output
       
        
        if (output5(1,i) ~= output2(1,i)) | (output5(2,i) ~= output2(2,i)) | (output5(3,i) ~= output2(3,i)) | (output5(4,i) ~= output2(4,i)) 
            
            %WEIGHTED RELAY AT DESTINATION COMBINER
            dist_s1(1,1) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA1)^2;
            dist_s1(1,2) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA2)^2;
            dist_s1(1,3) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA3)^2;
            dist_s1(1,4) = (h_mag(1,i)^2)*abs(source2(1,i) - pointA4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointA4)^2;
            dist_s1(1,5) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB1)^2;
            dist_s1(1,6) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB2)^2;
            dist_s1(1,7) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB3)^2;
            dist_s1(1,8) = (h_mag(1,i)^2)*abs(source2(1,i) - pointB4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointB4)^2;
            dist_s1(1,9) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC1)^2;
            dist_s1(1,10) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC2)^2;
            dist_s1(1,11) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC3)^2;
            dist_s1(1,12) = (h_mag(1,i)^2)*abs(source2(1,i) - pointC4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointC4)^2;
            dist_s1(1,13) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD1)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD1)^2;
            dist_s1(1,14) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD2)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD2)^2;
            dist_s1(1,15) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD3)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD3)^2;
            dist_s1(1,16) = (h_mag(1,i)^2)*abs(source2(1,i) - pointD4)^2 + (h_mag_yo(1,i)^2)*abs(relay2(1,i) - pointD4)^2;
 
            
            %NO WEIGHTINGS ON COMBINER
%             dist_s1(1,1) = abs(source2(1,i) - pointA1)^2 + abs(relay2(1,i) - pointA1)^2;
%             dist_s1(1,2) = abs(source2(1,i) - pointA2)^2 + abs(relay2(1,i) - pointA2)^2;
%             dist_s1(1,3) = abs(source2(1,i) - pointA3)^2 + abs(relay2(1,i) - pointA3)^2;
%             dist_s1(1,4) = abs(source2(1,i) - pointA4)^2 + abs(relay2(1,i) - pointA4)^2;
%             dist_s1(1,5) = abs(source2(1,i) - pointB1)^2 + abs(relay2(1,i) - pointB1)^2;
%             dist_s1(1,6) = abs(source2(1,i) - pointB2)^2 + abs(relay2(1,i) - pointB2)^2;
%             dist_s1(1,7) = abs(source2(1,i) - pointB3)^2 + abs(relay2(1,i) - pointB3)^2;
%             dist_s1(1,8) = abs(source2(1,i) - pointB4)^2 + abs(relay2(1,i) - pointB4)^2;
%             dist_s1(1,9) = abs(source2(1,i) - pointC1)^2 + abs(relay2(1,i) - pointC1)^2;
%             dist_s1(1,10) = abs(source2(1,i) - pointC2)^2 + abs(relay2(1,i) - pointC2)^2;
%             dist_s1(1,11) = abs(source2(1,i) - pointC3)^2 + abs(relay2(1,i) - pointC3)^2;
%             dist_s1(1,12) = abs(source2(1,i) - pointC4)^2 + abs(relay2(1,i) - pointC4)^2;
%             dist_s1(1,13) = abs(source2(1,i) - pointD1)^2 + abs(relay2(1,i) - pointD1)^2;
%             dist_s1(1,14) = abs(source2(1,i) - pointD2)^2 + abs(relay2(1,i) - pointD2)^2;
%             dist_s1(1,15) = abs(source2(1,i) - pointD3)^2 + abs(relay2(1,i) - pointD3)^2;
%             dist_s1(1,16) = abs(source2(1,i) - pointD4)^2 + abs(relay2(1,i) - pointD4)^2;
            
            
            [mini2, posn2] = min(dist_s1);
            if posn2 == 1
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 1; %A1
            end
            if posn2 == 2
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 0; %A2
            end
            if posn2 == 3
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 1; %A3
            end
            if posn2 == 4
                Combo2(1,i) = 1; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 0; %A4
            end
            if posn2 == 5
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 1;
            end
            if posn2 == 6
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 0;
            end
            if posn2 == 7
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 1;
            end
            if posn2 == 8
                Combo2(1,i) = 1; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 0;
            end
            if posn2 == 9
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 1;
            end
            if posn2 == 10
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 1; Combo2(4,i) = 0;
            end
            if posn2 == 11
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 1;
            end
            if posn2 == 12
                Combo2(1,i) = 0; Combo2(2,i) = 0; Combo2(3,i) = 0; Combo2(4,i) = 0;
            end
            if posn2 == 13
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 1;
            end
            if posn2 == 14
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 1; Combo2(4,i) = 0;
            end
            if posn2 == 15
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 1;
            end
            if posn2 == 16
                Combo2(1,i) = 0; Combo2(2,i) = 1; Combo2(3,i) = 0; Combo2(4,i) = 0;
            end
            Combo2(~Combo2)= -1;
            
            
        else 
            Combo2(1,i) = output5(1,i);
            Combo2(2,i) = output5(2,i);
            Combo2(3,i) = output5(3,i);
            Combo2(4,i) = output5(4,i);
        end
        
        if (data1(1,i) == Combo(1,i))      
            count17 = count17 + 1;              %count will get b-e-r
        end
        if (data1(2,i) == Combo(2,i))
            count17 = count17 + 1;              %count will get b-e-r
        end
        if (data1(3,i) == Combo(3,i))
            count17 = count17 + 1;              %count will get b-e-r
        end
        if (data1(4,i) == Combo(4,i))
            count17 = count17 + 1;              %count will get b-e-r
        end
        
        if (data1(1,i) ~= Combo(1,i)) | (data1(2,i) ~= Combo(2,i)) | (data1(3,i) ~= Combo(3,i)) | (data1(4,i) ~= Combo(4,i)) 
            count18 = count18 + 1;            %count 2 will get s-e-r
        end
        
        
        
        
        if (data2(1,i) == Combo2(1,i))      
            count19 = count19 + 1;              %count will get b-e-r
        end
        if (data2(2,i) == Combo2(2,i))
            count19 = count19 + 1;              %count will get b-e-r
        end
        if (data2(3,i) == Combo2(3,i))
            count19 = count19 + 1;              %count will get b-e-r
        end
        if (data2(4,i) == Combo2(4,i))
            count19 = count19 + 1;              %count will get b-e-r
        end
        
        if (data2(1,i) ~= Combo2(1,i)) | (data2(2,i) ~= Combo2(2,i)) | (data2(3,i) ~= Combo2(3,i)) | (data2(4,i) ~= Combo2(4,i)) 
            count20 = count20 + 1;            %count 2 will get s-e-r
        end
       end
     
     %ser_relay(1,z) = (count6 / size);
%      count_total_relay(1,z) = count5;
%      count_total_relay_2(1,z) = count6;

    ber_relay_s0(1,z) = ((size*4) - count7) / (size*4);
    ber_relay_s1(1,z) = ((size*4) - count9) / (size*4);
    count_total_relay_ber_s0(1,z) = count7;
    count_total_relay_ser_s0(1,z) = count8;
    count_total_relay_ber_s1(1,z) = count9;
    count_total_relay_ser_s1(1,z) = count10;
    
    
    
    ber_relay_s0_mp(1,z) = ((size*4) - count13) / (size*4);
    ber_relay_s1_mp(1,z) = ((size*4) - count15) / (size*4);
    count_total_relay_ber_s0_mp(1,z) = count13;
    count_total_relay_ser_s0_mp(1,z) = count14;
    count_total_relay_ber_s1_mp(1,z) = count15;
    count_total_relay_ser_s1_mp(1,z) = count16;
    
    
    
    %for fading - bit error-rate & symbol error-rate
    ber_s0(1,z) = ((size*4) - count) / (size*4);  
    ser_s0(1,z) = (count2 / size);
    count_total_s0(1,z) = count;   %two count vectors in case of errors
    count_total_s0_2(1,z) = count2;
    
    
    %for non-fading - bit error-rate & symbol error-rate
     ber_s1(1,z) = ((size*4) - count3) / (size*4);
     ser_s1(1,z) = (count4 / size);
     count_total_s1(1,z) = count3;
     count_total_s1_2(1,z) = count4;
     
    %for no diversity
     ber_no_diversity(1,z) = ((size*2) - count5) / (size*2);
     ser_no_diversity(1,z) = (count6 / size);
     count_total_no_diversity(1,z) = count5;
     count_total_no_diversity_2(1,z) = count6;
      
    ber_combo_s0(1,z) = ((size*4) - count17) / (size*4);
    ber_combo_s1(1,z) = ((size*4) - count19) / (size*4);
    count_total_combo_ber_s0(1,z) = count17;
    count_total_combo_ser_s0(1,z) = count18;
    count_total_combo_ber_s1(1,z) = count19;
    count_total_combo_ser_s1(1,z) = count20;
     
    z = z +1;
end

ber_combo(1,:) = (ber_combo_s0(1,:) + ber_combo_s1(1,:))/2;
ber_alamouti(1,:) = (ber_s0(1,:) + ber_s1(1,:))/2;
ber_relay(1,:) = (ber_relay_s0(1,:) + ber_relay_s1(1,:))/2;

db = [0:step:max_db];
%axis([0 15 10e-10 10e-1]);
% semilogy (db, ber_alamouti, 'r');
% title ('Alamouti Code Error Rate: QAM Relaying combiner - weighted Relay');
% xlabel ('Es/No (dB)');
% ylabel ('Error rate');
% 
% hold on;

%semilogy (db, ber_s1,'r');
% semilogy (db,ber_relay,'g');
hold on;
semilogy (db,ber_combo,'r');
legend('Relay re-transmission without diversity', 'Relay re-transmission using Alamouti diversity');
