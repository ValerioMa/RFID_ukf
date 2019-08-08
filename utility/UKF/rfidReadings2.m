function [h] = rfidReadings2(state, past_state,RFID)
%UPDATE_AGENT Summary of this function goes here
%   Detailed explanation goes here

measure_num = size(RFID,2);
particle_num = size(state,2);

% SLOW SOLUTION
h = zeros(measure_num,particle_num); 

x_p = state(1,:);
y_p = state(2,:);

x_pp = past_state(1,:);
y_pp = past_state(2,:);

for i=1:measure_num
    x_rfid = RFID(1,i);
    y_rfid = RFID(2,i);
    
    dx = x_rfid - x_p;
    dy = y_rfid - y_p;
    
    dxp = x_rfid - x_pp;
    dyp = y_rfid - y_pp;
    normaliz1 = sqrt(dx.*dx + dy.*dy);
    normaliz2 = sqrt(dxp.*dxp + dyp.*dyp);
        
    % FAST VERSION        
    h(i,:) = normaliz1 - normaliz2;        

end

end

