function [distance] = circdistance(angle1, angle2, diameter)

% convert degrees to rad
if max(abs(angle1)) > 2*pi
    angle1 = pi*angle1/180;
end

if max(abs(angle2)) > 2*pi
    angle2 = pi*angle2/180;
end

ndistance = 2*(1-cos(angle1-angle2)); % this formula gives a range from 0 to 4
maxdistance = diameter*pi/2; %circumference of semi-circule
distance = ndistance*maxdistance/4; %convert distance

end