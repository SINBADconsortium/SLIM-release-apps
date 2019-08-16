function j = jobstart(name,works)



if nargin < 2, works = 10;end

if works == 10
    torquename = 'torque3x4';
elseif work2 == 1
    torquename = 'torque1x4';
end


j = batch(name,'Profile',torquename,'Pool',works,'CaptureDiary',true);





