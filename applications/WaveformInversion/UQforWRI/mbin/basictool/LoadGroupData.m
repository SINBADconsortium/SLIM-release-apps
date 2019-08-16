function [V NormDx] = LoadGroupData(n)
%% function to load all existed 'x_.*dat';
% [V ] = LoadGroupData(n)
% n - size of the model
i = 0;
while exist(['./x_' num2str(i) '.dat'], 'file')
        V(:,:,i+1) = loaddat(['./x_' num2str(i) '.dat'],n);
        i = i + 1;
end

for i = 1:size(V,3)-1
        NormDx(i) = norm(vec(V(:,:,i+1) - V(:,:,i)));
end
