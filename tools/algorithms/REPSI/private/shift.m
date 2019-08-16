function [data] = shift(indata,nshi);
%function [data] = shift(indata,nshi);

%Executes shift;
sizeData                   = size(indata);
data                       = zeros(size(indata));
if (nshi > 0)
   data(nshi+1:sizeData(1),:) = indata(1:sizeData(1)-nshi,:);

   data(1:nshi,:)             = indata(sizeData(1)-nshi+1:sizeData(1),:);
else
   negshi = -nshi;
   data(1:sizeData(1)-negshi,:) = indata(negshi+1:sizeData(1),:);
   data(sizeData(1)-negshi+1:sizeData(1),:) = indata(1:negshi,:);
end
