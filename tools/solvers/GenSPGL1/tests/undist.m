function y = undist(x)
% Returns original value if not distributed, but gathers() if it is

if ~isa(x,'distributed')
    y = x;
else
    y = gather(x);
end
