function [out]=proj_card(in,c)
% computes projection onto the set with cardinality(nonzeros) = c
% c: integer, keep this many nonzeros
co=c;
c=round(c);
if norm(c-co)>eps; warning('provided a non integer cardinality constraint, rounded and continued');end

[~, ind] = sort(abs(in), 'ascend');
out = in;
out(ind(1:(end-c))) = 0;

