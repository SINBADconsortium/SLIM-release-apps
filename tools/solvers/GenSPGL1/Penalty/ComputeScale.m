function [dfmin dfval count] = ComputeScale(para)

inVal = para.scaleIn;

[dfmin, dfval, count] = scnewton(@studentsScale, inVal, 1e-6);


%[fminans] = fminsearch(@(x)studentsDF(x, para), inVal);



%dfmin
%dfval
%count
%fminans

end