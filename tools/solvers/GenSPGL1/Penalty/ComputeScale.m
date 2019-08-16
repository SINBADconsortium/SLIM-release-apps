function dfmin = ComputeScale(para)

inVal = para.scaleIn;

dfmin = scnewton(@studentsScale, inVal, 1e-6,para);


%[fminans] = fminsearch(@(x)studentsDF(x, para), inVal);



%dfmin
%dfval
%count
%fminans

end