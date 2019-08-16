function band=wextract(w,J,bandnumber)
%WEXTRACT Extraction of wavelet bands.
% 	band=wextract(w,J,bandnumber) extracts the band #   
% 	'bandnumber' of a wavelet transform with a depth J. The wavelet   
% 	coefficient vector has the following fine-to-coarse organization: 
% 	w=[wav1 wav2 ...  wavJ lowJ]
% 
% 	Input:
% 	w=[wav1 wav2	... wavJ lowJ] vector of wavelet coefficients
% 	J=depth of the decomposition, i.e., J wavelet bands + 1 lowpass
% 	bandnumber=band number (1<=bandnumber<=J+1)
% 	
% 	Output:
% 	band=wav{bandnumber} 	(with the convention wav{J+1}=lowJ)

M=length(w);
if bandnumber>J+1
	disp(' ')
	disp(['You are trying to access the ' num2str(bandnumber)...
	'th band of a wavelet transform that has only ' num2str(J+1) ' bands!'])
	disp(' ')
	band=[];
	return
end
if M~=2^round(log(M)/log(2))
	disp(' ')
	disp('The size of the input signal must be a power of two!')
	disp(' ')
	band=[];
	return
end

first=1;
s=M/2;
for j=1:(bandnumber-1)
	last=first+s-1;
	s=s/2;
	first=last+1;
end
last=first+s-1;

band=w(first:last);
