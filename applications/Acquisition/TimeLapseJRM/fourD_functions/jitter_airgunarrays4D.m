function [jitacq1, jitacq2] = jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams)

%-------------------------------------------------------------------------------------------------------------------------------------------
% jitter_airgunarrays4D outputs a structure array including the jittered acquisition parameters for baseline and monitor surveys
%
% Use:
%   [jitacq1, jitacq2] = jitter_airgunarrays4D(ns, ds, dt, rndfactor, p, rseed, boatspeed, tfireint_min, tdelay, delayboat, fig, figparams)
%
% Input:
%               ns - number of sources
%               ds - source sampling interval
%               dt - time sampling interval
%        rndfactor - [rndfactor(1) rndfactor(2)]
%                  - rndfactor(1): factor used to round the firing times (e.g., if firing times need to be rounded to the 
%                    third decimal place then rndfactor(1) = 1000)
%                  - rndfactor(2): factor used to round the jittered shot positions (e.g., if positions need to be rounded to the 
%                    second decimal place then rndfactor(2) = 100)
%                p - subsampling factor [NOTE: for data conventionally acquired at a source sampling of 50.0m 
%                    (i.e., airguns fire every 20.0 s), p = 2, 4, etc., when input data has ds = 25.0 m, 12.5 m, etc.]
%            rseed - random seed [NOTE: set four seeds, two for baseline survey and two for monitor survey.
%                    Also, random seed should be chosen in a way such that the jittered shot numbers lie in the interval [1,ns].]
%        boatspeed - speed of the boat (in meters/second)
%     tfireint_min - minimal interval between jittered firing times (in seconds), which cannot be violated for a
%                    pragmatic acquisition scenario (default is 10.0 s)
%           tdelay - time delay between airgun arrays on a boat (default is 10.0 s)
%        delayboat - time delay before the second boat starts firing (in seconds) [NOTE: when number of boats = 1, delayboat = 0]
%              fig - plot the acquisition scenario ('yes' or 'no')
%        figparams - parameters to set figure properties, such as, font size, font weight, etc.

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: April, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%-------------------------------------------------------------------------------------------------------------------------------------------


%----- ACQUISITION 1 (BASELINE SURVEY) -----%

% Jittered firing times and source positions for airgun array 1 on boat 1
rng(rseed(1));
dtfirejitb1arr1_acq1 = tfireint_min + rand(1,round(ns/p))*(2*tfireint_min);
tfirejitb1arr1_acq1 = cumsum(dtfirejitb1arr1_acq1); 
tfirejitb1arr1_acq1 = tfirejitb1arr1_acq1 - tdelay(1);
tfirejitb1arr1_acq1 = round(rndfactor(1)*tfirejitb1arr1_acq1)/rndfactor(1);
sjitb1arr1_acq1 = round(rndfactor(2)*boatspeed*tfirejitb1arr1_acq1)/rndfactor(2);  
gridsjitb1arr1IND_acq1 = round(sjitb1arr1_acq1/ds);

%disp('Acquisition 1 : airgun array 1')
%disp(['Minimum interval between jittered firing times: ' num2str(min(diff(tfirejitb1arr1_acq1))) ' s'])
%disp(['Maximum interval between jittered firing times: ' num2str(max(diff(tfirejitb1arr1_acq1))) ' s'])
if min(diff(tfirejitb1arr1_acq1)) < tfireint_min
   error('The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.')
end

%disp(['First jittered shot number: ' num2str(min(gridsjitb1arr1IND_acq1))])
%disp(['Last jittered shot number: ' num2str(max(gridsjitb1arr1IND_acq1))])
if min(gridsjitb1arr1IND_acq1) < 1 || max(gridsjitb1arr1IND_acq1) > ns
   error('The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.')
end
%fprintf('\n')


% Jittered firing times and source positions for airgun array 2 on boat 1
rng(rseed(2));
dtfirejitb1arr2_acq1 = tfireint_min + rand(1,round(ns/p))*(2*tfireint_min);
tfirejitb1arr2_acq1 = cumsum(dtfirejitb1arr2_acq1); 
tfirejitb1arr2_acq1 = round(rndfactor(1)*tfirejitb1arr2_acq1)/rndfactor(1);
sjitb1arr2_acq1 = round(rndfactor(2)*boatspeed*tfirejitb1arr2_acq1)/rndfactor(2);  
gridsjitb1arr2IND_acq1 = round(sjitb1arr2_acq1/ds);

%disp('Acquisition 1 : airgun array 2')
%disp(['Minimum interval between jittered firing times: ' num2str(min(diff(tfirejitb1arr2_acq1))) ' s'])
%disp(['Maximum interval between jittered firing times: ' num2str(max(diff(tfirejitb1arr2_acq1))) ' s'])
if min(diff(tfirejitb1arr2_acq1)) < tfireint_min
   error('The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.')
end

%disp(['First jittered shot number: ' num2str(min(gridsjitb1arr2IND_acq1))])
%disp(['Last jittered shot number: ' num2str(max(gridsjitb1arr2IND_acq1))])
if min(gridsjitb1arr2IND_acq1) < 1 || max(gridsjitb1arr2IND_acq1) > ns
   error('The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.')
end
%fprintf('\n')


% Jittered firing time grid indices 
tfirejitgrid_acq1 = 0 : dt : max(tfirejitb1arr1_acq1(end),tfirejitb1arr2_acq1(end));
tfirejitb1arr1gridIND_acq1 = zeros(1,length(tfirejitb1arr1_acq1)); 
tfirejitb1arr2gridIND_acq1 = zeros(1,length(tfirejitb1arr2_acq1)); 
for j = 1:length(tfirejitb1arr1_acq1)
    ind_arr1_acq1 = find(tfirejitgrid_acq1 <= tfirejitb1arr1_acq1(j));
    tfirejitb1arr1gridIND_acq1(j) = ind_arr1_acq1(end);
    ind_arr2_acq1 = find(tfirejitgrid_acq1 <= tfirejitb1arr2_acq1(j));
    tfirejitb1arr2gridIND_acq1(j) = ind_arr2_acq1(end);
end
tshiftb1arr1_acq1 = tfirejitgrid_acq1(tfirejitb1arr1gridIND_acq1) - tfirejitb1arr1_acq1;
tshiftb1arr2_acq1 = tfirejitgrid_acq1(tfirejitb1arr2gridIND_acq1) - tfirejitb1arr2_acq1;


%----- ACQUISITION 2 (MONITOR SURVEY) -----%

% Jittered firing times and source positions for airgun array 1 on boat 1
rng(rseed(3));
dtfirejitb1arr1_acq2 = tfireint_min + rand(1,round(ns/p))*(2*tfireint_min);
tfirejitb1arr1_acq2 = cumsum(dtfirejitb1arr1_acq2); 
tfirejitb1arr1_acq2 = tfirejitb1arr1_acq2 - tdelay(2);
tfirejitb1arr1_acq2 = round(rndfactor(1)*tfirejitb1arr1_acq2)/rndfactor(1);
sjitb1arr1_acq2 = round(rndfactor(2)*boatspeed*tfirejitb1arr1_acq2)/rndfactor(2);  
gridsjitb1arr1IND_acq2 = round(sjitb1arr1_acq2/ds);

%disp('Acquisition 2 : airgun array 1')
%disp(['Minimum interval between jittered firing times: ' num2str(min(diff(tfirejitb1arr1_acq2))) ' s'])
%disp(['Maximum interval between jittered firing times: ' num2str(max(diff(tfirejitb1arr1_acq2))) ' s'])
if min(diff(tfirejitb1arr1_acq2)) < tfireint_min
   error('The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.')
end

%disp(['First jittered shot number: ' num2str(min(gridsjitb1arr1IND_acq2))])
%disp(['Last jittered shot number: ' num2str(max(gridsjitb1arr1IND_acq2))])
if min(gridsjitb1arr1IND_acq2) < 1 || max(gridsjitb1arr1IND_acq2) > ns
   error('The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.')
end
%fprintf('\n')


% Jittered firing times and source positions for airgun array 2 on boat 1
rng(rseed(4));
dtfirejitb1arr2_acq2 = tfireint_min + rand(1,round(ns/p))*(2*tfireint_min);
tfirejitb1arr2_acq2 = cumsum(dtfirejitb1arr2_acq2) - tdelay(3); 
tfirejitb1arr2_acq2 = round(rndfactor(1)*tfirejitb1arr2_acq2)/rndfactor(1);
sjitb1arr2_acq2 = round(rndfactor(2)*boatspeed*tfirejitb1arr2_acq2)/rndfactor(2);  
gridsjitb1arr2IND_acq2 = round(sjitb1arr2_acq2/ds);

%disp('Acquisition 2 : airgun array 2')
%disp(['Minimum interval between jittered firing times: ' num2str(min(diff(tfirejitb1arr2_acq2))) ' s'])
%disp(['Maximum interval between jittered firing times: ' num2str(max(diff(tfirejitb1arr2_acq2))) ' s'])
if min(diff(tfirejitb1arr2_acq2)) < tfireint_min
   error('The minimum interval between jittered firing times cannot be less than tfireint_min. Please check the changes made in the code.')
end

%disp(['First jittered shot number: ' num2str(min(gridsjitb1arr2IND_acq2))])
%disp(['Last jittered shot number: ' num2str(max(gridsjitb1arr2IND_acq2))])
if min(gridsjitb1arr2IND_acq2) < 1 || max(gridsjitb1arr2IND_acq2) > ns
   error('The jittered shot numbers should be in the interval [1,ns]. Please change the random seed.')
end
%fprintf('\n')


% Jittered firing time grid indices 
tfirejitgrid_acq2 = 0 : dt : max(tfirejitb1arr1_acq2(end),tfirejitb1arr2_acq2(end));
tfirejitb1arr1gridIND_acq2 = zeros(1,length(tfirejitb1arr1_acq2)); 
tfirejitb1arr2gridIND_acq2 = zeros(1,length(tfirejitb1arr2_acq2)); 
for j = 1:length(tfirejitb1arr1_acq2)
    ind_arr1_acq2 = find(tfirejitgrid_acq2 <= tfirejitb1arr1_acq2(j));
    tfirejitb1arr1gridIND_acq2(j) = ind_arr1_acq2(end);
    ind_arr2_acq2 = find(tfirejitgrid_acq2 <= tfirejitb1arr2_acq2(j));
    tfirejitb1arr2gridIND_acq2(j) = ind_arr2_acq2(end);
end
tshiftb1arr1_acq2 = tfirejitgrid_acq2(tfirejitb1arr1gridIND_acq2) - tfirejitb1arr1_acq2;
tshiftb1arr2_acq2 = tfirejitgrid_acq2(tfirejitb1arr2gridIND_acq2) - tfirejitb1arr2_acq2;


% Plot the acquisition scheme
if not(exist('fig', 'var')) || isempty(fig)
    fig = 'no';
end

if strcmp(fig, 'yes')
   figure
   plot(sjitb1arr1_acq1, tfirejitb1arr1_acq1, 'o', sjitb1arr2_acq1, tfirejitb1arr2_acq1, 'b*'); 
   legend('Array 1', 'Array 2'); axis ij, axis('tight'); xlabel('Source position (m)'); ylabel('Recording time (s)');

%   xlabel('Source position (m)','FontSize',figparams.fsize_label,'FontName',figparams.fontname,'FontWeight',figparams.fontweight); 
%   ylabel('Recording time (s)','FontSize',figparams.fsize_label, 'FontName',figparams.fontname,'FontWeight',figparams.fontweight);
%   set(gca,'Fontsize',figparams.fsize_axis,'FontName',figparams.fontname,'FontWeight',figparams.fontweight,'TickDir','out')

   figure
   plot(sjitb1arr1_acq2, tfirejitb1arr1_acq2, 'o', sjitb1arr2_acq2, tfirejitb1arr2_acq2, 'b*'); 
   legend('Array 1', 'Array 2'); axis ij, axis('tight'); xlabel('Source position (m)'); ylabel('Recording time (s)');
   
%   xlabel('Source position (m)','FontSize',figparams.fsize_label,'FontName',figparams.fontname,'FontWeight',figparams.fontweight);
%   ylabel('Recording time (s)','FontSize',figparams.fsize_label, 'FontName',figparams.fontname,'FontWeight',figparams.fontweight);
%   set(gca,'Fontsize',figparams.fsize_axis,'FontName',figparams.fontname,'FontWeight',figparams.fontweight,'TickDir','out')
end


% Decide output variables
jitacq1.tfirejitb1arr1        = tfirejitb1arr1_acq1;
jitacq1.sjitb1arr1            = sjitb1arr1_acq1;
jitacq1.gridsjitb1arr1IND     = gridsjitb1arr1IND_acq1;
jitacq1.tfirejitb1arr1gridIND = tfirejitb1arr1gridIND_acq1;
jitacq1.tshiftb1arr1          = tshiftb1arr1_acq1;
jitacq1.tfirejitb1arr2        = tfirejitb1arr2_acq1;
jitacq1.sjitb1arr2            = sjitb1arr2_acq1;
jitacq1.gridsjitb1arr2IND     = gridsjitb1arr2IND_acq1;
jitacq1.tfirejitb1arr2gridIND = tfirejitb1arr2gridIND_acq1;
jitacq1.tshiftb1arr2          = tshiftb1arr2_acq1;
jitacq1.tfirejitgrid          = tfirejitgrid_acq1;

jitacq2.tfirejitb1arr1        = tfirejitb1arr1_acq2;
jitacq2.sjitb1arr1            = sjitb1arr1_acq2;
jitacq2.gridsjitb1arr1IND     = gridsjitb1arr1IND_acq2;
jitacq2.tfirejitb1arr1gridIND = tfirejitb1arr1gridIND_acq2;
jitacq2.tshiftb1arr1          = tshiftb1arr1_acq2;
jitacq2.tfirejitb1arr2        = tfirejitb1arr2_acq2;
jitacq2.sjitb1arr2            = sjitb1arr2_acq2;
jitacq2.gridsjitb1arr2IND     = gridsjitb1arr2IND_acq2;
jitacq2.tfirejitb1arr2gridIND = tfirejitb1arr2gridIND_acq2;
jitacq2.tshiftb1arr2          = tshiftb1arr2_acq2;
jitacq2.tfirejitgrid          = tfirejitgrid_acq2;

end   % function end

