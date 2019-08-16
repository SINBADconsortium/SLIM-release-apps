function [] = iwave_wrtrac_sim(model,label,data)

%% source/receiver/time
if isfield(model,'ysrc')
    osrc = [model.xsrc(1) model.ysrc(1)];
    dsrc = [model.xsrc(2) - model.xsrc(1) model.ysrc(2) - model.ysrc(1)];
    nsrc = [length(model.xsrc) length(model.ysrc)];
else
    osrc = [model.xsrc(1) 0];
    if length(model.xsrc)>1
        dsrc = [model.xsrc(2) - model.xsrc(1) 1];
    else
        dsrc = [1 1];
    end
    nsrc = [length(model.xsrc) 1];
end

nsrc = [1 1]; %%%%%%%%%%%%%%%%%%%%%%%%% temporary modification

if isfield(model,'yrec')
    orec = [model.xrec(1) model.yrec(1)];
    drec = [model.xrec(2) - model.xrec(1) model.yrec(2) - model.yrec(1)];
    nrec = [length(model.xrec) length(model.yrec)];
else
    orec = [model.xrec(1) 0];
    if length(model.xrec)>1
        drec = [model.xrec(2) - model.xrec(1) 1];
    else
        drec = [1 1];
    end
    nrec = [length(model.xrec) 1];
end

zsrc = model.zsrc(1);
zrec = model.zrec(1);

ot = model.t(1);
dt = model.t(2) - model.t(1);
nt = length(model.t);

%% headers
SuHeader.SegyFormatRevisionNumber=100;
SuHeader.DataSampleFormat=5;
SuHeader.Rev=GetSegyHeaderBasics;
SuHeader.ns=nt;
SuHeader.nsOrig=nt;
SuHeader.dt=dt*1e6;
SuHeader.dtOrig=dt*1e6;
SuHeader.FixedLengthTraceFlag=1;
SuHeader.NumberOfExtTextualHeaders=0;
it = 1;
segyid = fopen([label '_data.su'],'w','b');
for ks = 1:nsrc(2)
    for ls = 1:nsrc(1);
        for kr = 1:nrec(2)
            for lr = 1:nrec(1)
                SuTraceHeaders = InitSegyTraceHeader(nt,1e6*dt);
                
                SuTraceHeaders.ReceiverGroupElevation = -zrec;
                SuTraceHeaders.SourceSurfaceElevation   = -zsrc;
%                SuTraceHeaders.GroupX = orec(1) + (lr-1)*drec(1);
                SuTraceHeaders.GroupX = model.xrec(lr);
                SuTraceHeaders.GroupY = orec(2) + (kr-1)*drec(2);
                %SuTraceHeaders.SourceX = osrc(1) + (ls-1)*dsrc(1);
                %SuTraceHeaders.SourceY = osrc(2) + (ks-1)*dsrc(2);
                
                PutSegyTrace(segyid,data((it-1)*nt+1:it*nt),SuTraceHeaders,SuHeader);
                it = it + 1;
            end
        end
    end
end
fclose(segyid);
