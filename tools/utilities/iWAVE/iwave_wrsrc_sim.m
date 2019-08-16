function [] = iwave_wrsrc_sim(model,label)
% modify by vts

if isfield(model,'randseq')
    isrand = model.randseq;
else
    isrand = 0;
end

%
% if isfield(model,'encod')
%     isphase = model.encod;
% else
%     isphase = 0;
% end
%
% if isfield(model,'romb')
%     isromb = model.romb;
% else
%     isromb = 0;
% end

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
zsrc = model.zsrc(1);

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
segyid = fopen([label '_src.su'],'w','b');
for ks = 1:nsrc(2)
    for ls = 1:nsrc(1);
        SuTraceHeaders = InitSegyTraceHeader(nt,1e6*dt);
        
        SuTraceHeaders.ReceiverGroupElevation = -zsrc;
        if isrand < 1
            SuTraceHeaders.GroupX = osrc(1) + (ls-1)*dsrc(1);
            SuTraceHeaders.GroupY = osrc(2) + (ks-1)*dsrc(2);
        else
            SuTraceHeaders.GroupX = model.xsrc(ls);
            SuTraceHeaders.GroupY = osrc(2) + (ks-1)*dsrc(2);
        end
        
        if isfield(model,'srcdata')
            srcdata = model.srcdata(:,ls);
        else
            srcdata = getrick(model.f0,model.t0,dt,nt);
        end
        
        PutSegyTrace(segyid,srcdata,SuTraceHeaders,SuHeader);
        
        % if  isphase == 1
%             PutSegyTrace(segyid,getrickPhaseEncod(model.f0,model.t0,dt,nt,model.phase(:,ls)),SuTraceHeaders,SuHeader);
%         elseif isromb == 1
%             PutSegyTrace(segyid,getrickRomb(model.f0,model.t0,dt,nt,model.R{ls}),SuTraceHeaders,SuHeader);
%         else
%             PutSegyTrace(segyid,model.weights(ls)*getrick(model.f0,model.t0,dt,nt),SuTraceHeaders,SuHeader);
%         end
        it = it + 1;
    end
end
fclose(segyid);
