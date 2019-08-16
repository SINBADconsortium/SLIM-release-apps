function [] = iwave_wrdatahdr_sim(model,label)
% iwave_wrdatahdr_sim(model,label)
%% source/receiver/time

%if ~isfield(model,'ysrc')
  % model.ysrc = 0;
%end

if ~isfield(model,'yrec')
    model.yrec = 0;
end


ot = model.t(1);
dt = model.t(2) - model.t(1);
nt = length(model.t);

[orec,drec,nrec] = grid2odn(model.xrec,model.yrec,model.zrec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(model,'ysrc')
    osrc = [model.xsrc(1) model.ysrc(1) model.zsrc(1)];
    dsrc = [model.xsrc(2) - model.xsrc(1) model.ysrc(2) - model.ysrc(1)];
    nsrc = [length(model.xsrc) length(model.ysrc)];
else
    osrc = [model.xsrc(1) 0 model.zsrc(1)];
    if length(model.xsrc)>1
        dsrc = [model.xsrc(2) - model.xsrc(1) 1];
    else
        dsrc = [1 1];
    end
    nsrc = [length(model.xsrc) 1];
end

nsrc = [1 1]; %%%%%%%%%%%%%%%%%%%%%%%%% temporary modification
dsrc = [1 1]; %%%%%%%%%%%%%%%%%%%%%%%%% temporary modification
%[osrc,dsrc,nsrc] = grid2odn(1,model.ysrc,model.zsrc); % temporary modification for a single silm. shot

if ~isfield(model,'ysrc') %%%%%%% temp form
   model.ysrc = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% headers

[Gxr,Gyr] = ndgrid(model.xrec,model.yrec); % temporary modification for a single silm. shot


fid = fopen([label '.hdr'],'w');
fprintf(fid,'%f %f %f %f\n',[Gyr(:)';Gxr(:)']); % temporary modification for a single silm. shot

fclose(fid);

%cmd1 = sprintf('a2b < %s.hdr n1 = 6 > %s_hdr.bin',label,label);
%cmd2 = sprintf('sunull nt=%d ntr=%d dt=%f > %s',nt,nrec(1)*nrec(2)*nsrc(1)*nsrc(2),dt,[label '_null.su']);
%cmd3 = sprintf('sushw < %s_null.su infile=%s_hdr.bin key=sy,sx,gy,gx | sushw key=selev a=%f | sushw key=gelev a=%f > %s_hdr.su',label,label,-osrc(3),-orec(3),label);

fid = fopen([label '_mkhdr.sh'],'w');
%fprintf(fid,'!#/bin/bash\n');
%fprintf(fid,'module unload CWP/native/42\n');
%fprintf(fid,'module unload CWP/xdr/42\n');
%fprintf(fid,'module load CWP/xdr/42\n');
fprintf(fid,'a2b < %s.hdr n1 = 6 > %s_hdr.bin\n',label,label);
fprintf(fid,'sunull nt=%d ntr=%d dt=%f > %s\n',nt,nrec(1)*nrec(2),dt,[label '_null.su']);
fprintf(fid,'sushw < %s_null.su infile=%s_hdr.bin key=gy,gx | sushw key=selev a=%f | sushw key=gelev a=%f > %s_hdr.su\n',label,label,-osrc(3),-orec(3),label);
fprintf(fid,'%s\n',['rm -f ' label '_null.su ' label '.hdr ' label '_hdr.bin']);
fclose(fid);

%system(['module load CWP/xdr/42; ' cmd1 ';' cmd2 ';' cmd3]);
%system(['module load CWP/native/42; ' cmd1 ';' cmd2 ';' cmd3]);
%system(['rm ' label '_null.su ' label '.hdr ' label '_hdr.bin']);

