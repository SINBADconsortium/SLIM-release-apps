function [] = iwave_wrmodel(buoy,bulk,model,label)
%iwave_wrmodel(m,rho,model,label)
if isempty(bulk)
    bulk = 0*buoy + 1;
end

odnwrite([label '_buoy.rsf'],buoy,model.o,model.d,model.n);
odnwrite([label '_bulk.rsf'],bulk,model.o,model.d,model.n);

%rsf_write_all([label '_vp.rsf'],{'out=stdout'},reshape(1./sqrt(m(:)),model.n),model.d,model.o);
%rsf_write_all([label '_rho.rsf'],{'out=stdout'},reshape(rho(:),model.n),model.d,model.o);
