function output = DFd2DT(m,input,flag,model,options)
% Time domain modeling in the Born approximation by iwave.
% This is the Jacobian of Fd2DT(m,model,options)
%
% Use:
% output = DFd2DT(m,input,flag,model,options)


if isfield(model,'encod')
    isphase = model.encod;
else
    isphase = 0;
end

if isfield(model,'randseq')
    isrand = model.randseq;
else
    isrand = 0;
end

if isfield(model,'simultsrc')
    ismulti = model.simultsrc;
else
    ismulti = 0;
end

if isfield(model,'romb')
    isromb = model.romb;
else
    isromb = 0;
end

curdir = pwd;
%addpath(curdir);
if flag == 1
    if isfield(options,'label')
        label = [options.label,'_lin'];
    else
        label = 'lin_results_temp';
    end
else
    if isfield(options,'label')
        label = [options.label,'_adj'];
    else
        label = 'adj_results_temp';
    end
end

model.label = label;
if ~exist(label,'dir')
    mkdir(label);
end


ott = model.o;
dtt = model.d;
ntt = model.n;
model.xsrc = model.xsrc + model.d(2);
model.zsrc = model.zsrc + model.d(1);
model.xrec = model.xrec + model.d(2);
model.zrec = model.zrec + model.d(1);
if length(model.n) > 2
    if model.n(3) > 1
        ntt = model.n+2;
    else
        ntt(1:2) = model.n(1:2) + 2;
    end
else
    ntt = model.n + 2;
end
Ntt = prod(ntt);
[zt,xt] = odn2grid(ott,dtt,ntt);

Px  = opKron(opExtension(model.n(2),1),opExtension(model.n(1),1));
Px1 = opKron(opExtension(model.n(2),1,0),opExtension(model.n(1),1,0));

% Read bulk modulus and buoyancy
bulk = reshape(Px * m(1:length(m)/2),ntt);
buoy = reshape(Px * m(1+length(m)/2:end),ntt);

if flag == 1
    if options.bulkonly == 1
        dbuoy = zeros(ntt);
        dbulk = reshape(Px1 * input, ntt);
    else
        dbulk = reshape(Px1 * input(1:length(m)/2),ntt);
        dbuoy = reshape(Px1 * input(1+length(m)/2:end),ntt);
    end
    
else
    ddata = input;
end

% Write bulk modules and buoyancy into an rsf file
d   = model.d;
o   = model.o;
n   = model.n;
% WriteAllData([curdir '/' label '_buoy.rsf'],buoy,ntt,dtt,ott);
% WriteAllData([curdir '/' label '_bulk.rsf'],bulk,ntt,dtt,ott);
odnwrite_iwave([curdir '/' label '_buoy.rsf'],buoy,ott,dtt,ntt,'buoyancy');
odnwrite_iwave([curdir '/' label '_bulk.rsf'],bulk,ott,dtt,ntt,'bulkmod');
fname_buoy = [curdir '/' label '_buoy.rsf'];
fname_bulk = [curdir '/' label '_bulk.rsf'];

% Modify the rsf header file
fid = fopen(fname_buoy,'a');
str = ['data_type=buoyancy\n'];
fprintf(fid,str);
fclose(fid);

fid = fopen(fname_bulk,'a');
str = ['data_type=bulkmod\n'];
fprintf(fid,str);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag == 1
    
    % WriteAllData([curdir '/' label '_dbuoy.rsf'],dbuoy,ntt,dtt,ott);
    % WriteAllData([curdir '/' label '_dbulk.rsf'],dbulk,ntt,dtt,ott);
    odnwrite([curdir '/' label '_dbuoy.rsf'],dbuoy,ott,dtt,ntt);
    odnwrite([curdir '/' label '_dbulk.rsf'],dbulk,ott,dtt,ntt);
    
    fname_dbuoy = [curdir '/' label '_dbuoy.rsf'];
    fname_dbulk = [curdir '/' label '_dbulk.rsf'];
    % Modify the rsf header file
    fid = fopen(fname_dbuoy,'a');
    str = ['data_type=dbuoyancy\n'];
    fprintf(fid,str);
    fclose(fid);
    
    fid = fopen(fname_dbulk,'a');
    str = ['data_type=dbulkmod\n'];
    fprintf(fid,str);
    fclose(fid);
    lin_params = importdata(options.linpara);
    
    
    %% Linear part
    fname_buoy_temp  = ['../' label '_buoy.rsf'];
    fname_bulk_temp  = ['../' label '_bulk.rsf'];
    fname_dbuoy_temp = ['../' label '_dbuoy.rsf'];
    fname_dbulk_temp = ['../' label '_dbulk.rsf'];
    
    if ismulti ~= 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% modify by vts for simultaneous shots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Numsrc     = size(model.srcdata,3);
        
        
        spmd
            
            workidx                  = labindex;
            modellocal               = model;
            Part                     = codistributor1d.defaultPartition(Numsrc);
            Numsrc_local             = Part(workidx);
            nlab                     = numlabs;
            
            fname_buoy_temp          = [ label '_buoy.rsf'];
            fname_bulk_temp          = [ label '_bulk.rsf'];
            fname_dbuoy_temp         = [ label '_dbuoy.rsf'];
            fname_dbulk_temp         = [ label '_dbulk.rsf'];
            
            % Repartition of the supershots between the workers
            [label_temp]             = [label '_local' num2str(workidx)];
            id_localsrc              = sum(Part(1:workidx-1))+1:sum(Part(1:workidx));
            modellocal.srcdata       = modellocal.srcdata(:,:,id_localsrc);
            lin_params_locali        = importdata(options.linpara);
            lin_params_locali{end+1} = 'partask=1;';
            lin_params_locali{end+1} = 'srctype=''array'';';
            output_local             = [];
            
            for i = 1:Numsrc_local
                
                modellocali             = modellocal;
                modellocali.srcdata     = squeeze(modellocal.srcdata(:,:,i));
                options_locali          = options;
                label_locali            = [label_temp '_' num2str(i)];
                options_locali.label    = label_locali;
                
                fname_buoy_temp_locali  = [label_locali '_buoy.rsf'];
                fname_bulk_temp_locali  = [label_locali '_bulk.rsf'];
                fname_dbuoy_temp_locali = [label_locali '_dbuoy.rsf'];
                fname_dbulk_temp_locali = [label_locali '_dbulk.rsf'];
                
                
                [status,result]         =linwrapper(fname_buoy_temp,fname_bulk_temp,fname_dbulk_temp,fname_dbuoy_temp,modellocali,lin_params_locali,options_locali);
                assert(status==0,'iwave++ fails');
                data_label_locali       = [label_locali '_lin_data.su'];
                output_locali           = ReadSuFast(data_label_locali,length(model.t),'','b');
                output_local            = [output_local;output_locali(:)];
                
            end
            
            %Part   = codistributor1d.defaultPartition(Numsrc);
            Part   = Part * length(model.t) * length(model.xrec);
            codist = codistributor1d(1, Part, [sum(Part),1]);
            output = codistributed.build(output_local,codist,'noCommunication');
            
            
        end %end spmd
        output = gather(output);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
    elseif isrand ~=0
        Numsrc     = length(model.xsrc);
        
        
        spmd
            
            workidx                  = labindex;
            modellocal               = model;
            Part                     = codistributor1d.defaultPartition(Numsrc);
            Numsrc_local             = Part(workidx);
            nlab                     = numlabs;
            
            fname_buoy_temp          = [ label '_buoy.rsf'];
            fname_bulk_temp          = [ label '_bulk.rsf'];
            fname_dbuoy_temp         = [ label '_dbuoy.rsf'];
            fname_dbulk_temp         = [ label '_dbulk.rsf'];
            
            % Repartition of the supershots between the workers
            [label_temp]             = [label '_local' num2str(workidx)];
            id_localsrc              = sum(Part(1:workidx-1))+1:sum(Part(1:workidx));
            modellocal.xsrc          = modellocal.xsrc(id_localsrc);
            lin_params_locali        = importdata(options.linpara);
            lin_params_locali{end+1} = 'partask=1;';
            lin_params_locali{end+1} = 'srctype=''point'';';
            output_local             = [];
            
            for i = 1:Numsrc_local
                
                modellocali             = modellocal;
                modellocali.xsrc        = squeeze(modellocal.xsrc(i));
                options_locali          = options;
                label_locali            = [label_temp '_' num2str(i)];
                options_locali.label    = label_locali;
                
                fname_buoy_temp_locali  = [label_locali '_buoy.rsf'];
                fname_bulk_temp_locali  = [label_locali '_bulk.rsf'];
                fname_dbuoy_temp_locali = [label_locali '_dbuoy.rsf'];
                fname_dbulk_temp_locali = [label_locali '_dbulk.rsf'];
                
                
                [status,result]         =linwrapper(fname_buoy_temp,fname_bulk_temp,fname_dbulk_temp,fname_dbuoy_temp,modellocali,lin_params_locali,options_locali);
                assert(status==0,'iwave++ fails');
                data_label_locali       = [label_locali '_lin_data.su'];
                output_locali           = ReadSuFast(data_label_locali,length(model.t),'','b');
                output_local            = [output_local;output_locali(:)];
                
            end
            
            %Part   = codistributor1d.defaultPartition(Numsrc);
            Part   = Part * length(model.t) * length(model.xrec);
            codist = codistributor1d(1, Part, [sum(Part),1]);
            output = codistributed.build(output_local,codist,'noCommunication');
            
            
        end %end spmd
        output = gather(output);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
    elseif isphase ~= 0  || isromb ~= 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% modify by vts for phase encoded shots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Numsrc     = options.supershots;
        
        
        spmd
            
            workidx                  = labindex;
            modellocal               = model;
            Part                     = codistributor1d.defaultPartition(Numsrc);
            Numsrc_local             = Part(workidx);
            nlab                     = numlabs;
            
            fname_buoy_temp          = [ label '_buoy.rsf'];
            fname_bulk_temp          = [label '_bulk.rsf'];
            fname_dbuoy_temp         = [ label '_dbuoy.rsf'];
            fname_dbulk_temp         = [ label '_dbulk.rsf'];
            
            % Repartition of the supershots between the workers
            [label_temp]             = [label '_local' num2str(workidx)];
            if isphase ~= 0
                modellocal.phase         = model.phase(:,:,sum(Part(1:workidx-1))+1:sum(Part(1:workidx)));
            else
                modellocal.R             = model.R(sum(Part(1:workidx-1))+1:sum(Part(1:workidx)),:);
            end
            lin_params_locali        = importdata(options.linpara);
            lin_params_locali{end+1} = 'partask=1;';
            output_local             = [];
            
            for i = 1:Numsrc_local
                
                modellocali             = modellocal;
                if isphase ~= 0
                    modellocali.phase       = modellocal.phase(:,:,i);
                else
                    modellocali.R           = modellocal.R(i,:);
                end
                options_locali          = options;
                label_locali            = [label_temp '_' num2str(i)];
                options_locali.label    = label_locali;
                
                fname_buoy_temp_locali  = [label_locali '_buoy.rsf'];
                fname_bulk_temp_locali  = [label_locali '_bulk.rsf'];
                fname_dbuoy_temp_locali = [label_locali '_dbuoy.rsf'];
                fname_dbulk_temp_locali = [label_locali '_dbulk.rsf'];
                
                
                [status,result]         =linwrapper(fname_buoy_temp,fname_bulk_temp,fname_dbulk_temp,fname_dbuoy_temp,modellocali,lin_params_locali,options_locali);
                assert(status==0,'iwave++ fails');
                data_label_locali       = [label_locali '_lin_data.su'];
                output_locali           = ReadSuFast(data_label_locali,length(model.t),'','b');
                output_local            = [output_local;output_locali(:)];
                
            end
            
            %Part   = codistributor1d.defaultPartition(Numsrc);
            Part   = Part * length(model.t) * length(model.xrec);
            codist = codistributor1d(1, Part, [sum(Part),1]);
            output = codistributed.build(output_local,codist,'noCommunication');
            
            
        end %end spmd
        output = gather(output);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
        
        
        
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For classic use
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        data        = [label '_data.su'];
        nt          = length(model.t);
        [status,result] =linwrapper(fname_buoy_temp,fname_bulk_temp,fname_dbulk_temp,fname_dbuoy_temp,model,lin_params,options);
        assert(status==0,'iwave++ fails');
        output      = ReadSuFast(data,nt,'','b');
        output      = output(:);
        cd(curdir);
        
        if isfield(options,'delete')
            if options.delete == 1
                label;
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                % rmdir(label,'s');
                % delete([label '*']);
            end
        end
    end % end if condition for J
    
    
    
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          J'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% use su to write the ddata;
    if isrand == 0 && ismulti == 0 && isphase == 0 && isromb == 0
        iwave_wrtrac(model,label,ddata);
        datafile = [pwd '/' label '_data.su'];
        adj_params = importdata(options.adjpara);
        
        %% Adjoint part
        fname_buoy_temp = ['../' label '_buoy.rsf'];
        fname_bulk_temp = ['../' label '_bulk.rsf'];
        datafile_temp   = ['../' label '_data.su'];
        [status,result] = adj_wrapper(fname_buoy_temp,fname_bulk_temp,model,datafile_temp,adj_params,options);
        assert(status==0,'iwave++ fails');
        %fid             = fopen([label '_cammbulk.rsf@']);
        %dbulk           = fread(fid,Ntt,'single','ieee-le');
        cammbulkfile    = [label '_cammbulk.rsf'];
        clean_rsf_header(cammbulkfile);
        dbulk           = rsf_read_all(cammbulkfile);
        dbulk           = Px1' * dbulk(:);
        
        %fclose(fid);
        cammbuoyfile = [label '_cammbuoy.rsf'];
        if exist(cammbuoyfile,'file')
            % fid             = fopen([label '_cammbuoy.rsf@']);
            %              dbuoy           = fread(fid,Ntt,'single','ieee-le');
            % fclose(fid);
            clean_rsf_header(cammbuoyfile);
            dbuoy    = rsf_read_all(cammbuoyfile);
            dbuoy    = Px1' * dbuoy(:);
            
        else
            dbuoy = zeros(prod(model.n),1);
        end
        if options.bulkonly == 1
            output          = dbulk(:);
        else
            output          = [dbulk(:);dbuoy(:)];
        end
        cd(curdir);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmdir(label,'s');
                delete([label '*']);
            end
        end
        
    elseif isrand ~= 0 && ismulti == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modified by vts for random selection of shots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        adj_params = importdata(options.adjpara);
        Numsrc     = length(model.xsrc)*length(model.zsrc);
        ddata      = reshape(ddata,length(ddata)/Numsrc,Numsrc);
        ddata      = distributed(ddata);
        
        
        spmd
            
            workidx      = labindex;
            modellocal   = model;
            Part         = codistributor1d.defaultPartition(Numsrc);
            Numsrc_local = Part(workidx);
            dbulk_local  = zeros(prod(model.n),1);
            dbuoy_local  = dbulk_local;
            
            % Repartition of the sources between the workers
            [label_temp]             = [label '_local' num2str(workidx)];
            modellocal.xsrc          = model.xsrc(sum(Part(1:workidx-1))+1:sum(Part(1:workidx)));
            ddata_local              = getLocalPart(ddata);
            adj_params_locali        = importdata(options.adjpara);
            adj_params_locali{end+1} = 'partask=1;';
            adj_params_locali{end+1} = 'srctype=''point'';';
            
            for i = 1:length(modellocal.xsrc)
                modellocali            = modellocal;
                modellocali.xsrc       = modellocal.xsrc(i);
                options_locali         = options;
                label_locali           = [label_temp '_' num2str(i)]; % The adjpara.txt should be modified
                options_locali.label   = label_locali;
                ddata_locali           = ddata_local(:,i);
                
                iwave_wrtrac(modellocali,label_locali,ddata_locali);
                datafile_locali        = [pwd '/' label_locali '_data.su'];
                
                %% Adjoint part
                datafile_temp_locali   = [label_locali '_data.su'];
                
                [status,result]        = adj_wrapper(fname_buoy,fname_bulk,modellocali,datafile_temp_locali,adj_params_locali,options_locali);
                assert(status==0,'iwave++ fails');
                
                cammbulkfile_locali    = [label_locali '_adj_cammbulk.rsf'];
                clean_rsf_header(cammbulkfile_locali);
                dbulk_locali           = rsf_read_all(cammbulkfile_locali);
                dbulk_locali           = Px1' * dbulk_locali(:);
                dbulk_local            = dbulk_local + dbulk_locali;
                
                cammbuoyfile_locali = [label_locali '_adj_cammbuoy.rsf'];
                if exist(cammbuoyfile_locali,'file')
                    clean_rsf_header(cammbuoyfile_locali);
                    dbuoy_locali           = rsf_read_all(cammbuoyfile_locali);
                    dbuoy_locali           = Px1' * dbuoy_locali(:);
                    dbuoy_local            = dbuoy_local + dbuoy_locali;
                else
                    dbuoy                  = zeros(prod(modellocali.n),1);
                end
            end
            
            
            dbulk = pSPOT.utils.global_sum(dbulk_local);
            dbuoy = pSPOT.utils.global_sum(dbuoy_local);
            
        end %end spmd
        
        dbulk = dbulk{1};
        dbuoy = dbuoy{1};
        
        if options.bulkonly == 1
            output = dbulk(:);
        else
            output = [dbulk(:);dbuoy(:)];
        end
        
        cd(curdir);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
    elseif isphase ~= 0 || isromb ~=0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% modify by vts for phase encoded shots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        adj_params = importdata(options.adjpara);
        Numsrc     = options.supershots;
        ddata      = reshape(ddata,length(ddata)/Numsrc,Numsrc);
        ddata      = distributed(ddata);
        
        spmd
            
            workidx                  = labindex;
            modellocal               = model;
            Part                     = codistributor1d.defaultPartition(Numsrc);
            Numsrc_local             = Part(workidx);
            dbulk_local              = zeros(prod(model.n),1);
            dbuoy_local              = dbulk_local;
            
            
            % Repartition of the supershots between the workers
            [label_temp]             = [label '_local' num2str(workidx)];
            if isphase ~= 0
                modellocal.phase         = model.phase(:,:,sum(Part(1:workidx-1))+1:sum(Part(1:workidx)));
            else
                modellocal.R             = model.R(sum(Part(1:workidx-1))+1:sum(Part(1:workidx)),:);
            end
            ddata_local              = getLocalPart(ddata);
            adj_params_locali        = importdata(options.adjpara);
            adj_params_locali{end+1} = 'partask=1;';
            
            
            
            for i = 1:Numsrc_local
                modellocali            = modellocal;
                if isphase ~= 0
                    modellocali.phase      = modellocal.phase(:,:,i);
                else
                    modellocali.R           = modellocal.R(i,:);
                end
                options_locali         = options;
                label_locali           = [label_temp '_' num2str(i)];
                options_locali.label   = label_locali;
                ddata_locali           = ddata_local(:,i);
                
                iwave_wrtrac_sim(modellocali,label_locali,ddata_locali);
                datafile_locali        = [pwd '/' label_locali '_data.su'];
                
                %% Adjoint part
                fname_buoy_temp_locali = [label_locali '_buoy.rsf'];
                fname_bulk_temp_locali = [label_locali '_bulk.rsf'];
                datafile_temp_locali   = [label_locali '_data.su'];
                [status,result] = adj_wrapper(fname_buoy,fname_bulk,modellocali,datafile_temp_locali,adj_params_locali,options_locali);
                assert(status==0,'iwave++ fails');
                cammbulkfile_locali    = [label_locali '_adj_cammbulk.rsf'];
                clean_rsf_header(cammbulkfile_locali);
                dbulk_locali           = rsf_read_all(cammbulkfile_locali);
                dbulk_locali           = Px1' * dbulk_locali(:);
                dbulk_local            = dbulk_local + dbulk_locali;
                
                cammbuoyfile_locali = [label_locali '_adj_cammbuoy.rsf'];
                if exist(cammbuoyfile_locali,'file')
                    dbuoy_locali           = rsf_read_all(cammbuoyfile_locali);
                    clean_rsf_header(cammbuoyfile_locali);
                    dbuoy_locali           = Px1' * dbuoy_locali(:);
                    dbuoy_local            = dbuoy_local + dbuoy_locali;
                else
                    dbuoy                  = zeros(prod(modellocali.n),1);
                end
            end
            
            
            dbulk                    = pSPOT.utils.global_sum(dbulk_local);
            dbuoy                    = pSPOT.utils.global_sum(dbuoy_local);
            
            
        end %end spmd
        
        dbulk = dbulk{1};
        dbuoy = dbuoy{1};
        if options.bulkonly == 1
            output          = dbulk(:);
        else
            output          = [dbulk(:);dbuoy(:)];
        end
        cd(curdir);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% modify by vts for simultaneous shots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        adj_params = importdata(options.adjpara);
        Numsrc     = size(model.srcdata,3);
        ddata      = reshape(ddata,length(ddata)/Numsrc,Numsrc);
        ddata      = distributed(ddata);
        
        spmd
            
            workidx                  = labindex;
            modellocal               = model;
            Part                     = codistributor1d.defaultPartition(Numsrc);
            Numsrc_local             = Part(workidx);
            dbulk_local              = zeros(prod(model.n),1);
            dbuoy_local              = dbulk_local;
            
            
            % Repartition of the supershots between the workers
            [label_temp]             = [label '_local' num2str(workidx)];
            id_localsrc              = sum(Part(1:workidx-1))+1:sum(Part(1:workidx));
            modellocal.srcdata       = modellocal.srcdata(:,:,id_localsrc);
            ddata_local              = getLocalPart(ddata);
            adj_params_locali        = importdata(options.adjpara);
            adj_params_locali{end+1} = 'partask=1;';
            adj_params_locali{end+1} = 'srctype=''array'';';
            
            
            
            for i = 1:Numsrc_local
                modellocali            = modellocal;
                modellocali.srcdata    = squeeze(modellocal.srcdata(:,:,i));
                options_locali         = options;
                label_locali           = [label_temp '_' num2str(i)]; % The adjpara.txt should be modified
                options_locali.label   = label_locali;
                ddata_locali           = ddata_local(:,i);
                
                iwave_wrtrac_sim(modellocali,label_locali,ddata_locali);
                datafile_locali        = [pwd '/' label_locali '_data.su'];
                
                %% Adjoint part
                fname_buoy_temp_locali = [label_locali '_buoy.rsf'];
                fname_bulk_temp_locali = [label_locali '_bulk.rsf'];
                datafile_temp_locali   = [label_locali '_data.su'];
                [status,result] = adj_wrapper(fname_buoy,fname_bulk,modellocali,datafile_temp_locali,adj_params_locali,options_locali);
                assert(status==0,'iwave++ fails');
                cammbulkfile_locali    = [label_locali '_adj_cammbulk.rsf'];disp('before clean')
                clean_rsf_header(cammbulkfile_locali);disp('clean done')
                dbulk_locali           = rsf_read_all(cammbulkfile_locali);
                dbulk_locali           = Px1' * dbulk_locali(:);
                dbulk_local            = dbulk_local + dbulk_locali;
                
                cammbuoyfile_locali = [label_locali '_adj_cammbuoy.rsf'];
                if exist(cammbuoyfile_locali,'file')
                    clean_rsf_header(cammbuoyfile_locali);
                    dbuoy_locali           = rsf_read_all(cammbuoyfile_locali);
                    dbuoy_locali           = Px1' * dbuoy_locali(:);
                    dbuoy_local            = dbuoy_local + dbuoy_locali;
                else
                    dbuoy                  = zeros(prod(modellocali.n),1);
                end
            end
            
            
            dbulk                    = pSPOT.utils.global_sum(dbulk_local);
            dbuoy                    = pSPOT.utils.global_sum(dbuoy_local);
            
            
        end %end spmd
        
        dbulk = dbulk{1};
        dbuoy = dbuoy{1};
        if options.bulkonly == 1
            output          = dbulk(:);
        else
            output          = [dbulk(:);dbuoy(:)];
        end
        cd(curdir);
        
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
        
    end %end if isrand & ismulti
end %end if flag==1






