classdef FWIFreqData < FWIData
% Abstract frequency slice container that stores its own metadata
% 
% Curt Da Silva, 2016
% 
% Usage:
%  D is of type double
%   Data = FWIFreqData(D,model);
%   
%   Data = FWIFreqData(D,srcgrid,freqgrid);
% 
%  D is of type 'segyCon'
%   Data = FWIFreqData(D); 
%
% Input:
%   model    - regular grid model struct
%   srcgrid  - SourceGrid object
%   freqgrid - array of frequencies (Hz)
%   D        - data matrix (double) or data container (segyCon)
% 

    properties (SetAccess = protected)
        localData,srcgrid,freqgrid;
    end
    
    methods                         
        function Data = FWIFreqData(D,varargin)
            if isa(D,'double')            
                if isa(varargin{1},'struct')
                    model = varargin{1};
                    if isfield(model,'ysrc')
                        sg = RegularGrid(model.xsrc,model.ysrc,model.zsrc);
                        rg = RegularGrid(model.xrec,model.yrec,model.zrec);
                    else
                        sg = RegularGrid(model.zsrc,model.xsrc);
                        rg = RegularGrid(model.zrec,model.xrec);
                    end
                    srcgrid = DataGrid(sg,rg);
                    freqgrid = model.freq;
                elseif isa(varargin{1},'DataGrid')
                    srcgrid = varargin{1};
                    freqgrid = varargin{2};
                else
                    error('Unrecognized input');
                end
                Data.srcgrid = srcgrid;
                Data.localData = reshape(D,srcgrid.numrecs(),numel(D)/srcgrid.numrecs());
                Data.freqgrid = freqgrid;
            elseif isa(D,'segyCon')
                % Pull in the source/receiver locations
                % Make frequency slice
                % Store local copy of the data
                % Construct SourceGrid object
                % Once I figure out how segyCon works
                iF = varargin{1};
                hdr_labels = {'filename','mintrace','maxtrace','minsrcx','maxsrcx','minsrcy','maxsrcy','minrecx','maxrecx','minrecy','maxrecy','minnt','maxnt','mindt','maxdt'};
                hdr = C.header;
                find_str = @(x) find(cellfun(@(y) strcmp(x,y),hdr_labels));
                dt = hdr.metadata{1,find_str('mindt')};
                dt = dt*1e-6; % microsec -> sec
                nt = hdr.metadata{1,find_str('minnt')};
                t = (0:(nt-1))*dt;
                f = (0:(nt-1)/2)*1/(max(t));
                
                blkhdr_labels = {'tracenum','srcx','srcy','recx','recy'};
                find_blk_str = @(x) find(cellfun(@(y) strcmp(x,y),blkhdr_labels));
                isx = find_blk_str('srcx');
                isy = find_blk_str('srcy');
                irx = find_blk_str('recx');
                iry = find_blk_str('recy');                
                [~,~,iF] = intersect(f,iF);
                nf = length(iF);
                if nf==0, error(['Frequencies not found on grid with ' num2str((nt-1)/2) ' freqs and dt = ' num2str(1/max(t)) ]); end
                srcrec = [];
                localData = [];
                % Some metadata is improperly imported as int16s, resulting in overflow
                to_uint = @(x) reshape(cast(typecast(int16(vec(x)),'uint16'),'int32'),size(x));
                for k=1:size(C)
                    B = C.blocks(k);
                    blkhdr = B.header;
                    blk_metadata = blkhdr.metadata;        
                    if min(vec(blk_metadata(:,[isx isy irx iry]))) < 0                        
                        blk_metadata = to_uint(blk_metadata);
                    end
                    srcrec = [srcrec; [blk_metadata(:,isx),blk_metadata(:,isy),blk_metadata(:,irx),blk_metadata(:,iry) ]];            
                    B = B(1:size(B,1),1:size(B,2));
                    B = fft(B,[],1)/sqrt(size(B,1));
                    localData = [localData, B(iF,:)];
                end
                usxsy = unique(srcrec(:,1:2),'rows');
                ns = size(usxsy,1);
                recs = cell(ns,1);
                freq_slices = cell(ns,nf);
                for k=1:ns
                    I = find(srcrec(:,1)==usxsy(k,1) & srcrec(:,2)==usxsy(k,2));
                    recs{k} = [srcrec(I,3),srcrec(I,4)];    
                    for i=1:nf
                        freq_slices{i,k} = localData(i,I);
                    end
                end                                                              
            end
        end
        
        function Data = getData(obj,index)
            Data = obj.localData(:,index);
        end
        function setData(obj,index,data)
            obj.localData(:,index) = data;
        end
        function grid = getSourceGrid(obj)
            grid = obj.srcgrid;
        end
    end
end