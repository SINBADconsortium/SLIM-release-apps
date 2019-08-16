function test_suite = test_oppKron
test_suite=buildFunctionHandleTestSuite(localfunctions);
end
function test_oppKron_5D_skipped
end

%
% function test_oppKron_5D
% %%
%     warning('off','dataCon:RedistributingX');
%     dims    = 5;
%     lim     = 10;
%     DIMDIST = randi(dims);
%     for i=1:dims
%         m    = randi(lim);
%         n    = randi(lim);
%         A{i} = opGaussian(m,n);
%     end
%     % Build oppKron
%     K1 = oppKron(A{:});
%     K2 = oppKron2Lo( opKron(A{1:end-1}),A{end},1);
%     
%     % Construct N-D array x
%     xsize = cellfun(@size,K1.children,'UniformOutput',0);
%     for i = 1:length(xsize)
%         xgsize{i} = xsize{i}(2);
%     end
%     xgsize = fliplr(xgsize);
%     x = randn(xgsize{:});
%     
%     spmd
%         x = codistributed(x,codistributor1d(DIMDIST));
%     end
%     dx = dataContainer(x);
%     dx = ivec(dx);
%     
%     % fprintf('oppKron2Lo : '),tic,
%     y2 = K2*x(:); % toc
%     % fprintf('oppKron    : '),tic,
%     y1 = K1*dx; % toc
%     
%     y1 = double(unDistriCon(vec(y1)));
%     assertElementsAlmostEqual(y1,y2);
%     
%     warning('on','dataCon:RedistributingX');
%     
% end
% 
% 
% function bleh
% %%
% 
% clear all
% clc
% q  = randn(2,1,2,1,2);
% OP = opGaussian(3,2);
% w  = spot.utils.nDimsMultiply(OP,q);
% 
% %%
% q  = randn(2,1,2,1,2);
% for i = 1:2
%     for j = 1:1
%         for k = 1:2;
%             r(:,:,i,j,k) = randn(3,1);
%         end
%     end
% end
% end
% 
