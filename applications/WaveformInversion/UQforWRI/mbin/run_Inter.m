sk = 10;
dmm = SIGMA_mp(:) .* ones(205*451,1);
S1=opSmooth(205,sk);S2=opSmooth(451,sk);S=opKron(S2,S1);
ms=S*m;
ms=reshape(ms,205,451);
ms=repmat(ms(:,1),1,451);
dm=(m-ms(:))/2;
%dm = S * dmm *4;


for i = 1:length(hh)
       fprintf('%d\n', i);
       mt    = m + hh(i) * dm;
       fq(i) = fhq(mt);
       ft(i) = fhh(mt);
       fp(i) = fpp(mt-m);
end
%figure;plot(hh(:),[ft(:) fq(:)]);legend('ft','fq')
figure;plot(hh(:),[ft(:)-fp(:) fq(:)-fp(:)]);legend('ft','fq')
title(['\lambda =' num2str(lambda(1))]);
figure;imagesc(reshape(m-dm,model.n));
   
