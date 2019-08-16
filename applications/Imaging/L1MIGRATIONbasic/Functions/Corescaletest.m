load model

C = opMeCurvelet2d(128,384,4,16,0,'ME')

S = opCoreScale(128,384,4,16,0,5,'ME')

cc = C' * marmhard(:);

cc = S * cc;

img = C * cc;


figure;imagesc(reshape(abs(img),128,384))