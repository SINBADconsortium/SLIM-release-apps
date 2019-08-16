function codegen_helm3d_7p_mvp
% Matlab function to compile


% various coder options, see the documentation for details
cfg = coder.config('mex');
cfg.MATLABSourceComments = true;
cfg.EnableDebugging = false;
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = false;
cfg.EnableOpenMP = true; % won't do anything with old clang
%cfg.BuildConfiguration = 'Faster Runs';
%cfg.PostCodeGenCommand = 'buildInfo.addLinkFlags(''-fopenmp'')'; 
%cfg.CustomLibrary = '-fopenmp';

% setup type + size of arguments
% here wn is a n x 1 complex vector, where n can be anything
wn = coder.typeof(complex(double(0),double(0)),[Inf 1]);

% h, n are 3 x 1 double vector
h = coder.typeof(double(0),[3 1]);
n = coder.typeof(double(0),[3 1]);

freq = coder.typeof(complex(double(0),double(0)),1,1);

% npml is a 2 x 3 double matrix
npml = coder.typeof(double(0),[2 3]);

% x is a n x 1 complex vector, n can be anything
x = coder.typeof(complex(double(0),double(0)),[Inf 1]);

% mode is a double scalar
mode = coder.typeof(double(0),1,1);

% all of the arguments to this function
args = {wn,h,n,npml,freq,x,mode};

% run Matlab Coder
codegen('helm3d_7pt_mvp','-config',cfg,'-args',args);
end