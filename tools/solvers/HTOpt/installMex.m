curDir = pwd;
cd([curDir '/objective']);

mex -O mxlsmisfitht.c CFLAGS="\$CFLAGS -std=c99"

cd([curDir]);

disp(['Mex file compiled successfully']);