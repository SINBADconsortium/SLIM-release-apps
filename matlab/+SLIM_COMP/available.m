function available
pckgdir=fileparts(which('SLIM_COMP.available'));
pckgpardir=fileparts(pckgdir);
dlist=dir(fullfile(pckgdir,'+*'));

    list={};
    list=addfiles(pckgdir,list);

    fprintf('MATLAP path modules available in\n %s directory:\n',pckgdir);
    for i=1:length(list)
        pname=packagename(pckgpardir,list{i});
        fprintf('\t%s\n',pname);
    end
end

function list=addfiles(mydir,list)
    mlist=dir(fullfile(mydir,'*.m'));
    for i=1:length(mlist)
        if isempty(strfind(mlist(i).name,'available.m'))
            list{end+1}=fullfile(mydir,mlist(i).name);
        end
    end
    dlist=dir(fullfile(mydir,'+*'));
    for i=1:length(dlist)
        list=addfiles(fullfile(mydir,dlist(i).name),list);
    end
end

function pname=packagename(ppdir,fname)
    dummy=strrep(fname,ppdir,'');
    dummy=dummy(3:end);
    dummy=strrep(dummy,'/+','.');
    dummy=strrep(dummy,'/','.');
    dummy=dummy(1:end-2);
    pname=dummy;
end
