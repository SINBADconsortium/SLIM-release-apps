classdef (Abstract) FWIData < handle
        
    methods (Abstract)
        Data = getData(obj,index);
        setData(obj,index,data);
        grid = getSourceGrid(obj);
    end
end