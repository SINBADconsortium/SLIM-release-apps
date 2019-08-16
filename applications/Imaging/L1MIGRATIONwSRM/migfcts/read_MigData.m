function [data, source, dt, model_bg, model_true, model_pert]=read_MigData(options)
% Syntax:
% [data, source, dt, model_bg, model_true, model_pert]=read_MigData(options)
%
% Note:
% This function is exclusively called by Migration_with_SRM. The only purpose to
% separate it is to make the main driver more human-readable.
%
% Descriptions:
% Read data, sources, and velocity models for migration. There is a placeholder
% for density models too for future extensions. ALL velocities will be of unit
% "km/s" and ALL slowness will be of unit "s/km" for output.
%
% Input list:
% options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
%
% Output list:
% data: data you migrate
% source: a struct that contains three fields. Source.Q contains the point 
% 		source, source.P contains areal source P, and source.QmP contains areal
%		source Q-P, which is set to empty for now since Q is 1D and P is 3D.
% dt: sampling interval in time
% model_bg: background model. a struct that contains two fields: velocity and
%		density that is empty for now.
% model_true: true model. a struct that contains two fields: velocity and
%		density that is empty for now. When input field "options.vel_true" is
%		empty, the velocity field is empty too.
% model_pert: model perturbation. a struct that contains two fields: slowness2
%		that means slowness squares and density that is empty for now. When
%		input field "options.model_pert" is empty, the slowness2 field is empty
%		too.
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

test_model_flag = options.test_model;
linear_flag = options.make_linear_data;
model_flag = options.model_data;
separate_flag = options.separate_input;
input_data_name = options.input_data;
data_variable_name = options.data_variable_name;
input_point_source_name = options.input_point_source;
Q_variable_name = options.Q_variable_name;
input_areal_source_name = options.input_areal_source;
P_variable_name = options.P_variable_name;
solmat_name = options.EPSI_solmat;
data_type = options.data_type;
source_type =options.source_type;
velmodel_name = options.vel_bg;
truevelmodel_name = options.vel_true;
model_pert_name = options.model_pert;
nt = options.nt;
endian = options.input_endian;
nrec = length(options.x_rcv_grid);
nshot = length(options.x_src);
nz = options.nz;
nx = options.nx;
unit = options.unit_vel;

disp('===Now start checking/reading input files.')
%% read data and source files
if not(separate_flag)
    % reading data from EPSI solution file
    % check input
    disp('Separate input flag is off.')
    disp('Will read data from EPSI solution file.')
    if isempty(solmat_name)
        error('EPSI solution file is not designated.')
    end
    if not(isempty(input_data_name))
        disp(['Conflict found: You have also designated separate data ' ...
              'input. It will be ignored.'])
    end
    if not(isempty(input_point_source_name))
        disp(['Conflict found: You have also designated separate point source ' ...
              'input. It will be ignored.'])
    end
    if not(isempty(input_areal_source_name))
        disp(['Conflict found: You have also designated separate areal source ' ...
              'input. It will be ignored.'])
    end
    % load data
    if linear_flag || test_model_flag || model_flag
        data = [];
    else
        if isempty(data_type)
            error('Fatal: data type unknown.')
        end
        switch data_type
          case 1
            load(solmat_name,'primaryIR')
            data = primaryIR;
            clear primaryIR
          case 2
            load(solmat_name,'primary')
            data = primary;
            clear primary
          case 3
            load(solmat_name,'multiple')
            data = multiple;
            clear multiple
          case 4
            load(solmat_name,'primary')
            load(solmat_name,'multiple')
            data = primary + multiple;
            clear primary multiple
        end
    end
    
    % load source
    % switch is used for more efficient memory usage by selectively
    % loading necessary source files
    if isempty(source_type)
        error('Fatal: source type not designated.')
    end
    switch source_type
      case 1
        source.I = 1;
        source.Q = [];
        source.P = [];
      case 2
        source.I = [];
        load(solmat_name,'q_est')
        [nt_conv, nb_q] = size(q_est);
        nt_Q = nt_conv/2;
        source.Q = q_est(1:nt_Q,nb_q)+q_est(nt_Q+1:nt_conv,nb_q);
        clear q_est
        source.P = [];
      case 3
        source.I = [];
        load(solmat_name,'initial_data')
        source.P = initial_data;
        clear initial_data;
        source.Q = [];
      case 4
        source.I = [];
        load(solmat_name,'q_est')
        [nt_conv, nb_q] = size(q_est);
        nt_Q = nt_conv/2;
        source.Q = q_est(1:nt_Q,nb_q)+q_est(nt_Q+1:nt_conv,nb_q);
        load(solmat_name,'initial_data')
        source.P = initial_data;
        clear q_est initial_data
    end

    % load dt
    load(solmat_name, 'dt')
    if ~exist('dt','var')
        dt = [];
    end
else
    % read separate input files
    disp('Separate input flag is on.');
    disp('Will read data from separate files.')
    if not(isempty(solmat_name))
        disp(['Conflict found: You have also designated EPSI ' ...
              'solution file as input. It will be ignored.'])
    end

    % load input data
    if linear_flag || test_model_flag || model_flag
        data = [];
    else
        if not(isempty(input_data_name))
            disp('Separate data input found.')
            disp(['Note: data_variable_name needs to be designated if ' ...
                  'input is matlab file.'])
        else
            error('Fatal: separate data input NOT designated.')
        end
        data = readDataType(input_data_name, data_variable_name, endian, nt);
        if numel(data) ~= nt*nrec*nshot
            error('Input data inconsistent with source/receiver geometry.')
        end
        data = reshape(data, nt, nrec, nshot);
    end
    % load source file
    % switch for more efficient memory usage by selectively
    % loading necessary source files
    if isempty(source_type)
        error('Fatal: source type not designated.')
    end
    switch source_type
      case 1
        source.I = 1;
        source.Q = [];
        source.P = [];
      case 2
        if not(isempty(input_point_source_name))
            disp('Separate point source input found.')
            disp(['Note: Q_variable_name needs to be designated if ' ...
                  'input is matlab file.'])
        end
        source.I = [];
        source.Q = readDataType(input_point_source_name, Q_variable_name, endian, nt);
        if (size(source.Q,1) == 2*nt) 
            source.Q = source.Q(1:nt,end)+source.Q(nt+1:2*nt,end);
        end
        source.P = [];
      case 3
        if not(isempty(input_areal_source_name))
            disp('Separate areal source input found.')
            disp(['Note: P_variable_name needs to be designated if ' ...
                  'input is matlab file.'])
        end
        source.I = [];
        source.P = readDataType(input_areal_source_name, P_variable_name, endian, nt);
        if numel(source.P) ~= nt*nrec*nshot
            error('Input areal source inconsistent with source/receiver geometry.')
        end
        source.P = reshape(source.P, nt, nrec, nshot);
        source.Q = [];
      case 4
        if not(isempty(input_point_source_name))
            disp('Separate point source input found.')
            disp(['Note: Q_variable_name needs to be designated if ' ...
                  'input is matlab file.'])
        end
        if not(isempty(input_areal_source_name))
            disp('Separate areal source input found.')
            disp(['Note: P_variable_name needs to be designated if ' ...
                  'input is matlab file.'])
        end
        source.I = [];
        source.Q = readDataType(input_point_source_name, Q_variable_name, endian, nt);
        if (size(source.Q,1) == 2*nt) 
            source.Q = source.Q(1:nt,end)+source.Q(nt+1:2*nt,end);
        end
        source.P = readDataType(input_areal_source_name, P_variable_name, endian, nt);
        if numel(source.P) ~= nt*nrec*nshot
            error('Input areal source inconsistent with source/receiver geometry.')
        end
        source.P = reshape(source.P, nt, nrec, nshot);
    end
    % load dt
    dt = [];
end
%% load models: velocity and density(kept here for potential future use)

model_bg.vel = readDataType(velmodel_name, 'vel_bg', endian, nz);
model_bg.vel = reshape(model_bg.vel, nz, nx);
model_bg.den = [];

% load true models if they are present
if ~isempty(truevelmodel_name)
    model_true.vel = readDataType(truevelmodel_name, 'vel_true', endian, nz);
    model_true.vel = reshape(model_true.vel, nz, nx);
else
    model_true.vel = [];
end
model_true.den = [];

if ~isempty(model_pert_name)
    model_pert.slowness2 = readDataType(model_pert_name, 'model_pert', endian, nz);
else
    model_pert.slowness2 = [];
end
model_pert.den = [];

% change to 'km/s'
if strcmp(unit,'m/s')
    model_bg.vel = model_bg.vel/1000;
    model_true.vel = model_true.vel/1000;
    model_pert.slowness2 = model_pert.slowness2*1e6;
else
    if ~strcmp(unit, 'km/s')
        error('Unit of velocity can only be m/s or km/s.')
    end
end

disp('Data successfully read.')