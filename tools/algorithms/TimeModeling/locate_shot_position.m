function fwdpara_dis = locate_shot_position(fwdpara_dis)
%% Shot position after domain decomposition
% compute physical position and index
% Author : Xiangli edited by Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015
ny = length(fwdpara_dis.mmy_loc);

if ny == 1
	src_idx  = find(fwdpara_dis.slx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.slx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
					fwdpara_dis.slz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.slz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
	fwdpara_dis.slx_loc =  fwdpara_dis.slx(src_idx);
	fwdpara_dis.sly_loc =  fwdpara_dis.sly(src_idx);
	fwdpara_dis.slz_loc =  fwdpara_dis.slz(src_idx);		
else
	src_idx  = find(fwdpara_dis.slx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.slx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
					fwdpara_dis.sly>=fwdpara_dis.mmy_loc(fwdpara_dis.np_extra) & fwdpara_dis.sly<fwdpara_dis.mmy_loc(end-fwdpara_dis.np_extra) & ...
					fwdpara_dis.slz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.slz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
	fwdpara_dis.slx_loc =  fwdpara_dis.slx(src_idx);
	fwdpara_dis.sly_loc =  fwdpara_dis.sly(src_idx);
	fwdpara_dis.slz_loc =  fwdpara_dis.slz(src_idx);
end
fwdpara_dis.src_idx = src_idx;