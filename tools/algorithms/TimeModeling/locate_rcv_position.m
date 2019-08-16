function fwdpara_dis = locate_rcv_position(fwdpara_dis)
%% Receiver position after domain decomposition
% compute physical position and index
% Author : XiangLi edited by Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015
ny = length(fwdpara_dis.mmy_loc);
if ny == 1
	rec_idx  = find(fwdpara_dis.rlx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
					fwdpara_dis.rlz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
					
	fwdpara_dis.rlx_loc =  fwdpara_dis.rlx(rec_idx);
	fwdpara_dis.rly_loc =  fwdpara_dis.rly(rec_idx);
	fwdpara_dis.rlz_loc =  fwdpara_dis.rlz(rec_idx);
	fwdpara_dis.rec_idx =  rec_idx;
else
	rec_idx  = find(fwdpara_dis.rlx>=fwdpara_dis.mmx_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlx<fwdpara_dis.mmx_loc(end-fwdpara_dis.np_extra) & ...
					fwdpara_dis.rly>=fwdpara_dis.mmy_loc(fwdpara_dis.np_extra) & fwdpara_dis.rly<fwdpara_dis.mmy_loc(end-fwdpara_dis.np_extra) & ...
					fwdpara_dis.rlz>=fwdpara_dis.mmz_loc(fwdpara_dis.np_extra) & fwdpara_dis.rlz<fwdpara_dis.mmz_loc(end-fwdpara_dis.np_extra));
	fwdpara_dis.rlx_loc =  fwdpara_dis.rlx(rec_idx);
	fwdpara_dis.rly_loc =  fwdpara_dis.rly(rec_idx);
	fwdpara_dis.rlz_loc =  fwdpara_dis.rlz(rec_idx);
	fwdpara_dis.rec_idx =  rec_idx;
end	