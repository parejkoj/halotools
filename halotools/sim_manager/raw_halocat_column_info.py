# -*- coding: utf-8 -*-
"""
Collection of column info dtypes used to read halo catalog ASCII
"""
import numpy as np 

################################################################################################


################################################################################################
######## Headers copied-and-pasted directly from ASCII downloaded on the indicated date ########
######## dtypes written by hand according to the header ########

################################################################################################
##### Consuelo Rockstar header of file downloaded from 
##### http://www.slac.stanford.edu/~behroozi/Consuelo_Catalogs on 19 July, 2015
#
header_slac_consuelo_rockstar_july19_2015='scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Rs_Klypin(35) Mvir_all(36) M200b(37) M200c(38) M500c(39) M2500c(40) Xoff(41) Voff(42) Spin_Bullock(43) b_to_a(44) c_to_a(45) A[x](46) A[y](47) A[z](48) b_to_a(500c)(49) c_to_a(500c)(50) A[x](500c)(51) A[y](500c)(52) A[z](500c)(53) T/|U|(54) M_pe_Behroozi(55) M_pe_Diemer(56) Halfmass_Radius(57) Macc(58) Mpeak(59) Vacc(60) Vpeak(61) Halfmass_Scale(62) Acc_Rate_Inst(63) Acc_Rate_100Myr(64) Acc_Rate_1*Tdyn(65) Acc_Rate_2*Tdyn(66) Acc_Rate_Mpeak(67) Mpeak_Scale(68) Acc_Scale(69) First_Acc_Scale(70) First_Acc_Mvir(71) First_Acc_Vmax(72) Vmax@Mpeak(73)'
#
dtype_slac_consuelo_rockstar_july19_2015 = np.dtype([
    ('halo_scale_factor', 'f4'), 
    ('halo_id', 'i8'), 
    ('halo_scale_factor_desc', 'f4'), 
    ('halo_id_desc', 'i8'), 
    ('halo_num_prog', 'i4'), 
    ('halo_pid', 'i8'), 
    ('halo_upid', 'i8'), 
    ('halo_pid_desc', 'i8'), 
    ('halo_phantom', 'i4'),  
    ('halo_mvir_sam', 'f4'), 
    ('halo_mvir', 'f4'), 
    ('halo_rvir', 'f4'),  
    ('halo_rs', 'f4'), 
    ('halo_vrms', 'f4'), 
    ('halo_mmp', 'i4'), 
    ('halo_scale_factor_lastmm', 'f4'), 
    ('halo_vmax', 'f4'), 
    ('halo_x', 'f4'), 
    ('halo_y', 'f4'), 
    ('halo_z', 'f4'), 
    ('halo_vx', 'f4'), 
    ('halo_vy', 'f4'), 
    ('halo_vz', 'f4'), 
    ('halo_jx', 'f4'), 
    ('halo_jy', 'f4'), 
    ('halo_jz', 'f4'), 
    ('halo_spin', 'f4'), 
    ('halo_id_breadth_first', 'i8'), 
    ('halo_id_depth_first', 'i8'), 
    ('halo_id_tree_root', 'i8'), 
    ('halo_id_orig', 'i8'), 
    ('halo_snap_num', 'i4'), 
    ('halo_id_next_coprog_depthfirst', 'i8'), 
    ('halo_id_last_prog_depthfirst', 'i8'), 
    ('halo_id_last_mainleaf_depthfirst', 'i8'),
    ('halo_rs_klypin', 'f4'), 
    ('halo_mvir_all', 'f4'), 
    ('halo_m200b', 'f4'), 
    ('halo_m200c', 'f4'), 
    ('halo_m500c', 'f4'), 
    ('halo_m2500c', 'f4'), 
    ('halo_xoff', 'f4'), 
    ('halo_voff', 'f4'), 
    ('halo_spin_bullock', 'f4'), 
    ('halo_b_to_a', 'f4'), 
    ('halo_c_to_a', 'f4'), 
    ('halo_axisA_x', 'f4'), 
    ('halo_axisA_y', 'f4'), 
    ('halo_axisA_z', 'f4'), 
    ('halo_b_to_a_500c', 'f4'), 
    ('halo_c_to_a_500c', 'f4'), 
    ('halo_axisA_x_500c', 'f4'), 
    ('halo_axisA_y_500c', 'f4'), 
    ('halo_axisA_z_500c', 'f4'), 
    ('halo_t_by_u', 'f4'), 
    ('halo_mass_pe_behroozi', 'f4'), 
    ('halo_mass_pe_diemer', 'f4'), 
    ('halo_halfmass_radius', 'f4'), 
    ('halo_macc', 'f4'), 
    ('halo_mpeak', 'f4'), 
    ('halo_vacc', 'f4'), 
    ('halo_vpeak', 'f4'), 
    ('halo_halfmass_scale', 'f4'), 
    ('halo_dmvir_dt_inst', 'f4'), 
    ('halo_dmvir_dt_100myr', 'f4'), 
    ('halo_dmvir_dt_tdyn', 'f4'), 
    ('halo_dmvir_dt_2dtyn', 'f4'), 
    ('halo_dmvir_dt_mpeak', 'f4'), 
    ('halo_scale_factor_mpeak', 'f4'), 
    ('halo_scale_factor_lastacc', 'f4'), 
    ('halo_scale_factor_firstacc', 'f4'), 
    ('halo_mvir_firstacc', 'f4'), 
    ('halo_vmax_firstacc', 'f4'), 
    ('halo_vmax_mpeak', 'f4')
    ])


#
################################################################################################
##### Bolshoi BDM header of file downloaded from 
##### http://www.slac.stanford.edu/~behroozi/Bolshoi_Catalogs_BDM on 19 July, 2015
#
header_slac_bolshoi_bdm_july19_2015 = 'scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Xoff 2K/Ep-1 Rrms Axba Axca Xax Yax Zax Macc Mpeak Vacc Vpeak'
#
dtype_slac_bolshoi_bdm_july19_2015 = np.dtype([
    ('halo_scale_factor', 'f4'), 
    ('halo_id', 'i8'), 
    ('halo_scale_factor_desc', 'f4'), 
    ('halo_id_desc', 'i8'), 
    ('halo_num_prog', 'i4'), 
    ('halo_pid', 'i8'), 
    ('halo_upid', 'i8'), 
    ('halo_pid_desc', 'i8'), 
    ('halo_phantom', 'i4'),  
    ('halo_mvir_sam', 'f4'), 
    ('halo_mvir', 'f4'), 
    ('halo_rvir', 'f4'),  
    ('halo_rs', 'f4'), 
    ('halo_vrms', 'f4'), 
    ('halo_mmp', 'i4'), 
    ('halo_scale_factor_lastmm', 'f4'), 
    ('halo_vmax', 'f4'), 
    ('halo_x', 'f4'), 
    ('halo_y', 'f4'), 
    ('halo_z', 'f4'), 
    ('halo_vx', 'f4'), 
    ('halo_vy', 'f4'), 
    ('halo_vz', 'f4'), 
    ('halo_jx', 'f4'), 
    ('halo_jy', 'f4'), 
    ('halo_jz', 'f4'), 
    ('halo_spin', 'f4'), 
    ('halo_id_breadth_first', 'i8'), 
    ('halo_id_depth_first', 'i8'), 
    ('halo_id_tree_root', 'i8'), 
    ('halo_id_orig', 'i8'), 
    ('halo_snap_num', 'i4'), 
    ('halo_id_next_coprog_depthfirst', 'i8'), 
    ('halo_id_last_prog_depthfirst', 'i8'), 
    ('halo_xoff', 'f4'), 
    ('halo_2K/Ep-1', 'f4'), 
    ('halo_Rrms', 'f4'), 
    ('halo_b_to_a', 'f4'), 
    ('halo_c_to_a', 'f4'), 
    ('halo_axisA_x', 'f4'), 
    ('halo_axisA_y', 'f4'), 
    ('halo_axisA_z', 'f4'), 
    ('halo_macc', 'f4'), 
    ('halo_mpeak', 'f4'), 
    ('halo_vacc', 'f4'), 
    ('halo_vpeak', 'f4')
    ]) 


################################################################################################
##### Bolshoi rockstar header of file downloaded from 
##### http://www.slac.stanford.edu/~behroozi/Bolshoi_Catalogs/ on 19 July, 2015
#
##### Bolshoi-Planck rockstar header of file downloaded from 
# http://www.slac.stanford.edu/~behroozi/BPlanck_Hlists/ on 19 July, 2015
#
##### Multidark rockstar header of file downloaded from 
# http://slac.stanford.edu/~behroozi/MultiDark_Hlists_Rockstar on 19 July, 2015
#
##### All three headers were identical, so only the one appears below
#
header_slac_bolshoi_rockstar_july19_2015 = 'scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) pid(5) upid(6) desc_pid(7) phantom(8) sam_mvir(9) mvir(10) rvir(11) rs(12) vrms(13) mmp?(14) scale_of_last_MM(15) vmax(16) x(17) y(18) z(19) vx(20) vy(21) vz(22) Jx(23) Jy(24) Jz(25) Spin(26) Breadth_first_ID(27) Depth_first_ID(28) Tree_root_ID(29) Orig_halo_ID(30) Snap_num(31) Next_coprogenitor_depthfirst_ID(32) Last_progenitor_depthfirst_ID(33) Last_mainleaf_depthfirst_ID(34) Rs_Klypin(35) Mvir_all(36) M200b(37) M200c(38) M500c(39) M2500c(40) Xoff(41) Voff(42) Spin_Bullock(43) b_to_a(44) c_to_a(45) A[x](46) A[y](47) A[z](48) b_to_a(500c)(49) c_to_a(500c)(50) A[x](500c)(51) A[y](500c)(52) A[z](500c)(53) T/|U|(54) M_pe_Behroozi(55) M_pe_Diemer(56) Macc(57) Mpeak(58) Vacc(59) Vpeak(60) Halfmass_Scale(61) Acc_Rate_Inst(62) Acc_Rate_100Myr(63) Acc_Rate_1*Tdyn(64) Acc_Rate_2*Tdyn(65) Acc_Rate_Mpeak(66) Mpeak_Scale(67) Acc_Scale(68) First_Acc_Scale(69) First_Acc_Mvir(70) First_Acc_Vmax(71) Vmax@Mpeak(72)'
header_slac_multidark_rockstar_july19_2015 = header_slac_bolshoi_rockstar_july19_2015
header_slac_bolplanck_rockstar_july19_2015 = header_slac_bolshoi_rockstar_july19_2015
#
dtype_slac_bolshoi_rockstar_july19_2015 = np.dtype([
    ('halo_scale_factor', 'f4'), 
    ('halo_id', 'i8'), 
    ('halo_scale_factor_desc', 'f4'), 
    ('halo_id_desc', 'i8'), 
    ('halo_num_prog', 'i4'), 
    ('halo_pid', 'i8'), 
    ('halo_upid', 'i8'), 
    ('halo_pid_desc', 'i8'), 
    ('halo_phantom', 'i4'), 
    ('halo_mvir_sam', 'f4'), 
    ('halo_mvir', 'f4'), 
    ('halo_rvir', 'f4'), 
    ('halo_rs', 'f4'), 
    ('halo_vrms', 'f4'), 
    ('halo_mmp', 'i4'), 
    ('halo_scale_factor_lastmm', 'f4'), 
    ('halo_vmax', 'f4'),
    ('halo_x', 'f4'), 
    ('halo_y', 'f4'), 
    ('halo_z', 'f4'), 
    ('halo_vx', 'f4'), 
    ('halo_vy', 'f4'), 
    ('halo_vz', 'f4'), 
    ('halo_jx', 'f4'), 
    ('halo_jy', 'f4'), 
    ('halo_jz', 'f4'), 
    ('halo_spin', 'f4'), 
    ('halo_id_breadth_first', 'i8'), 
    ('halo_id_depth_first', 'i8'), 
    ('halo_id_tree_root', 'i8'), 
    ('halo_id_orig', 'i8'), 
    ('halo_snap_num', 'i4'), 
    ('halo_id_next_coprog_depthfirst', 'i8'), 
    ('halo_id_last_prog_depthfirst', 'i8'), 
    ('halo_id_last_mainleaf_depthfirst', 'i8'),
    ('halo_rs_klypin', 'f4'),
    ('halo_mvir_all', 'f4'), 
    ('halo_m200b', 'f4'), 
    ('halo_m200c', 'f4'), 
    ('halo_m500c', 'f4'), 
    ('halo_m2500c', 'f4'), 
    ('halo_xoff', 'f4'), 
    ('halo_voff', 'f4'), 
    ('halo_spin_bullock', 'f4'), 
    ('halo_b_to_a', 'f4'), 
    ('halo_c_to_a', 'f4'), 
    ('halo_axisA_x', 'f4'), 
    ('halo_axisA_y', 'f4'), 
    ('halo_axisA_z', 'f4'), 
    ('halo_b_to_a_500c', 'f4'), 
    ('halo_c_to_a_500c', 'f4'), 
    ('halo_axisA_x_500c', 'f4'), 
    ('halo_axisA_y_500c', 'f4'), 
    ('halo_axisA_z_500c', 'f4'), 
    ('halo_t_by_u', 'f4'), 
    ('halo_mass_pe_behroozi', 'f4'), 
    ('halo_mass_pe_diemer', 'f4'), 
    ('halo_macc', 'f4'), 
    ('halo_mpeak', 'f4'), 
    ('halo_vacc', 'f4'), 
    ('halo_vpeak', 'f4'), 
    ('halo_halfmass_scale', 'f4'), 
    ('halo_dmvir_dt_inst', 'f4'), 
    ('halo_dmvir_dt_100myr', 'f4'), 
    ('halo_dmvir_dt_tdyn', 'f4'), 
    ('halo_dmvir_dt_2dtyn', 'f4'), 
    ('halo_dmvir_dt_mpeak', 'f4'), 
    ('halo_scale_factor_mpeak', 'f4'), 
    ('halo_scale_factor_lastacc', 'f4'), 
    ('halo_scale_factor_firstacc', 'f4'), 
    ('halo_mvir_firstacc', 'f4'), 
    ('halo_vmax_firstacc', 'f4'), 
    ('halo_vmax_mpeak', 'f4')
    ])
dtype_slac_bolplanck_rockstar_july19_2015 = dtype_slac_bolshoi_rockstar_july19_2015
dtype_slac_multidark_rockstar_july19_2015 = dtype_slac_bolshoi_rockstar_july19_2015


################################################################################################
################################################################################################





