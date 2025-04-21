import numpy as np


def Mx2_DS_Match(reco_track, minerva_data, reco_vtx,mc=True,verbose=False):


    # Hard-coded values:
    trk_vtx_dist = 5 
    tpc_trk_min_l = 5
    minerva_trk_min_l = 10
    tpc_fv_dist = 60
    min_disp = 15


    if(mc):
        mnvoffsetx = 0
        mnvoffsety = 0 
    else:
        mnvoffsetx = -10
        mnvoffsety = 5

    # Variables to return 
    ds_dp = 0 
    bm_index = None
    exit_minerva = False
    deltaX = None
    deltaY = None


    tpc_track_s = reco_track[0]
    tpc_track_e = reco_track[1]


    ixn_vx = reco_vtx[0]
    ixn_vy = reco_vtx[1]
    ixn_vz = reco_vtx[2]

    dvtx_x = np.abs(tpc_track_s[0] - ixn_vx)
    dvtx_y = np.abs(tpc_track_s[1] - ixn_vy)
    dvtx_z = np.abs(tpc_track_s[2] - ixn_vz)
    dvtx = np.sqrt(dvtx_x**2 + dvtx_y**2 + dvtx_z**2)

    if(dvtx>trk_vtx_dist):
        return bm_index,ds_dp,deltaX,deltaY,exit_minerva
    else: 

        # Data of tpc track 
        tpc_dx = tpc_track_e[0] - tpc_track_s[0]  
        tpc_dy = tpc_track_e[1] - tpc_track_s[1]
        tpc_dz = tpc_track_e[2] - tpc_track_s[2]
        tpc_trk_length = np.sqrt(np.square(tpc_dx) + np.square(tpc_dy) + np.square(tpc_dz))
        tpc_dir_vec = np.array([tpc_dx,tpc_dy,tpc_dz])/tpc_trk_length
        if(verbose): print(f"2x2 track length {tpc_trk_length:.2f} cm")

        if(tpc_dir_vec[2]<0):

            if(verbose):
                print("Flipping track direction")
                print(f"Original start {tpc_track_s}")
                print(f"Original end {tpc_track_e}")
            tpc_dir_vec[2] = -1*tpc_dir_vec[2]
            temp = tpc_track_s
            tpc_track_s = tpc_track_e
            tpc_track_e = temp 
            if(verbose):
                print(f"New start {tpc_track_s}")
                print(f"New end {tpc_track_e}")


        if(tpc_trk_length> tpc_trk_min_l):
            if(np.abs(tpc_track_s[2]) > tpc_fv_dist or np.abs(tpc_track_e[2]) > tpc_fv_dist):
                nminerva_tracks = np.sum(minerva_data['rec.nd.minerva.ixn.tracks..length'])
                for itrack in range(nminerva_tracks):
                    mn_track_sx = minerva_data['rec.nd.minerva.ixn.tracks.start.x'][itrack]
                    mn_track_sy = minerva_data['rec.nd.minerva.ixn.tracks.start.y'][itrack]
                    mn_track_sz = minerva_data['rec.nd.minerva.ixn.tracks.start.z'][itrack]

                    mn_track_ex = minerva_data['rec.nd.minerva.ixn.tracks.end.x'][itrack]
                    mn_track_ey = minerva_data['rec.nd.minerva.ixn.tracks.end.y'][itrack]
                    mn_track_ez = minerva_data['rec.nd.minerva.ixn.tracks.end.z'][itrack]

                    if(mn_track_ez>300):
                        exit_minerva=True


                    # Only DS track
                    if(mn_track_sz > 0 and mn_track_ez > 0 and
                    ((tpc_track_s[2] > tpc_fv_dist ) or (tpc_track_e[2] > tpc_fv_dist ))):
                        mn_dx = mn_track_ex - mn_track_sx
                        mn_dy = mn_track_ey - mn_track_sy
                        mn_dz = mn_track_ez - mn_track_sz
                        mn_length = np.sqrt(np.square(mn_dx) + 
                                            np.square(mn_dy) + 
                                            np.square(mn_dz))    

                        if(mn_length<minerva_trk_min_l):
                            continue 

                        else: 
                            mn_dir = np.array([
                                mn_dx,
                                mn_dy,
                                mn_dz
                            ])/mn_length
                            dot_product = np.dot(tpc_dir_vec,mn_dir)  
                            extrapdz = mn_track_sz - tpc_track_e[2]
                            extrapdy = tpc_dy/tpc_dz*(extrapdz) + tpc_track_e[1] - mn_track_sy
                            extrapdx = tpc_dx/tpc_dz*(extrapdz) + tpc_track_e[0] - mn_track_sx
                            diffextrap = np.sqrt(np.square(extrapdy - mn_track_sy))

                            ang_xz_diff = np.abs(np.arctan(mn_dx/mn_dz) - np.arctan(tpc_dx/tpc_dz))
                            ang_yz_diff = np.abs(np.arctan(mn_dy/mn_dz) - np.arctan(tpc_dy/tpc_dz))

                            if(ds_dp < dot_product and 
                               np.abs(extrapdy - mnvoffsety) < min_disp and
                               np.abs(extrapdx - mnvoffsetx) < min_disp and 
                               ang_xz_diff < 0.06 and
                               ang_yz_diff < 0.06):
                                    ds_dp = dot_product
                                    bm_index = itrack
                                    deltaX = extrapdx 
                                    deltaY = extrapdy

        return bm_index,ds_dp,deltaX,deltaY,exit_minerva