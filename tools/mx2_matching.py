import numpy as np


def Mx2_DS_Match(reco_track, minerva_data):


    # Variables to return 
    ds_dp = 0 
    bm_index = None
    exit_minerva = False
    deltaX = None
    deltaY = None


    tpc_track_s = reco_track[0]
    tpc_track_e = reco_track[1]

    # Data of tpc track 
    tpc_dx = tpc_track_e[0] - tpc_track_s[0]  
    tpc_dy = tpc_track_e[1] - tpc_track_s[1]
    tpc_dz = tpc_track_e[2] - tpc_track_s[2]
    tpc_trk_length = np.sqrt(np.square(tpc_dx) + np.square(tpc_dy) + np.square(tpc_dz))
    tpc_dir_vec = np.array([tpc_dx,tpc_dy,tpc_dz])/tpc_trk_length
    print(f"2x2 track length {tpc_trk_length:.2f} cm")

    if(tpc_dir_vec[2]<0):
        print("Flipping track direction")
        print(f"Original start {tpc_track_s}")
        print(f"Original end {tpc_track_e}")
        tpc_dir_vec[2] = -1*tpc_dir_vec[2]
        temp = tpc_track_s
        tpc_track_s = tpc_track_e
        tpc_track_e = temp 
        print(f"New start {tpc_track_s}")
        print(f"New end {tpc_track_e}")

    # Hard coded 5 cm (?)
    if(tpc_trk_length>5):
        if(tpc_track_s[2] > 58 or tpc_track_e[2] > 58):
            temp_punchthrough = -1
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
                   ((tpc_track_s[2] > 58) or (tpc_track_e[2] > 58))):
                    mn_dx = mn_track_ex - mn_track_sx
                    mn_dy = mn_track_ey - mn_track_sy
                    mn_dz = mn_track_ez - mn_track_sz
                    mn_length = np.sqrt(np.square(mn_dx) + 
                                        np.square(mn_dy) + 
                                        np.square(mn_dz))    

                    if(mn_length<10):
                        continue 
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

                    if(ds_dp < dot_product and np.abs(extrapdy)< 10):
                        ds_dp = dot_product
                        bm_index = itrack
                        deltaX = extrapdx 
                        deltaY = extrapdy

    return bm_index,ds_dp,deltaX,deltaY,exit_minerva