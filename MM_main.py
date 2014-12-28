from MM_AR.GPS_point_cand_point.GPS_point import ObsPoint 
from MM_AR.HMM_procedure.HMM_KN_Imp_Transition_Probabilities \
    import TransitionWeight
import pickle
import numpy as np
from MM_AR.HMM_procedure.ReadFile import ReadDataPoint
import time
import sys


def Viterbi(datafile, lon_col_id, lat_col_id, timestamp_col_id, 
            gps_mean, gps_std_dev, 
            road_net_shp, road_net_multigraph_pickled, beta): 
    '''
    Viterbi decoder of optimum sequence of the candidate road segments.
    Input file should be a csv file with at least three columns for 
    Latitude, Longitude and Time stamp.  

    Parameters
    ==========
    lon_col_id, lat_col_id, timestamp_col_id  are the column numbers counting 
    from leftmost column as zero
    
    Future work
    =================================
    1. If the path finder process breaks for two consecutive points then it 
    does not remove them from the sequence as Newson and Krumm did manually.
    
    2. Emission probability does not incorporate the general definition
    for emission probability by Oran and Jaillet (2013)
     
    3. Performance not optimized
     
    '''   
     
    MultiGraph = pickle.load(open(road_net_multigraph_pickled, 'rb'))
    
    PointGenerator = ReadDataPoint(datafile, lon_col_id, 
                                   lat_col_id, timestamp_col_id)
    record = 0
    points = []    
    try:
        for GPSrecord in PointGenerator:
    
            if record == 0:   
                print "IF starts"
                print 'rec ', record     
                p1 = ObsPoint( gps_lon = GPSrecord[0], 
                               gps_lat = GPSrecord[1],
                               timestamp = GPSrecord[2], 
                               gps_mean = gps_mean, 
                               gps_std_dev = gps_std_dev, 
                               road_net_shp = road_net_shp )
                
                if len(p1.candidate_points) == 0: #catches first valid point
                    record = 0 
                    continue
                
                points.append(p1)
                
                lnHeadProbVec = np.empty((1,1)) #initiate as list. from 2nd GPS 
                                                #point, becomes NumPy array
                WaveHead = []            
                for cand_pt_0 in points[record].candidate_points:
                    np.append(lnHeadProbVec,
                              [[np.log(cand_pt_0.cand_pt_emission_prob)]], 
                              axis = 0)
                    
                    WaveHead.append([(cand_pt_0.cand_pt_easting, 
                                     cand_pt_0.cand_pt_northing)])
                    
            
            elif record > 0:
                
                if record > 1:
                    del(points[0])
                
                p2 = ObsPoint( gps_lon = GPSrecord[0], 
                               gps_lat = GPSrecord[1],
                               timestamp = GPSrecord[2], 
                               gps_mean = gps_mean, 
                               gps_std_dev = gps_std_dev, 
                               road_net_shp = road_net_shp )  
                points.append(p2)
                row_len = len(points[0].candidate_points)
                col_len = len(points[1].candidate_points)
                
            
                #check and remove the points with no candidate points or 
                #if it is the same as the immediately previous point
                if points[0] == points[1] or col_len == 0:                   
                    record += 1
                    points[1] = points[0] 
                    continue
                
                
                TransitionObjMatrix = np.empty_like([[0]*col_len]*row_len, 
                                                    TransitionWeight)
                TransWeightMatrix = np.empty([row_len,col_len], dtype = float)
                EmissionProbMatrix = np.empty(col_len, dtype=float)
                #array(cand_pt_t_minus_1  X cand_pt_t)
                i = 0
                for cand_pt_t in points[1].candidate_points:
                     
                    j = 0
                    EmissionProbMatrix[i] = np.log(
                        cand_pt_t.cand_pt_emission_prob)
                    for cand_pt_t_minus_1 in points[0].candidate_points:
                        TransitionObjMatrix[j][i] = TransitionWeight(
                                                    cand_pt_t_minus_1, 
                                                    cand_pt_t, MultiGraph, 
                                                    road_net_shp, beta)                                               
                        TransWeightMatrix[j][i] =   \
                                        TransitionObjMatrix[j][i].transition_wt 
                        #(jXi)matrix
                        
                        j += 1
                    i += 1
                try:
                    
                    sum_row_wt = np.sum(TransWeightMatrix, axis=1)
                    for i in range(row_len):
                        if sum_row_wt[i-1] == 0:
                            sum_row_wt[i-1] = 1.0
                    ln_sum_row_wt = np.log(sum_row_wt)
                                                    #(1xj)),   log(sum(row))
                
                    print 'record# before crashpoint',record
                    p = np.log(TransWeightMatrix).transpose()
                    lnTransProbMatrix = (p - ln_sum_row_wt).transpose()
                    
                
                
                lnHTE = lnHeadProbVec + lnTransProbMatrix + EmissionProbMatrix
                lnHeadProbVec = np.max(lnHTE, axis=0).reshape((
                                                            lnHTE.shape[1],1)) 
                
                #find the row-column positions of the maximum values for 
                #max-probable cand_pt_t
                
                WaveHead_temp = []
                lnHTE_trnsps = np.transpose(lnHTE)
                for col_argmax in xrange(col_len):
                    row_argmax = np.argmax(lnHTE_trnsps[col_argmax])
                    bridge = TransitionObjMatrix[row_argmax][col_argmax].\
                                                            shortest_path_nodes
                
                    if record == 1:
                        WaveHead_temp.append(bridge)
                    elif lnHeadProbVec[col_argmax] == 0:
                        WaveHead_temp.append('out of network')
                    else:
                        WaveHead_temp.append(WaveHead[row_argmax]+bridge[:])
                
                WaveHead = WaveHead_temp
                
                
            record += 1
            
        max_prob = np.max(lnHeadProbVec)
        max_prob_path = WaveHead[np.argmax(lnHeadProbVec)]
        return max_prob, max_prob_path
    
    except UnboundLocalError:
        print 'Stationary object'
        return None
        


if __name__ == '__main__':
    pass
