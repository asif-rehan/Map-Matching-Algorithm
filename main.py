from MM_AR.GPS_point_cand_point.GPS_point import ObsPoint 
from MM_AR.HMM_procedure.HMM_KN_Imp_Transition_Probabilities import TransitionWeight
import pickle, numpy as np
from MM_AR.HMM_procedure.ReadFile import ReadDataPoint
import time
import sys

sys.stdout= open('output_{0}_{1}.txt'.format(time.strftime('%Y-%m-%dT%H.%M.%S', 
                                            time.localtime()), 
                                             'red_20141014'), 'w')

start_time = time.time()
def main(datafile, gps_mean, gps_std_dev, road_net_shp, beta):    
    print "MAIN RUNS"
    MultiGraph_pickled = "MM_AR/Preprocess/MultiGraph.p" 
    MultiGraph = pickle.load(open(MultiGraph_pickled, 'rb'))
    PointGenerator = ReadDataPoint(datafile)

    record = 0
    points = []    
    
    for GPSrecord in PointGenerator:
        print 'rec ', record
        if record == 0:   
            print "IF starts"
            print 'rec ', record     
            p1 = ObsPoint( gps_lon = GPSrecord[0], 
                           gps_lat = GPSrecord[1],
                           timestamp = GPSrecord[2], 
                           gps_mean = gps_mean, 
                           gps_std_dev = gps_std_dev, 
                           road_net_shp = road_net_shp )
            
            points.append(p1)
            
            #print 'print len(points)',len(points)
            #print 'len p1.candidate points', len(p1.candidate_points)
            
            lnHeadProbVec = np.empty((1,1)) #initiate as np.array
                                            
            WaveHead = []            
            for cand_pt_0 in points[record].candidate_points:
                np.append(lnHeadProbVec,
                          [[np.log(cand_pt_0.cand_pt_emission_prob)]], 
                          axis = 0)
                WaveHead.append(None)
                print (cand_pt_0.cand_pt_easting, cand_pt_0.cand_pt_northing)
            print "lnHeadProbVec", lnHeadProbVec
            print "WaveHead", WaveHead
            #print 'IF END'
        
        elif record>0:
            print "ELIF START"            
            
            if record > 1:
                del(points[0])
            
            p2 = ObsPoint( gps_lon = GPSrecord[0], 
                           gps_lat = GPSrecord[1],
                           timestamp = GPSrecord[2], 
                           gps_mean = gps_mean, 
                           gps_std_dev = gps_std_dev, 
                           road_net_shp = road_net_shp )  
            points.append(p2)
            for p in points:
                print 'GPS point', p.gps_easting, p.gps_northing
            row_len = len(points[0].candidate_points)
            col_len = len(points[1].candidate_points)
            print "(row_len, col_len) = (", row_len, col_len, ")"
            
            for pos in xrange(row_len):
                print ('points[0]',
                points[0].candidate_points[pos].cand_pt_easting,
                        points[0].candidate_points[pos].cand_pt_northing) 
            for pos in xrange(col_len):
                print ('points[1]',
                points[1].candidate_points[pos].cand_pt_easting,
                        points[1].candidate_points[pos].cand_pt_northing)
            
            if points[0] == points[1] or col_len == 0:
                record += 1
                continue
            
            #print (p2.__dict__)          
            #print 'len p2.candidate points', len(p2.candidate_points)
            #print 'print len(points)', "*"*len(points)
           
            TransitionObjMatrix = np.empty_like([[0]*col_len]*row_len, TransitionWeight)
            TransWeightMatrix = np.empty([row_len,col_len], dtype = float)
            EmissionProbMatrix = np.empty(col_len, dtype=float)
            #array(cand_pt_t_minus_1  X cand_pt_t)
            print "len(points[1].candidate_points) = ",len(points[1].candidate_points)
            print "len(points[0].candidate_points) = ",len(points[0].candidate_points)
            i = 0
            for cand_pt_t in points[1].candidate_points:
                 
                j = 0
                EmissionProbMatrix[i] = np.log(cand_pt_t.cand_pt_emission_prob)
                for cand_pt_t_minus_1 in points[0].candidate_points:
                    TransitionObjMatrix[j][i] = TransitionWeight(cand_pt_t_minus_1, 
                                    cand_pt_t, MultiGraph, road_net_shp, beta)                                               
                    TransWeightMatrix[j][i] = TransitionObjMatrix[j][i].transition_wt 
                    #(jXi)matrix
                    #print "j = ", j
                    j += 1
                #print "i = ", i
                i += 1
            print "TransWeightMatrix",TransWeightMatrix
            print record, "EmissionProbMatrix",EmissionProbMatrix
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
                
            except RuntimeWarning:
                print RuntimeWarning
                print 'p=', p, 'record=', record
            
            print "log_sum_row_wt", ln_sum_row_wt
            print "TransitionProbMatrix", lnTransProbMatrix
            
            
            lnHTE = lnHeadProbVec + lnTransProbMatrix + EmissionProbMatrix
            print 'lnHTE', lnHTE
            lnHeadProbVec = np.max(lnHTE, axis = 0).reshape((lnHTE.shape[1],1)) 
            #head prob values
            #find the row-column positions of the maximum values for 
            #max-probable cand_pt_t
            print "lnHeadProbVec", lnHeadProbVec
            
            WaveHead_temp = []
            lnHTE_trnsps = np.transpose(lnHTE)
            for col_argmax in xrange(col_len):
                row_argmax = np.argmax(lnHTE_trnsps[col_argmax])
                bridge = TransitionObjMatrix[row_argmax][col_argmax].\
                                                        shortest_path_nodes
                print 'bridge', bridge
                if record == 1:
                    WaveHead_temp.append(bridge)
                elif lnHeadProbVec[col_argmax] == 0:
                    WaveHead_temp.append('out of network')
                else:
                    WaveHead_temp.append(WaveHead[row_argmax]+bridge[:])
            
            WaveHead = WaveHead_temp
            print "WaveHead", WaveHead
            
        record += 1
        
    max_prob = np.max(lnHeadProbVec)
    max_prob_path = WaveHead[np.argmax(lnHeadProbVec)]
    print " max_prob  ", max_prob
    print " max_prob_path  ", max_prob_path
    return max_prob, max_prob_path
    #assumption: no tied value for likelihood calculation

main(datafile= "MM_AR/Relevant_files/phnGPS_orng.csv", 
     gps_mean = 0, gps_std_dev=50, 
     road_net_shp = "MM_AR/Relevant_files/LineString_Road_Network_UTM.shp",
beta=1)

print "--- {0} seconds ---".format(time.time() - start_time)
sys.stdout.close()
sys.stderr = sys.__stderr__
    