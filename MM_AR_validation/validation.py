'''
Created on Dec 7, 2014

@author: asr13006
'''
import pandas as pd
import MM_main
from dateutil.parser import parse
import pytz
import datetime
import fiona
from MM_AR.GPS_point_cand_point.GPS_point import lat_long_to_UTM_point
import  matplotlib.pyplot as plt

class validation(object):
    '''Validation class for MM algo'''
    def __init__(self, raw_shuttle_gps_file):
        self.csv_file = raw_shuttle_gps_file
        
    def timezone_converter(self, input_dt, target_tz='US/Eastern', 
                                            current_tz='UTC'):
        '''
        input_dt is an datetime.datetime object
        assumes timezone by default. See parameter values
        '''
        current_tz = pytz.timezone(current_tz)
        target_tz = pytz.timezone(target_tz)
        target_dt = current_tz.localize(input_dt).astimezone(target_tz)
        return target_tz.normalize(target_dt) 

    def create_csv_freq_x(self):
        data_pd = pd.read_csv(self.csv_file)
        return data_pd
    
    def create_parsed_local_time_col(self, data_pd):
        '''
        call the three column file as lean csv file
        '''
        selected_columns_df = data_pd[['latitude','longitude',  
                                       'utc_time_stamp']]
        selected_columns_df = selected_columns_df.sort('utc_time_stamp')
        selected_columns_df.reset_index(drop=True, inplace=True)
        selected_columns_df['local_time_parsed'] = ''
        
        for idx in range(len(selected_columns_df)):
            utc_parsed = parse(selected_columns_df[
                                'utc_time_stamp'].iloc[idx])             
            selected_columns_df.loc[idx, 'local_time_parsed'] = \
                self.timezone_converter(utc_parsed)
        return selected_columns_df 
    
    def create_time_delta_col(self, selected_columns_df):
        selected_columns_df['time_delta'] = ''
        selected_columns_df['time_delta'] = \
            selected_columns_df['local_time_parsed'].diff()
        return selected_columns_df
    
    def get_unique_dates(self, datetime_col):
        '''
        return the unique dates in the dataframe
        '''
        return datetime_col.map(pd.Timestamp.date).unique()
    
    def chunk_data(self, dataframe, start_dt_local, duration):
        '''
        returns data chunk for the given day
        duration in minutes
        
        Parameters
        ==========
        start_dt is datetime.datetime object
        '''
        end_dt_local = start_dt_local + duration
        start_dt_utc = self.timezone_converter(start_dt_local, 
                                target_tz='UTC', 
                                current_tz='US/Eastern') 
        end_dt_utc = self.timezone_converter(end_dt_local, 
                                target_tz='UTC', 
                                current_tz='US/Eastern') 
        
        start_dt_str = start_dt_utc.replace(tzinfo=None).isoformat()
        end_dt_str = end_dt_utc.replace(tzinfo=None).isoformat()
        
        dataframe = dataframe[dataframe['utc_time_stamp'] > start_dt_str]
        chunk = dataframe[dataframe['utc_time_stamp'] < end_dt_str]
        return chunk
        
    def resampler_for_desired_freq(self, chunk, frequency_second):
        '''
        Resamples from the time-chunked dataset to obtain desired frequency
        
        Parameters
        ==========
        frequency in seconds
        start_time = to be input as 'YYYY:mm:ddThh:mm:ss' 
        '''
        chunk_durn_sec = parse(chunk['utc_time_stamp'].max()) -  \
                    parse(chunk['utc_time_stamp'].min())
        start_index = chunk.index[0]
        median_freq = chunk['time_delta'].median()
        jump = frequency_second / (median_freq.astype('timedelta64[s]'))

        up_to_row = start_index + (chunk_durn_sec.total_seconds() /  
                                                        jump)

        select_rows = range(start_index, chunk.index.max(), int(jump))
        resampled = chunk.ix[select_rows]
        return resampled
        
    def mapmatch_resampled_chunk(self, resampled_chunk, gps_mean, gps_std_dev, 
                                road_net_shp,road_net_multigraph_pickled,
                                beta):
        csv_file = 'Resampled_chunk_{}'.format(self.csv_file)
        resampled_chunk.to_csv(csv_file)
        try:
            map_matched = MM_main.Viterbi(csv_file, 
                            lon_col_id=2, 
                            lat_col_id=1, 
                            timestamp_col_id=3, 
                            gps_mean=gps_mean, 
                            gps_std_dev= gps_std_dev, 
                            road_net_shp=road_net_shp, 
                            road_net_multigraph_pickled=  \
                                                    road_net_multigraph_pickled,
                            beta=beta)
            return map_matched[1]
        except TypeError:
            return None 
        
    def plot_roadnetwork(self, road_shp_file, ax, fig):
        '''
        plots the road network with matplotlib
        fig = matplotlib.pyplot.figure()
        ax = matplotlib.pyplot.axes()
        '''
        with fiona.open(road_shp_file) as road_shp:
            for item in road_shp:                
                line = item['geometry']['coordinates']
                (line_easting, line_northing) = zip(*line)
                ax.plot(line_easting, line_northing, 'k')
                        
    def plot_map_matched_path_points(self, map_matched_seq, ax, fig, label=None):
        '''
        fig = matplotlib.pyplot.figure()
        ax = matplotlib.pyplot.axes()
        '''
        (easting, northing) = zip(*map_matched_seq)
        ax.plot(easting, northing, 'yd', label=label)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=5, fontsize='x-small')
        
    def plot_raw_gps_seq(self, chunk, ax, fig, legend='ro', label=None):
        '''
        fig = matplotlib.pyplot.figure()
        ax = matplotlib.pyplot.axes()
        '''
        utm_pts = map(lat_long_to_UTM_point, 
                      chunk['longitude'],
                      chunk['latitude'])
        (easting, northing) = zip(*utm_pts)
        ax.plot(easting, northing, legend, label=label)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=5, fontsize='x-small')
        return min(easting), max(easting), min(northing), max(northing)
    
    def validation_main(self, start_time_local, end_time_local,  
                        road_net_shp, road_net_multigraph_pickled, 
                        gps_std_dev=30, beta=1, duration_minute=60, 
                        freq_sec=60):
        '''start_time_local and end_time_local are datetime.time object'''
        dur_td = datetime.timedelta(minutes=duration_minute)
        
        df = self.create_csv_freq_x()
        df = self.create_parsed_local_time_col(df)
        df = self.create_time_delta_col(df) 
        dates = self.get_unique_dates(df['local_time_parsed'])
        for tgt_date in dates:
            start_datetime = datetime.datetime(tgt_date.year, 
                                               tgt_date.month,
                                               tgt_date.day,
                                               start_time_local.hour,
                                               start_time_local.minute,
                                               start_time_local.second)
            end_datetime = datetime.datetime(tgt_date.year, 
                                               tgt_date.month,
                                               tgt_date.day,
                                               end_time_local.hour,
                                               end_time_local.minute,
                                               end_time_local.second)
            chunk_start_dt = start_datetime
           
            while chunk_start_dt <= end_datetime :
                chunk_name = '{file}_{year}-{month}-{day}-'  \
                            '{hour}-{minute}-{second}_'  \
                            'dur{dur_min}min_freq{freq_sec}sec'.format(
                                            file=self.csv_file[:-4],
                                            year=chunk_start_dt.year,
                                            month=chunk_start_dt.month,
                                            day=chunk_start_dt.day,
                                            hour=chunk_start_dt.hour,
                                            minute=chunk_start_dt.minute,
                                            second=chunk_start_dt.second,
                                            dur_min=str(duration_minute),
                                            freq_sec=str(freq_sec))
                
                chunk = self.chunk_data(df, chunk_start_dt, dur_td)
                resampled_chunk = self.resampler_for_desired_freq(chunk, 
                                                                  freq_sec)

                map_matched_seq = self.mapmatch_resampled_chunk(
                                            resampled_chunk, 
                                            gps_mean=0, 
                                            gps_std_dev=gps_std_dev, 
                                            road_net_shp=road_net_shp, 
                                            road_net_multigraph_pickled=  \
                                                road_net_multigraph_pickled,
                                            beta=beta)
                fig = plt.figure()
                ax = plt.axes()
                self.plot_roadnetwork(road_net_shp, ax, fig)
                (east_min, east_max, north_min, north_max) =  \
                        self.plot_raw_gps_seq(chunk, ax, fig, legend='cx', 
                                              label='chunk')
               
                self.plot_map_matched_path_points(map_matched_seq, ax, fig,
                                                  label='map-matched')
                self.plot_raw_gps_seq(resampled_chunk, ax, fig, legend='rH',
                                      label='re-sampled chunk')
                plt.xlim(east_min-200, 200+east_max)
                plt.ylim(north_min-200, 200+north_max)
                fig.suptitle(chunk_name, fontsize=10, fontweight='bold')
                fig.savefig('plot_{}'.format(chunk_name))
                plt.close()
                
                #while loop control condition
                chunk_start_dt = min(chunk_start_dt + dur_td, end_datetime)
                        
        return None
    
if __name__ == '__main__':
    import sys
    lab_rat = '032629uconn201.csv'
    val = validation(lab_rat)
    start_time_local = datetime.time(15,30,0)
    end_time_local = datetime.time(16,30,0,)
    road_net_shp = '../MM_AR/Relevant_files/LineString_Road_Network_UTM.shp'
    road_net_multigraph_pickled = "../MM_AR/Preprocess/MultiGraph.p"
    
    for lab_rat_id in [1, 2]:
        sys.stdout= open('output_032629uconn20{}.txt'.format(lab_rat_id), 'w')
        for freq_sec in [30, 60, 90, 180, 210, 240] :   
            val.validation_main(start_time_local, end_time_local,  
                                road_net_shp, road_net_multigraph_pickled,
                                gps_std_dev=30, beta=1, duration_minute=30, 
                                freq_sec=freq_sec)
        sys.stdout.close()