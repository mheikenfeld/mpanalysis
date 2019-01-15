import numpy as np
import pandas as pd
import iris
import os

def get_timing(track,process=None,measure=None,value=None):
    columns = [('cell', int),
               ('time_start', np.timedelta64 ),
               ('time_end', np.timedelta64),
               (f'{process}_{measure}', np.timedelta64)]
    # create the dataframe from a dict
    df_out=pd.DataFrame({k: pd.Series(dtype=t) for k, t in columns})
    #df_out=pd.DataFrame()
    for cell,track_cell in track.groupby('cell'):
        index=None
        if measure=='max':
            index=track_cell[process].idxmax()
            # reject a maximum smaller than the threshold value
            if track_cell.loc[index,process]<value:
                index=None
            
        elif measure=='min':
            index=track_cell[process].idxmin()
        elif measure=='first':
            if any(track_cell.index[track_cell[process] > value]):
                index=track_cell.index[track_cell[process] > value][0]        
        elif measure=='last':
            if any(track_cell.index[track_cell[process] > value]):
                index=track_cell.index[track_cell[process] > value][-1]        
        lifetime=track_cell['time_cell'].max()
        time_start=track_cell['time'].head(1).values[0]
        time_end=track_cell['time'].tail(1).values[0]
        if index:
            time_process=track_cell.loc[index,'time_cell']
        else:
            time_process=None
            
        df_out=df_out.append({'cell':cell,
                                'time_start': time_start, 
                                'time_end': time_end, 
                                'lifetime': lifetime, 
                                f'{process}_{measure}': time_process},ignore_index=True)
    df_out.reset_index(drop=True,inplace=True)
#     print(df_out)
    return df_out

def latent_heating(track):
    L_sl=334e3
    L_lv=2.26476e6
    L_sv=L_sl+L_lv
    latent_heating=L_lv*abs(track['Condensation']) \
                   -L_lv*abs(track['Evaporation']) \
                   +L_sl*abs(track['Freezing']) \
                   -L_sl*abs(track['Melting']) \
                   +L_sv*abs(track['Deposition']) \
                   -L_sv*abs(track['Sublimation'])

    return latent_heating

def make_timings(track,threshold_value=None):
    track['total_latent_heating']=latent_heating(track)
    timing_latent_heating_max=get_timing(track,process='total_latent_heating',measure='max',value=threshold_value)
    timing_Freezing_first=get_timing(track,process='Freezing',measure='first',value=threshold_value)
    timing_Condensation_max=get_timing(track,process='Condensation',measure='max',value=threshold_value)
    timing_Freezing_max=get_timing(track,process='Freezing',measure='max',value=threshold_value)
    timing_out=timing_Freezing_first[['time_start','time_end','lifetime','cell']]
    
    timing_out['Freezing_first']=timing_Freezing_first['Freezing_first']
    timing_out['Freezing_time_relative']=timing_Freezing_first['Freezing_first'].dt.total_seconds()/timing_Freezing_first['lifetime'].dt.total_seconds()
    timing_out['Condensation_max_relative']=timing_Condensation_max['Condensation_max'].dt.total_seconds()/timing_Condensation_max['lifetime'].dt.total_seconds()
    timing_out['Freezing_max_relative']=timing_Freezing_max['Freezing_max'].dt.total_seconds()/timing_Freezing_max['lifetime'].dt.total_seconds()
    timing_out['Freezing_max']=timing_Freezing_max['Freezing_max']
    timing_out['Condensation_max']=timing_Condensation_max['Condensation_max']        
    timing_out['Heating_max']=timing_latent_heating_max['total_latent_heating_max']        
    timing_out['Heating_max_relative']=timing_latent_heating_max['total_latent_heating_max'].dt.total_seconds()/timing_latent_heating_max['lifetime'].dt.total_seconds()

    timing_out['category']=None
    timing_out.at[timing_out['Freezing_time_relative']==0,'category']='nowarm'
    timing_out.at[timing_out['Freezing_first'].isnull(),'category']='warm'
    timing_out.at[(timing_out['Freezing_time_relative']>0) & (~timing_out['Freezing_first'].isnull()),'category']='warmcold'

    return timing_out

