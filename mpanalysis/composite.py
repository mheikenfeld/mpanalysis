import numpy as np
import pandas as pd
import iris
import os

def get_timing(track,process=None,measure=None,value=None):
    columns = [('cell', int),
               ('time_start', np.timedelta64(1,'ns') ),
               ('time_end', np.timedelta64(1,'ns')),
               (f'{process}_{measure}', np.timedelta64(1,'ns'))]
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

def make_timings(track,threshold_value_processes=None,threshold_value_mass=None):
    track['total_latent_heating']=latent_heating(track)
    timing_latent_heating_max=get_timing(track,process='total_latent_heating',measure='max',value=threshold_value_processes)
    timing_Freezing_first=get_timing(track,process='Freezing',measure='first',value=threshold_value_processes)
    timing_Condensation_max=get_timing(track,process='Condensation',measure='max',value=threshold_value_processes)
    timing_Freezing_max=get_timing(track,process='Freezing',measure='max',value=threshold_value_processes)
    timing_Rain_max=get_timing(track,process='Rain formation',measure='max',value=threshold_value_processes)
    
    timing_totalwater_max=get_timing(track,process='total_water_mass',measure='max',value=threshold_value_mass)
    timing_liquidwater_max=get_timing(track,process='liquid_water_mass',measure='max',value=threshold_value_mass)
    timing_frozenwater_max=get_timing(track,process='frozen_water_mass',measure='max',value=threshold_value_mass)

    timing_out=timing_Freezing_first[['time_start','time_end','lifetime','cell']]
    
    timing_out['Freezing_first']=timing_Freezing_first['Freezing_first']
    timing_out['Freezing_time_relative']=timing_Freezing_first['Freezing_first'].dt.total_seconds()/timing_Freezing_first['lifetime'].dt.total_seconds()
    timing_out['Condensation_max_relative']=timing_Condensation_max['Condensation_max'].dt.total_seconds()/timing_Condensation_max['lifetime'].dt.total_seconds()
    timing_out['Freezing_max_relative']=timing_Freezing_max['Freezing_max'].dt.total_seconds()/timing_Freezing_max['lifetime'].dt.total_seconds()
    timing_out['Freezing_max']=timing_Freezing_max['Freezing_max']
    timing_out['Condensation_max']=timing_Condensation_max['Condensation_max']   
    timing_out['Heating_max']=timing_latent_heating_max['total_latent_heating_max']        
    timing_out['Heating_max_relative']=timing_latent_heating_max['total_latent_heating_max'].dt.total_seconds()/timing_latent_heating_max['lifetime'].dt.total_seconds()
    
    timing_out['total_max']=timing_totalwater_max['total_water_mass_max']
    timing_out['liquid_max']=timing_liquidwater_max['liquid_water_mass_max']
    timing_out['frozen_max']=timing_frozenwater_max['frozen_water_mass_max']
    return timing_out

def make_values(track,timings,threshold_value_processes=None,threshold_value_mass=None):
    values=timings[['cell','time_start','time_end']].copy()
    values['Freezing_max']=None
    values['Condensation_max']=None
    values['Rain_max']=None
    values['Heating_max']=None
    values['total_max']=None
    values['liquid_max']=None
    values['frozen_max']=None

    
    for i in values.index:
        cell=values.at[i,'cell']
        values.at[i,'Freezing_max']=track.loc[track['cell']==cell,'Freezing'].max()
        values.at[i,'Condensation_max']=track.loc[track['cell']==cell,'Condensation'].max()
        values.at[i,'Rain_max']=track.loc[track['cell']==cell,'Rain formation'].max()
        values.at[i,'Heating_max']=track.loc[track['cell']==cell,'total_latent_heating'].max()
        values.at[i,'total_max']=track.loc[track['cell']==cell,'total_water_mass'].max()
        values.at[i,'liquid_max']=track.loc[track['cell']==cell,'liquid_water_mass'].max()
        values.at[i,'frozen_max']=track.loc[track['cell']==cell,'frozen_water_mass'].max()
        values.at[i,'frozen_fraction_max']=track.loc[track['cell']==cell,'frozen_water_fraction'].max()
        values.at[i,'frozen_fraction_min']=track.loc[track['cell']==cell,'frozen_water_fraction'].min()

    return values


def categorise(track,timing_in,values_in):
    timing_out=timing_in
    category_out=timing_in[['cell','time_start','time_end']].copy()
#     display(track.head(5))
#     display(timing_in.head(5))
    category_out['category']=None
    Freezing_at_first=timing_in['Freezing_time_relative']==0
    no_Freezing=timing_in['Freezing_first'].isnull()
    no_ice=values_in['frozen_fraction_max']<0.03
    Freezing_later=timing_out['Freezing_time_relative']>0 & ~timing_in['Freezing_first'].isnull()
    no_warm=values_in['frozen_fraction_min']>0.2
#     print(Freezing_at_first)
    for i in timing_out.index:
        if no_warm[i]:
            category_out.at[i,'category']='cold'
        elif no_Freezing[i] and no_ice[i]:
            category_out.at[i,'category']='warm'
        elif Freezing_later[i] and not no_ice[i]:
            category_out.at[i,'category']='warmcold'
        else:
            category_out.at[i,'category']='rest'
    return category_out

def composite_Processes_time(timing,model,case,cells='all',type_mask='_TWC',shiftby='initiation',n_left=100,n_right=100):
    list_processes=['Condensation','Evaporation','Freezing','Melting','Deposition','Sublimation','Rain formation']
    from copy import deepcopy
    if cells=='all':
        cells=timing['cell'].values.astype(int)

    cell_0=timing['cell'].values.astype(int)[0]
    savedir_0=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}')
    Processes_sum_0=iris.load(os.path.join(savedir_0,f'Processes_lumped_profile{type_mask}.nc'))
    Processes_sum_0=Processes_sum_0.extract(list_processes)

    Processes_sum_Composite=iris.cube.CubeList()
    
    n_total=n_left+n_right
    
    Processes_shifted={}
    for process in Processes_sum_0:
        Processes_shifted[process.name()]=[]
        process_composite=iris.cube.Cube(np.zeros((n_total,process.shape[1])),long_name=process.name(),units=process.units)
        time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
        process_composite.add_dim_coord(time_coord,data_dim=0)
        process_composite.add_dim_coord(process.coord('geopotential_height'),data_dim=1)
        Processes_sum_Composite.append(process_composite)
    Processes_mean_Composite=deepcopy(Processes_sum_Composite)
    
    #print(cells)
    if list(cells):
        for i_cell,cell in enumerate(cells):
    #         print('###############:')
    #         print('cell:',cell)
            savedir_i=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}')

            filename_profile=os.path.join(savedir_i,f'Processes_lumped_profile{type_mask}.nc')
            filename_track=os.path.join(savedir_i,f'track_processes_lumped_integrated_TWC{type_mask}.h5')
            if (glob.glob(filename_profile) and glob.glob(filename_track)):
                Processes_sum=iris.load(filename_profile)
                Processes_sum=Processes_sum.extract(list_processes)

                Track_Processes_integrated=pd.read_hdf(filename_track,'table')

                timing_cell=timing[timing['cell']==cell]
                time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
                time_maximum_heating=timing_cell['Heating_max'].dt.total_seconds().values[0]
                lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
                time_axis=Track_Processes_integrated['time_cell'].dt.total_seconds().values

                length=len(time_axis)

                if shiftby is 'initiation':
                    shift=0
                elif shiftby is 'glaciation':
                    if not np.isnan(time_freezing):        
                        shift=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
                    else:
                        shift=length
                elif shiftby is 'maximum_heating':
                    if not np.isnan(time_maximum_heating):        
                        shift=np.argwhere(abs(time_axis-time_maximum_heating)<1)[0][0]
                    else:
                        shift=length

                pad_left=n_left-shift
                pad_right=n_right-(len(time_axis)-shift)
                cut_left=0
                cut_right=0
                if pad_left<0:
                    cut_left=-pad_left
                    pad_left=0

                if pad_right<0:
                    cut_right=-pad_right
                    pad_right=0

                for process in Processes_sum_Composite:
                    process_i=Processes_sum.extract_strict(process.name()).data
                    process_i=process_i[0+cut_left:process_i.shape[0]-cut_right]
                    process_i[process_i>1e20]=0
                    process_shifted=np.pad(process_i,((pad_left,pad_right),(0,0)), 'constant', constant_values=((0, 0),(0,0)))
                    Processes_shifted[process.name()].append(process_shifted)

        for process in Processes_sum_Composite:
            process.data=np.sum(Processes_shifted[process.name()],axis=0)

        for process in Processes_mean_Composite:
            process.data=np.mean(Processes_shifted[process.name()],axis=0)

    return Processes_sum_Composite,Processes_mean_Composite

def composite_Hydrometeors_time(timing,model,case,cells='all',type_mask='_TWC',shiftby='initiation',n_left=100,n_right=100):
    from copy import deepcopy
    if cells=='all':
        cells=timing['cell'].values.astype(int)

    cell_0=timing['cell'].values.astype(int)[0]
    savedir_0=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}')
    Hydrometeors_sum_0=iris.load(os.path.join(savedir_0,f'Hydrometeors_profile{type_mask}.nc'))

    Hydrometeors_sum_Composite=iris.cube.CubeList()
    
    n_total=n_left+n_right
    
    Hydrometeors_shifted={}
    for hydrometeor in Hydrometeors_sum_0:
        Hydrometeors_shifted[hydrometeor.name()]=[]
        hydrometeor_composite=iris.cube.Cube(np.zeros((n_total,hydrometeor.shape[1])),long_name=hydrometeor.name(),units=hydrometeor.units)
        time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
        hydrometeor_composite.add_dim_coord(time_coord,data_dim=0)
        hydrometeor_composite.add_dim_coord(hydrometeor.coord('geopotential_height'),data_dim=1)
        Hydrometeors_sum_Composite.append(hydrometeor_composite)
    Hydrometeors_mean_Composite=deepcopy(Hydrometeors_sum_Composite)
    
    if list(cells):
        for i_cell,cell in enumerate(cells):
            savedir_i=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}')

            filename_profile=os.path.join(savedir_i,f'Hydrometeors_profile{type_mask}.nc')
            filename_track=os.path.join(savedir_i,f'track_hydrometeros_integrated{type_mask}.h5')
            if (glob.glob(filename_profile) and glob.glob(filename_track)):
                Hydrometeors_sum=iris.load(filename_profile)

                Track_Hydrometeors_integrated=pd.read_hdf(filename_track,'table')

                timing_cell=timing[timing['cell']==cell]
                time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
                lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
                time_maximum_heating=timing_cell['Heating_max'].dt.total_seconds().values[0]
                time_axis=Track_Hydrometeors_integrated['time_cell'].dt.total_seconds().values

                length=len(time_axis)

                if shiftby is 'initiation':
                    shift=0
                elif shiftby is 'glaciation':
                    if not np.isnan(time_freezing):        
                        shift=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
                    else:
                        shift=length
                elif shiftby is 'maximum_heating':
                    if not np.isnan(time_maximum_heating):        
                        shift=np.argwhere(abs(time_axis-time_maximum_heating)<1)[0][0]
                    else:
                        shift=length
                        
                pad_left=n_left-shift
                pad_right=n_right-(len(time_axis)-shift)
                cut_left=0
                cut_right=0
                if pad_left<0:
                    cut_left=-pad_left
                    pad_left=0

                if pad_right<0:
                    cut_right=-pad_right
                    pad_right=0

                for hydrometeor in Hydrometeors_sum_Composite:
                    hydrometeor_i=Hydrometeors_sum.extract_strict(hydrometeor.name()).data
                    hydrometeor_i=hydrometeor_i[0+cut_left:hydrometeor_i.shape[0]-cut_right]
                    hydrometeor_i[hydrometeor_i>1e20]=0
                    hydrometeor_shifted=np.pad(hydrometeor_i,((pad_left,pad_right),(0,0)), 'constant', constant_values=((0, 0),(0,0)))
                    Hydrometeors_shifted[hydrometeor.name()].append(hydrometeor_shifted)

        for hydrometeor in Hydrometeors_sum_Composite:
            hydrometeor.data=np.sum(Hydrometeors_shifted[hydrometeor.name()],axis=0)

        for hydrometeor in Hydrometeors_mean_Composite:
            hydrometeor.data=np.mean(Hydrometeors_shifted[hydrometeor.name()],axis=0)

    return Hydrometeors_sum_Composite,Hydrometeors_mean_Composite


def composite_Precipitation_time(timing,model,case,cells='all',type_mask='_TWC',shiftby='initiation',n_left=100,n_right=100):
    from copy import deepcopy
    
    if cells=='all':
        cells=timing['cell'].values.astype(int)

    cell_0=timing['cell'].values.astype(int)[0]
    
    list_precip_types= ['surface_precipitation','surface_precipitation_accumulated']#,'surface_precipitation_instantaneous']
    savedir_0=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}')
    Precipitation_0=pd.read_hdf(os.path.join(savedir_0,f'track_precipitation{type_mask}.h5','table'))
    
    if model=='WRF':
        Precipitation_0['surface_precipitation']=1/3600*Precipitation_0['surface_precipitation_average']

    Precipitation_Composite_sum=pd.DataFrame()
    Precipitation_shifted={}

    n_total=n_left+n_right
    for precip_type in list_precip_types:
        Precipitation_shifted[precip_type]=[]
        Precipitation_Composite_sum[precip_type]=np.zeros((n_total))
        Precipitation_Composite_sum['time']=np.arange(-n_left,n_right)
    Precipitation_Composite_mean=deepcopy(Precipitation_Composite_sum)
    for i_cell,cell in enumerate(cells):
        savedir_i=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}')
#        filename_profile=os.path.join(savedir_i,f'Hydrometeors_profile_TWC_{cell}.nc')
        filename_track=os.path.join(savedir_i,f'track_precipitation{type_mask}.h5')
        if (glob.glob(filename_track)):
            #precip_sum=iris.load(filename_profile)
            precip_sum=pd.read_hdf(filename_track,'table')
            if model=='WRF':
                precip_sum['surface_precipitation']=1/3600*precip_sum['surface_precipitation_average']

            timing_cell=timings[model][case][timings[model][case]['cell']==cell]
            time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values
            lifetime=timing_cell['lifetime'].dt.total_seconds().values
            time_maximum_heating=timing_cell['Heating_max'].dt.total_seconds().values[0]
            time_axis=precip_sum['time_cell'].dt.total_seconds().values

            length=len(time_axis)
            
            if shiftby is 'initiation':
                shift=0
            elif shiftby is 'glaciation':
                if not np.isnan(time_freezing):        
                    shift=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
                else:
                    shift=length
            elif shiftby is 'maximum_heating':
                if not np.isnan(time_maximum_heating):        
                    shift=np.argwhere(abs(time_axis-time_maximum_heating)<1)[0][0]
                else:
                    shift=length
            
            pad_left=n_left-shift
            pad_right=n_right-(len(time_axis)-shift)
            cut_left=0
            cut_right=0
            if pad_left<0:
                cut_left=-pad_left
                pad_left=0

            if pad_right<0:
                cut_right=-pad_right
                pad_right=0
            
            for precip_type in list_precip_types:
                precip_i=precip_sum[precip_type].values
                precip_i=precip_i[0+cut_left:precip_i.shape[0]-cut_right]
                precip_shifted=np.pad(precip_i,(pad_left,pad_right), 'constant', constant_values=(0, 0))
                Precipitation_shifted[precip_type].append(precip_shifted)
            
    for precip_type in list_precip_types:
        Precipitation_Composite_sum[precip_type]=np.sum(Precipitation_shifted[precip_type],axis=0)
        Precipitation_Composite_mean[precip_type]=np.mean(Precipitation_shifted[precip_type],axis=0)

    return Precipitation_Composite_sum, Precipitation_Composite_mean
