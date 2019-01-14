import numpy as np
import pandas as pd

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

def make_timings(track,threshold_value=None):
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
    
    timing_out['category']=None
    timing_out.at[timing_out['Freezing_time_relative']==0,'category']='nowarm'
    timing_out.at[timing_out['Freezing_first'].isnull(),'category']='warm'
    timing_out.at[(timing_out['Freezing_time_relative']>0) & (~timing_out['Freezing_first'].isnull()),'category']='warmcold'

    return timing_out


# def composite_Processes_slices(timing,model,case,cells='all',shiftby='initiation'):
#     list_processes=['Condensation','Evaporation','Freezing','Melting','Deposition','Sublimation','Rain formation']
#     from copy import deepcopy
#     if cells=='all':
#         cells=timing['cell'].values.astype(int)

#     cell_0=timing['cell'].values.astype(int)[0]
#     Processes_0=iris.load(os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}/Processes_lumped_along_{cell_0}.nc'))
    
#     Processes_base_Composite=iris.cube.CubeList()
#     n_left=100
#     n_right=100
#     n_total=n_left+n_right
    
#     Processes_shifted_along={}
#     Processes_shifted_across={}
    
#     for process in Processes_0:
#         Processes_shifted_along[process.name()]=[]
#         Processes_shifted_across[process.name()]=[]
#         process_composite=iris.cube.Cube(np.zeros((n_total,process.coord('geopotential_height').points.shape[0],process.coord('x_dx').points.shape[0])),
#                                           long_name=process.name(),units=process.units)
# #        process_composite=iris.cube.Cube(np.zeros((n_total,process.coord('model_level_number').points.shape[0],process.coord('x_dx').points.shape[0])),
# #                                         long_name=process.name(),units=process.units)

#         time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
#         process_composite.add_dim_coord(time_coord,data_dim=0)
#         process_composite.add_dim_coord(process.coord('geopotential_height'),data_dim=1)
# #        process_composite.add_dim_coord(process.coord('model_level_number'),data_dim=1)

#         process_composite.add_dim_coord(process.coord('x_dx'),data_dim=2)
#         Processes_base_Composite.append(process_composite)
        
#     Processes_along_mean_Composite=deepcopy(Processes_base_Composite)
#     Processes_across_mean_Composite=deepcopy(Processes_base_Composite)
#     Processes_along_sum_Composite=deepcopy(Processes_base_Composite)
#     Processes_across_sum_Composite=deepcopy(Processes_base_Composite)
    
#     #print(cells)
#     if list(cells):
#         for i_cell,cell in enumerate(cells):
# #             print('###############:')
# #             print('cell:',cell)
#             filename_along=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Processes_lumped_along_{cell}.nc')
#             filename_across=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Processes_lumped_across_{cell}.nc')
#             filename_track=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/track_processes_lumped_integrated_TWC_{cell}.h5')
#             if (glob.glob(filename_along) and glob.glob(filename_across) and glob.glob(filename_track)):
#                 Processes_along=iris.load(filename_along)
#                 Processes_across=iris.load(filename_across)

#                 Track_Processes_integrated=pd.read_hdf(filename_track,'table')

#                 timing_cell=timing[timing['cell']==cell]
#                 time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
#                 lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
#                 time_axis=Track_Processes_integrated['time_cell'].dt.total_seconds().values
#                 length=len(time_axis)
                
#                 if shiftby is 'initiation':
#                     shift=0
#                 elif shiftby is 'glaciation':
#                     if not np.isnan(time_freezing):        
#                         index_freezing=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
#                         shift=index_freezing
#                     else:
#                         shift=length
#                         index_freezing=np.nan
#                 pad_left=n_left-shift
#                 pad_right=n_right-(len(time_axis)-shift)
#                 cut_left=0
#                 cut_right=0
#                 if pad_left<0:
#                     cut_left=-pad_left
#                     pad_left=0

#                 if pad_right<0:
#                     cut_right=-pad_right
#                     pad_right=0

#                 for process in Processes_along_sum_Composite:
#                     process_along_i=Processes_along.extract_strict(process.name()).data
#                     process_along_i=process_along_i[0+cut_left:process_along_i.shape[0]-cut_right]
#                     process_along_i[process_along_i>1e20]=0
#                     process_along_shifted=np.pad(process_along_i,((pad_left,pad_right),(0,0),(0,0)),
#                                                  'constant', constant_values=((0, 0),(0,0),(0,0)))
#                     Processes_shifted_along[process.name()].append(process_along_shifted)

#                     process_across_i=Processes_across.extract_strict(process.name()).data
#                     process_across_i=process_across_i[0+cut_left:process_across_i.shape[0]-cut_right]
#                     process_across_i[process_across_i>1e20]=0
#                     process_across_shifted=np.pad(process_across_i,((pad_left,pad_right),(0,0),(0,0)), 
#                                                   'constant', constant_values=((0, 0),(0,0),(0,0)))
#                     Processes_shifted_across[process.name()].append(process_across_shifted)

#         for process in Processes_along_sum_Composite:
#             process.data=np.sum(Processes_shifted_along[process.name()],axis=0)

#         for process in Processes_across_sum_Composite:
#             process.data=np.sum(Processes_shifted_across[process.name()],axis=0)

#         for process in Processes_along_mean_Composite:
#             process.data=np.mean(Processes_shifted_along[process.name()],axis=0)
        
#         for process in Processes_across_mean_Composite:
#             process.data=np.mean(Processes_shifted_across[process.name()],axis=0)



#     return Processes_along_sum_Composite,Processes_across_sum_Composite,Processes_along_mean_Composite,Processes_across_mean_Composite

# def composite_Hydrometeors_mass_slices(timing,model,case,cells='all',shiftby='initiation'):

#     from copy import deepcopy
#     if cells=='all':
#         cells=timing['cell'].values.astype(int)

#     cell_0=timing['cell'].values.astype(int)[0]
#     Hydrometeors_0=iris.load(os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}/Hydrometeors_mass_along_{cell_0}.nc'))
    
#     Hydrometeors_base_Composite=iris.cube.CubeList()
#     Hydrometeors_base_Composite=iris.cube.CubeList()

#     n_left=100
#     n_right=100
#     n_total=n_left+n_right
    
#     Hydrometeors_shifted_along={}
#     Hydrometeors_shifted_across={}
    
#     for hydrometeor in  Hydrometeors_0:
#         Hydrometeors_shifted_along[hydrometeor.name()]=[]
#         Hydrometeors_shifted_across[hydrometeor.name()]=[]
#         hydrometeor_composite=iris.cube.Cube(np.zeros((n_total,hydrometeor.coord('geopotential_height').points.shape[0],hydrometeor.coord('x_dx').points.shape[0])),
#                                              long_name=hydrometeor.name(),units=hydrometeor.units)
# #        hydrometeor_composite=iris.cube.Cube(np.zeros((n_total,hydrometeor.coord('model_level_number').points.shape[0],hydrometeor.coord('x_dx').points.shape[0])),
# #                                         long_name=hydrometeor.name(),units=hydrometeor.units)

#         time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
#         hydrometeor_composite.add_dim_coord(time_coord,data_dim=0)
#         hydrometeor_composite.add_dim_coord(hydrometeor.coord('geopotential_height'),data_dim=1)
# #        hydrometeor_composite.add_dim_coord(hydrometeor.coord('model_level_number'),data_dim=1)

#         hydrometeor_composite.add_dim_coord(hydrometeor.coord('x_dx'),data_dim=2)
#         Hydrometeors_base_Composite.append(hydrometeor_composite)
        
#     Hydrometeors_along_mean_Composite=deepcopy(Hydrometeors_base_Composite)
#     Hydrometeors_across_mean_Composite=deepcopy(Hydrometeors_base_Composite)
#     Hydrometeors_along_sum_Composite=deepcopy(Hydrometeors_base_Composite)
#     Hydrometeors_across_sum_Composite=deepcopy(Hydrometeors_base_Composite)
    
#     #print(cells)
#     if list(cells):
#         for i_cell,cell in enumerate(cells):
# #             print('###############:')
# #             print('cell:',cell)
#             filename_along=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Hydrometeors_mass_along_{cell}.nc')
#             filename_across=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Hydrometeors_mass_across_{cell}.nc')
#             filename_track=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/track_processes_lumped_integrated_TWC_{cell}.h5')
#             if (glob.glob(filename_along) and glob.glob(filename_across) and glob.glob(filename_track)):
#                 Hydrometeors_along=iris.load(filename_along)
#                 Hydrometeors_across=iris.load(filename_across)

#                 Track_Hydrometeors_integrated=pd.read_hdf(filename_track,'table')

#                 timing_cell=timing[timing['cell']==cell]
#                 time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
#                 lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
#                 time_axis=Track_Hydrometeors_integrated['time_cell'].dt.total_seconds().values
#         #         print('time_freezing:',time_freezing)
#         #         print('lifetime:',lifetime)
#         #         print('time_axis:',time_axis)

#                 length=len(time_axis)
#         #         print('length:',length)            

#                 if shiftby is 'initiation':
#                     shift=0
#                 elif shiftby is 'glaciation':
#                     if not np.isnan(time_freezing):        
#                         index_freezing=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
#                         shift=index_freezing
#                     else:
#                         shift=length
#                         index_freezing=np.nan
                        
#                 pad_left=n_left-shift
#                 pad_right=n_right-(len(time_axis)-shift)
#                 cut_left=0
#                 cut_right=0
#                 if pad_left<0:
#                     cut_left=-pad_left
#                     pad_left=0

#                 if pad_right<0:
#                     cut_right=-pad_right
#                     pad_right=0

#                 for hydrometeor in Hydrometeors_along_sum_Composite:
#                     hydrometeor_along_i=Hydrometeors_along.extract_strict(hydrometeor.name()).data
#                     hydrometeor_along_i=hydrometeor_along_i[0+cut_left:hydrometeor_along_i.shape[0]-cut_right]
#                     hydrometeor_along_i[hydrometeor_along_i>1e20]=0
#                     hydrometeor_along_shifted=np.pad(hydrometeor_along_i,((pad_left,pad_right),(0,0),(0,0)),
#                                                  'constant', constant_values=((0, 0),(0,0),(0,0)))
#                     Hydrometeors_shifted_along[hydrometeor.name()].append(hydrometeor_along_shifted)

#                     hydrometeor_across_i=Hydrometeors_across.extract_strict(hydrometeor.name()).data
#                     hydrometeor_across_i=hydrometeor_across_i[0+cut_left:hydrometeor_across_i.shape[0]-cut_right]
#                     hydrometeor_across_i[hydrometeor_across_i>1e20]=0
#                     hydrometeor_across_shifted=np.pad(hydrometeor_across_i,((pad_left,pad_right),(0,0),(0,0)), 
#                                                   'constant', constant_values=((0, 0),(0,0),(0,0)))
#                     Hydrometeors_shifted_across[hydrometeor.name()].append(hydrometeor_across_shifted)

#         for hydrometeor in Hydrometeors_along_sum_Composite:
#             hydrometeor.data=np.sum(Hydrometeors_shifted_along[hydrometeor.name()],axis=0)

#         for hydrometeor in Hydrometeors_across_sum_Composite:
#             hydrometeor.data=np.sum(Hydrometeors_shifted_across[hydrometeor.name()],axis=0)

#         for hydrometeor in Hydrometeors_along_mean_Composite:
#             hydrometeor.data=np.mean(Hydrometeors_shifted_along[hydrometeor.name()],axis=0)
        
#         for hydrometeor in Hydrometeors_across_mean_Composite:
#             hydrometeor.data=np.mean(Hydrometeors_shifted_across[hydrometeor.name()],axis=0)



#     return Hydrometeors_along_sum_Composite,Hydrometeors_across_sum_Composite,Hydrometeors_along_mean_Composite,Hydrometeors_across_mean_Composite


# def composite_Hydrometeors_number_slices(timing,model,case,cells='all',shiftby='initiation'):

#     from copy import deepcopy
#     if cells=='all':
#         cells=timing['cell'].values.astype(int)

#     cell_0=timing['cell'].values.astype(int)[0]
#     Hydrometeors_0=iris.load(os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}/Hydrometeors_number_along_{cell_0}.nc'))
        
#     Hydrometeors_base_Composite=iris.cube.CubeList()
#     Hydrometeors_base_Composite=iris.cube.CubeList()

#     n_left=100
#     n_right=100
#     n_total=n_left+n_right
    
#     Hydrometeors_shifted_along={}
#     Hydrometeors_shifted_across={}
    
#     for hydrometeor in  Hydrometeors_0:
#         Hydrometeors_shifted_along[hydrometeor.name()]=[]
#         Hydrometeors_shifted_across[hydrometeor.name()]=[]
#         hydrometeor_composite=iris.cube.Cube(np.zeros((n_total,hydrometeor.coord('geopotential_height').points.shape[0],hydrometeor.coord('x_dx').points.shape[0])),
#                                              long_name=hydrometeor.name(),units=hydrometeor.units)
# #        hydrometeor_composite=iris.cube.Cube(np.zeros((n_total,hydrometeor.coord('model_level_number').points.shape[0],hydrometeor.coord('x_dx').points.shape[0])),
# #                                         long_name=hydrometeor.name(),units=hydrometeor.units)

#         time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
#         hydrometeor_composite.add_dim_coord(time_coord,data_dim=0)
#         hydrometeor_composite.add_dim_coord(hydrometeor.coord('geopotential_height'),data_dim=1)
# #        hydrometeor_composite.add_dim_coord(hydrometeor.coord('model_level_number'),data_dim=1)

#         hydrometeor_composite.add_dim_coord(hydrometeor.coord('x_dx'),data_dim=2)
#         Hydrometeors_base_Composite.append(hydrometeor_composite)
        
#     Hydrometeors_along_mean_Composite=deepcopy(Hydrometeors_base_Composite)
#     Hydrometeors_across_mean_Composite=deepcopy(Hydrometeors_base_Composite)
#     Hydrometeors_along_sum_Composite=deepcopy(Hydrometeors_base_Composite)
#     Hydrometeors_across_sum_Composite=deepcopy(Hydrometeors_base_Composite)
    
#     #print(cells)
#     if list(cells):
#         for i_cell,cell in enumerate(cells):
# #             print('###############:')
# #             print('cell:',cell)
#             filename_along=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Hydrometeors_number_along_{cell}.nc')
#             filename_across=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Hydrometeors_number_across_{cell}.nc')
#             filename_track=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/track_processes_lumped_integrated_TWC_{cell}.h5')
#             if (glob.glob(filename_along) and glob.glob(filename_across) and glob.glob(filename_track)):
#                 Hydrometeors_along=iris.load(filename_along)
#                 Hydrometeors_across=iris.load(filename_across)

#                 Track_Hydrometeors_integrated=pd.read_hdf(filename_track,'table')

#                 timing_cell=timing[timing['cell']==cell]
#                 time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
#                 lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
#                 time_axis=Track_Hydrometeors_integrated['time_cell'].dt.total_seconds().values
#         #         print('time_freezing:',time_freezing)
#         #         print('lifetime:',lifetime)
#         #         print('time_axis:',time_axis)

#                 length=len(time_axis)
#         #         print('length:',length)            
#                 if shiftby is 'initiation':
#                     shift=0
#                 elif shiftby is 'glaciation':
#                     if not np.isnan(time_freezing):        
#                         index_freezing=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
#                         shift=index_freezing
#                     else:
#                         shift=length
#                         index_freezing=np.nan
#                 pad_left=n_left-shift
#                 pad_right=n_right-(len(time_axis)-shift)
#                 cut_left=0
#                 cut_right=0
#                 if pad_left<0:
#                     cut_left=-pad_left
#                     pad_left=0

#                 if pad_right<0:
#                     cut_right=-pad_right
#                     pad_right=0

#                 for hydrometeor in Hydrometeors_along_sum_Composite:
#                     hydrometeor_along_i=Hydrometeors_along.extract_strict(hydrometeor.name()).data
#                     hydrometeor_along_i=hydrometeor_along_i[0+cut_left:hydrometeor_along_i.shape[0]-cut_right]
#                     hydrometeor_along_i[hydrometeor_along_i>1e20]=0
#                     hydrometeor_along_shifted=np.pad(hydrometeor_along_i,((pad_left,pad_right),(0,0),(0,0)),
#                                                  'constant', constant_values=((0, 0),(0,0),(0,0)))
#                     Hydrometeors_shifted_along[hydrometeor.name()].append(hydrometeor_along_shifted)

#                     hydrometeor_across_i=Hydrometeors_across.extract_strict(hydrometeor.name()).data
#                     hydrometeor_across_i=hydrometeor_across_i[0+cut_left:hydrometeor_across_i.shape[0]-cut_right]
#                     hydrometeor_across_i[hydrometeor_across_i>1e20]=0
#                     hydrometeor_across_shifted=np.pad(hydrometeor_across_i,((pad_left,pad_right),(0,0),(0,0)), 
#                                                   'constant', constant_values=((0, 0),(0,0),(0,0)))
#                     Hydrometeors_shifted_across[hydrometeor.name()].append(hydrometeor_across_shifted)

#         for hydrometeor in Hydrometeors_along_sum_Composite:
#             hydrometeor.data=np.sum(Hydrometeors_shifted_along[hydrometeor.name()],axis=0)

#         for hydrometeor in Hydrometeors_across_sum_Composite:
#             hydrometeor.data=np.sum(Hydrometeors_shifted_across[hydrometeor.name()],axis=0)

#         for hydrometeor in Hydrometeors_along_mean_Composite:
#             hydrometeor.data=np.mean(Hydrometeors_shifted_along[hydrometeor.name()],axis=0)
        
#         for hydrometeor in Hydrometeors_across_mean_Composite:
#             hydrometeor.data=np.mean(Hydrometeors_shifted_across[hydrometeor.name()],axis=0)



#     return Hydrometeors_along_sum_Composite,Hydrometeors_across_sum_Composite,Hydrometeors_along_mean_Composite,Hydrometeors_across_mean_Composite

# def composite_Processes_time(timing,model,case,cells='all',shiftby='initiation'):
#     list_processes=['Condensation','Evaporation','Freezing','Melting','Deposition','Sublimation','Rain formation']
#     from copy import deepcopy
#     if cells=='all':
#         cells=timing['cell'].values.astype(int)

#     cell_0=timing['cell'].values.astype(int)[0]
        
#     Processes_sum_0=iris.load(os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}/Processes_lumped_profile_w_TWC_{cell_0}.nc'))
#     Processes_sum_0=Processes_sum_0.extract(list_processes)

#     Processes_sum_Composite=iris.cube.CubeList()
    
#     n_left=100
#     n_right=100
#     n_total=n_left+n_right
    
#     Processes_shifted={}
#     for process in Processes_sum_0:
#         Processes_shifted[process.name()]=[]
#         process_composite=iris.cube.Cube(np.zeros((n_total,process.shape[1])),long_name=process.name(),units=process.units)
#         time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
#         process_composite.add_dim_coord(time_coord,data_dim=0)
#         process_composite.add_dim_coord(process.coord('geopotential_height'),data_dim=1)
#         Processes_sum_Composite.append(process_composite)
#     Processes_mean_Composite=deepcopy(Processes_sum_Composite)
    
#     #print(cells)
#     if list(cells):
#         for i_cell,cell in enumerate(cells):
#     #         print('###############:')
#     #         print('cell:',cell)
#             filename_profile=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Processes_lumped_profile_TWC_{cell}.nc')
#             filename_track=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/track_processes_lumped_integrated_TWC_{cell}.h5')
#             if (glob.glob(filename_profile) and glob.glob(filename_track)):
#                 Processes_sum=iris.load(filename_profile)
#                 Processes_sum=Processes_sum.extract(list_processes)

#                 Track_Processes_integrated=pd.read_hdf(filename_track,'table')

#                 timing_cell=timing[timing['cell']==cell]
#                 time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
#                 lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
#                 time_axis=Track_Processes_integrated['time_cell'].dt.total_seconds().values

#                 length=len(time_axis)

#                 if shiftby is 'initiation':
#                     shift=0
#                 elif shiftby is 'glaciation':
#                     if not np.isnan(time_freezing):        
#                         index_freezing=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
#                         shift=index_freezing
#                     else:
#                         shift=length
#                         index_freezing=np.nan
#                 pad_left=n_left-shift
#                 pad_right=n_right-(len(time_axis)-shift)
#                 cut_left=0
#                 cut_right=0
#                 if pad_left<0:
#                     cut_left=-pad_left
#                     pad_left=0

#                 if pad_right<0:
#                     cut_right=-pad_right
#                     pad_right=0

#                 for process in Processes_sum_Composite:
#                     process_i=Processes_sum.extract_strict(process.name()).data
#                     process_i=process_i[0+cut_left:process_i.shape[0]-cut_right]
#                     process_i[process_i>1e20]=0
#                     process_shifted=np.pad(process_i,((pad_left,pad_right),(0,0)), 'constant', constant_values=((0, 0),(0,0)))
#                     Processes_shifted[process.name()].append(process_shifted)

#         for process in Processes_sum_Composite:
#             process.data=np.sum(Processes_shifted[process.name()],axis=0)

#         for process in Processes_mean_Composite:
#             process.data=np.mean(Processes_shifted[process.name()],axis=0)

#     return Processes_sum_Composite,Processes_mean_Composite


# def composite_Hydrometeors_time(timing,model,case,cells='all',shiftby='initiation'):
#     from copy import deepcopy
#     if cells=='all':
#         cells=timing['cell'].values.astype(int)

#     cell_0=timing['cell'].values.astype(int)[0]
        
#     Hydrometeors_sum_0=iris.load(os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}/Hydrometeors_profile_TWC_{cell_0}.nc'))

#     Hydrometeors_sum_Composite=iris.cube.CubeList()
    
#     n_left=100
#     n_right=100
#     n_total=n_left+n_right
    
#     Hydrometeors_shifted={}
#     for hydrometeor in Hydrometeors_sum_0:
#         Hydrometeors_shifted[hydrometeor.name()]=[]
#         hydrometeor_composite=iris.cube.Cube(np.zeros((n_total,hydrometeor.shape[1])),long_name=hydrometeor.name(),units=hydrometeor.units)
#         time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
#         hydrometeor_composite.add_dim_coord(time_coord,data_dim=0)
#         hydrometeor_composite.add_dim_coord(hydrometeor.coord('geopotential_height'),data_dim=1)
#         Hydrometeors_sum_Composite.append(hydrometeor_composite)
#     Hydrometeors_mean_Composite=deepcopy(Hydrometeors_sum_Composite)
    
#     if list(cells):
#         for i_cell,cell in enumerate(cells):
#             filename_profile=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/Hydrometeors_profile_TWC_{cell}.nc')
#             filename_track=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}/track_hydrometeros_integrated_TWC_{cell}.h5')
#             if (glob.glob(filename_profile) and glob.glob(filename_track)):
#                 Hydrometeors_sum=iris.load(filename_profile)

#                 Track_Hydrometeors_integrated=pd.read_hdf(filename_track,'table')

#                 timing_cell=timing[timing['cell']==cell]
#                 time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
#                 lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
#                 time_axis=Track_Hydrometeors_integrated['time_cell'].dt.total_seconds().values

#                 length=len(time_axis)

#                 if not np.isnan(time_freezing):        
#                     index_freezing=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
#                     shift=index_freezing
#                 else:
#                     shift=length
#                     index_freezing=np.nan
#                 pad_left=n_left-shift
#                 pad_right=n_right-(len(time_axis)-shift)
#                 cut_left=0
#                 cut_right=0
#                 if pad_left<0:
#                     cut_left=-pad_left
#                     pad_left=0

#                 if pad_right<0:
#                     cut_right=-pad_right
#                     pad_right=0

#                 for hydrometeor in Hydrometeors_sum_Composite:
#                     hydrometeor_i=Hydrometeors_sum.extract_strict(hydrometeor.name()).data
#                     hydrometeor_i=hydrometeor_i[0+cut_left:hydrometeor_i.shape[0]-cut_right]
#                     hydrometeor_i[hydrometeor_i>1e20]=0
#                     hydrometeor_shifted=np.pad(hydrometeor_i,((pad_left,pad_right),(0,0)), 'constant', constant_values=((0, 0),(0,0)))
#                     Hydrometeors_shifted[hydrometeor.name()].append(hydrometeor_shifted)

#         for hydrometeor in Hydrometeors_sum_Composite:
#             hydrometeor.data=np.sum(Hydrometeors_shifted[hydrometeor.name()],axis=0)

#         for hydrometeor in Hydrometeors_mean_Composite:
#             hydrometeor.data=np.mean(Hydrometeors_shifted[hydrometeor.name()],axis=0)

#     return Hydrometeors_sum_Composite,Hydrometeors_mean_Composite


# def composite_Precipitation_time(timing,model,case,cells='all',shiftby='initiation'):
#     from copy import deepcopy
#     if cells=='all':
#         cells=timing['cell'].values.astype(int)

#     cell_0=timing['cell'].values.astype(int)[0]
    
#     list_precip_types= ['surface_precipitation','surface_precipitation_accumulated','surface_precipitation_instantaneous']
    
#     Precipitation_0=pd.read_hdf(os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell_0}/track_precipitation_TWC_{cell_0}.h5'),'table')
# #    Hydrometeors_sum_0=Hydrometeors_sum_0.extract(list_Hydrometeors)

#     Precipitation_Composite_sum=pd.DataFrame()
#     Precipitation_shifted={}

#     n_left=100
#     n_right=100
#     n_total=n_left+n_right
#     for precip_type in list_precip_types:
#         Precipitation_shifted[precip_type]=[]
#         Precipitation_Composite_sum[precip_type]=np.zeros((n_total))
#         Precipitation_Composite_sum['time']=np.arange(-n_left,n_right)
#         #precip_out=Precipitation_Composite_sum[precip_type].values
#         #precip=iris.cube.Cube(np.zeros((n_total,hydrometeor.shape[1])),long_name=hydrometeor.name(),units=hydrometeor.units)
#         #time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
#         #hydrometeor_composite.add_dim_coord(time_coord,data_dim=0)
#         #hydrometeor_composite.add_dim_coord(hydrometeor.coord('geopotential_height'),data_dim=1)
#         #Hydrometeors_sum_Composite.append(hydrometeor_composite)
#     Precipitation_Composite_mean=deepcopy(Precipitation_Composite_sum)
#     for i_cell,cell in enumerate(cells):
# #        print('###############:')
# #        print('cell:',cell)
#         savedir_i=os.path.join(acpc_workspace,f'Analysis/mheiken/Tracking/Save_1min/{version}/{model}/{case}/cells/{cell}')
# #        filename_profile=os.path.join(savedir_i,f'Hydrometeors_profile_TWC_{cell}.nc')
#         filename_track=os.path.join(savedir_i,f'track_precipitation_TWC_{cell}.h5')
#         if (glob.glob(filename_track)):
#             #precip_sum=iris.load(filename_profile)
#             precip_sum=pd.read_hdf(filename_track,'table')
#             timing_cell=timings[model][case][timings[model][case]['cell']==cell]
#             time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values
#             lifetime=timing_cell['lifetime'].dt.total_seconds().values
#             time_axis=precip_sum['time_cell'].dt.total_seconds().values

#             length=len(time_axis)
# #            print(time_freezing)
#             if shiftby is 'initiation':
#                 shift=0
#             elif shiftby is 'glaciation':
#                 if not np.isnan(time_freezing):        
#                     index_freezing=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
#                     shift=index_freezing
#                 else:
#                     shift=length
#                     index_freezing=np.nan
            
#             pad_left=n_left-shift
#             pad_right=n_right-(len(time_axis)-shift)
#             cut_left=0
#             cut_right=0
#             if pad_left<0:
#                 cut_left=-pad_left
#                 pad_left=0

#             if pad_right<0:
#                 cut_right=-pad_right
#                 pad_right=0
            
#             for precip_type in list_precip_types:
#                 precip_i=precip_sum[precip_type].values
#                 precip_i=precip_i[0+cut_left:precip_i.shape[0]-cut_right]
#                 precip_shifted=np.pad(precip_i,(pad_left,pad_right), 'constant', constant_values=(0, 0))
#                 Precipitation_shifted[precip_type].append(precip_shifted)
            
#     for precip_type in list_precip_types:
#         Precipitation_Composite_sum[precip_type]=np.sum(Precipitation_shifted[precip_type],axis=0)
#         Precipitation_Composite_mean[precip_type]=np.mean(Precipitation_shifted[precip_type],axis=0)

#     return Precipitation_Composite_sum, Precipitation_Composite_mean