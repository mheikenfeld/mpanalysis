import os
import iris
import pandas as pd
import glob
import numpy as np
from mpdiag import water_content_from_hydrometeor

def load_track_processes_mass(savedir,microphysics_scheme):
    '''
    Load data from cellwise preprocessed data for calculations around the composites and categorisation
    '''
    list_processes=['Condensation','Evaporation','Freezing','Melting','Deposition','Sublimation','Rain formation']
    Track=pd.read_hdf(os.path.join(savedir,'Track.h5'),'table')
    Track_filtered=Track
    cells=Track_filtered['cell'].dropna().unique()

    Track_processes=pd.DataFrame()
    Track_hydrometeors=pd.DataFrame()
    for cell in cells:
        filename=glob.glob(os.path.join(savedir,f'cells/{int(cell)}/track_processes_lumped_integrated_TWC.h5'))
        if filename:
            track_cell_processes=pd.read_hdf(filename[0],'table')
            Track_processes=Track_processes.append(track_cell_processes,ignore_index=False)
            for process in list_processes:
                Track_processes[process]=abs(Track_processes[process])
        filename=glob.glob(os.path.join(savedir,f'cells/{int(cell)}/track_hydrometeors_integrated_TWC.h5'))
        if filename:
            track_cell_hydrometeors=pd.read_hdf(filename[0],'table')
            Track_hydrometeors=Track_hydrometeors.append(track_cell_hydrometeors)
    Track_hydrometeors['total_water_mass']=water_content_from_hydrometeor(Track_hydrometeors,type='total',microphysics_scheme=microphysics_scheme)
    Track_hydrometeors['liquid_water_mass']=water_content_from_hydrometeor(Track_hydrometeors,type='liquid',microphysics_scheme=microphysics_scheme)
    Track_hydrometeors['frozen_water_mass']=water_content_from_hydrometeor(Track_hydrometeors,type='frozen',microphysics_scheme=microphysics_scheme)
    Track_hydrometeors['frozen_water_fraction']=Track_hydrometeors['frozen_water_mass']/Track_hydrometeors['total_water_mass']

    Track_all=Track_processes.merge(Track_hydrometeors[['feature','total_water_mass','liquid_water_mass','frozen_water_mass','frozen_water_fraction']],on='feature',how='left')
    #Track_all=Track_processes
    return Track_all


def get_timing(track,process=None,measure=None,value=None):
    '''
    Calculate timings of important stages in the evolution of the tracked cells
    '''
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
            print(index)
            print(type(index))
            # reject a maximum smaller than the threshold value
            if np.isnan(index):
                index=None
            elif track_cell.loc[index,process]<value:
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
    return df_out

def latent_heating(track):
    '''
    Calculate total latent heating from microphysical process rates
    '''
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
    '''
    Calculate the timing of specific stages of the cloud evolution
    '''
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
    '''
    Evaluate properties at specific times of the evolution of the cloud for further analyses
    '''
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
    '''
    Categorise clouds into different types
    '''
    timing_out=timing_in
    category_out=timing_in[['cell','time_start','time_end']].copy()
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

def timing_index(timing_cell,track,shiftby):
    time_freezing=timing_cell['Freezing_first'].dt.total_seconds().values[0]
    time_maximum_heating=timing_cell['Heating_max'].dt.total_seconds().values[0]
    time_maximum_mass=timing_cell['total_max'].dt.total_seconds().values[0]
#         time_maximum_rain=timing_cell['Rain_max'].dt.total_seconds().values[0]

    lifetime=timing_cell['lifetime'].dt.total_seconds().values[0]
    time_axis=track['time_cell'].dt.total_seconds().values
        
    length=len(time_axis)

    if shiftby is 'initiation':
        shift=0
    elif shiftby is 'glaciation':
        if not np.isnan(time_freezing):        
            shift=np.argwhere(abs(time_axis-time_freezing)<1)[0][0]
        else:
            shift=length-1
    elif shiftby is 'maximum_heating':
        if not np.isnan(time_maximum_heating):        
            shift=np.argwhere(abs(time_axis-time_maximum_heating)<1)[0][0]
        else:
            shift=length-1

    elif shiftby is 'maximum_mass':
        if not np.isnan(time_maximum_heating):        
            shift=np.argwhere(abs(time_axis-time_maximum_heating)<1)[0][0]
        else:
            shift=length-1
    else:
        raise ValueError(f'shiftby {shiftby} unknown, must be initiation,glaciation,maximum_heating,maximum_mass or maximum_rain')
    return shift


def make_shift(timing_cell,track,shiftby,n_left,n_right):
    '''
    Shift individual clouds for cloud composites
    '''
    shift=timing_index(timing_cell,track,shiftby)
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
    return shift, pad_left, pad_right, cut_left, cut_right

def composite_cells_time(cell_dict,track_dict,timing,shiftby='initiation',n_left=100,n_right=100):
    '''
    Create cloud composites for a dictionary containing individual vertically summed cloud profiles
    '''
    from copy import deepcopy
    n_total=n_left+n_right
    dim_shift=0
    cell_dict_shifted={}
    cubes_shifted={}
    cubelist_sum_Composite=iris.cube.CubeList()
    cubelist_mean_Composite=iris.cube.CubeList()
    cells=list(cell_dict.keys())
    for cube in cell_dict[cells[0]]:
        cubes_shifted[cube.name()]=[]
        cube_composite=iris.cube.Cube(np.zeros((n_total,cube.shape[1])),long_name=cube.name(),units=cube.units)
        time_coord=iris.coords.DimCoord(np.arange(-n_left,n_right),long_name='time',units='minute')
        cube_composite.add_dim_coord(time_coord,data_dim=0)
        cube_composite.add_dim_coord(cube.coord('geopotential_height'),data_dim=1)
        cubelist_sum_Composite.append(cube_composite)
    cubelist_mean_Composite=deepcopy(cubelist_sum_Composite)
    
    for cell in cells:
        cubes_in=cell_dict[cell]
        track_integrated=track_dict[cell]
        shift, pad_left, pad_right, cut_left, cut_right=make_shift(timing_cell=timing.loc[timing['cell']==cell],track=track_integrated,
                                                                       shiftby=shiftby,
                                                                       n_left=n_left,n_right=n_right)        
        for cube in cubelist_sum_Composite:
            cube_i=cubes_in.extract_strict(cube.name()).data
            cube_i=cube_i[0+cut_left:cube_i.shape[0]-cut_right]
            # Why this?
            cube_i[cube_i>1e20]=0
            cube_shifted=np.pad(cube_i.data,((pad_left,pad_right),(0,0)), 'constant', constant_values=((0, 0),(0,0)))
            cubes_shifted[cube.name()].append(cube_shifted)
                
    for cube in cubelist_sum_Composite:
        cube.data=np.sum(cubes_shifted[cube.name()],axis=dim_shift)

    for cube in cubelist_mean_Composite:
        cube.data=np.mean(cubes_shifted[cube.name()],axis=dim_shift)

    return cubelist_sum_Composite,cubelist_mean_Composite

def composite_Processes_time(savedir,timing,cells='all',type_mask='_TWC',shiftby='initiation',n_left=100,n_right=100):
    '''
    Composite calculation for microphysical processes including loading data from files into dictionary
    '''
    from copy import deepcopy
    if cells=='all':
        cells=timing['cell'].values.astype(int)
    cell_dict={}
    track_dict={}
    
    for cell in cells:    
        savedir_cell=os.path.join(savedir,'cells',f'{cell}')    
        filename_profile=os.path.join(savedir_cell,f'Processes_lumped_profile{type_mask}.nc')
        filename_track=os.path.join(savedir_cell,f'track_processes_lumped_integrated{type_mask}.h5')
        if (glob.glob(filename_profile) and glob.glob(filename_track)):
            Processes_sum=iris.load(filename_profile)
            cell_dict[cell]=iris.load(filename_profile)
            track_dict[cell]=pd.read_hdf(filename_track,'table')

    Processes_sum_Composite,Processes_mean_Composite=composite_cells_time(cell_dict,track_dict,timing,shiftby=shiftby,n_left=n_left,n_right=n_right)
    return Processes_sum_Composite,Processes_mean_Composite

def composite_Hydrometeors_mass_time(savedir,timing,cells='all',type_mask='_TWC',shiftby='initiation',n_left=100,n_right=100):
    '''
    Composite calculation for hydrometeor mixing ratios including loading data from files into dictionary
    '''

    from copy import deepcopy
    if cells=='all':
        cells=timing['cell'].values.astype(int)
    cell_dict={}
    track_dict={}

    for cell in cells:  
        savedir_cell=os.path.join(savedir,'cells',f'{cell}')    
        filename_profile=os.path.join(savedir_cell,f'Hydrometeors_profile{type_mask}.nc')
        filename_track=os.path.join(savedir_cell,f'track_hydrometeors_integrated{type_mask}.h5')
        if (glob.glob(filename_profile) and glob.glob(filename_track)):
            Hydrometeors_sum=iris.load(filename_profile)

            cell_dict[cell]=iris.load(filename_profile)
            track_dict[cell]=pd.read_hdf(filename_track,'table')
    print('data loaded')

    Hydrometeors_sum_Composite,Hydrometeors_mean_Composite=composite_cells_time(cell_dict,track_dict,timing,shiftby=shiftby,n_left=n_left,n_right=n_right)
    return Hydrometeors_sum_Composite,Hydrometeors_mean_Composite

def composite_Precipitation_time_(savedir,timing,cells='all',type_mask='_TWC',shiftby='initiation',n_left=100,n_right=100):
    '''
    Composite calculation for microphysical processes including loading data from files into dictionary
    '''
    from copy import deepcopy
    
    if cells=='all':
        cells=timing['cell'].values.astype(int)

    cell_0=timing['cell'].values.astype(int)[0]
    
    list_precip_types= ['surface_precipitation','surface_precipitation_accumulated']#,'surface_precipitation_instantaneous']
    savedir_0=os.path.join(savedir,'cells',f'{cell_0}')
    Precipitation_0=pd.read_hdf(os.path.join(savedir_0,f'track_precipitation{type_mask}.h5'),'table')
    print(Precipitation_0)

    Precipitation_Composite_sum=pd.DataFrame()
    Precipitation_shifted={}

    n_total=n_left+n_right
    for precip_type in list_precip_types:
        Precipitation_shifted[precip_type]=[]
        Precipitation_Composite_sum[precip_type]=np.zeros((n_total))
        Precipitation_Composite_sum['time']=np.arange(-n_left,n_right)
    Precipitation_Composite_mean=deepcopy(Precipitation_Composite_sum)
    for i_cell,cell in enumerate(cells):
        savedir_i=os.path.join(savedir,'cells',f'{cell}')
#        filename_profile=os.path.join(savedir_i,f'Hydrometeors_profile_TWC_{cell}.nc')
        filename_track=os.path.join(savedir_i,f'track_precipitation{type_mask}.h5')
        if (glob.glob(filename_track)):
            precip_sum=pd.read_hdf(filename_track,'table')
            shift, pad_left, pad_right, cut_left, cut_right=make_shift(timing_cell=timing.loc[timing['cell']==cell],track=precip_sum,
                                                                           shiftby=shiftby,
                                                                           n_left=n_left,n_right=n_right)
            
            for precip_type in list_precip_types:
                precip_i=precip_sum[precip_type].values
                precip_i=precip_i[0+cut_left:precip_i.shape[0]-cut_right]
                precip_shifted=np.pad(precip_i,(pad_left,pad_right), 'constant', constant_values=(0, 0))
                Precipitation_shifted[precip_type].append(precip_shifted)
            
    for precip_type in list_precip_types:
        Precipitation_Composite_sum[precip_type]=np.sum(Precipitation_shifted[precip_type],axis=0)
        Precipitation_Composite_mean[precip_type]=np.mean(Precipitation_shifted[precip_type],axis=0)

    return Precipitation_Composite_sum, Precipitation_Composite_mean

def composite_cells_slice(cell_dict,track_dict,timing,shiftby='initiation'):
    '''
    Create cloud composites for a dictionary containing individual vertically summed cloud profiles
    '''
    from copy import deepcopy    
    dim_shift=0
    cell_dict_shifted={}
    cubes_shifted={}
    cubelist_sum_Composite=iris.cube.CubeList()
    cubelist_mean_Composite=iris.cube.CubeList()
    cells=list(cell_dict.keys())
    for cube in cell_dict[cells[0]]:
        cubes_shifted[cube.name()]=[]
        cube_composite=iris.cube.Cube(np.zeros((cube.shape[1],cube.shape[2])),long_name=cube.name(),units=cube.units)
        cube_composite.add_dim_coord(cube.coord('geopotential_height'),0)
        cube_composite.add_dim_coord(cube.coord('x_dx'),1)

        cubelist_sum_Composite.append(cube_composite)
    cubelist_mean_Composite=deepcopy(cubelist_sum_Composite)
    
    for cell in cells:
        cubes_in=cell_dict[cell]
        track_integrated=track_dict[cell]
        shift=timing_index(timing_cell=timing.loc[timing['cell']==cell],track=track_integrated,shiftby=shiftby)  
        for cube in cubelist_sum_Composite:
            cube_i=cubes_in.extract_strict(cube.name()).data
            cube_shifted=cube_i[shift]
            cubes_shifted[cube.name()].append(cube_shifted)
                
    for cube in cubelist_sum_Composite:
        cube.data=np.sum(cubes_shifted[cube.name()],axis=dim_shift)

    for cube in cubelist_mean_Composite:
        cube.data=np.mean(cubes_shifted[cube.name()],axis=dim_shift)

    return cubelist_sum_Composite,cubelist_mean_Composite

def composite_Processes_slice(savedir,timing,cells='all',shiftby='initiation'):
    '''
    Composite calculation for microphysical processes including loading data from files into dictionary
    '''
    from copy import deepcopy
    if cells=='all':
        cells=timing['cell'].values.astype(int)
    cell_dict_along={}
    cell_dict_across={}
    track_dict={}
    
    for cell in cells:
        savedir_cell=os.path.join(savedir,'cells',f'{cell}')    
        filename_along=os.path.join(savedir_cell,f'Processes_lumped_along.nc')
        filename_across=os.path.join(savedir_cell,f'Processes_lumped_across.nc')

        filename_track=os.path.join(savedir_cell,f'track_processes_lumped_integrated_TWC.h5')
        if (glob.glob(filename_along) and glob.glob(filename_across) and glob.glob(filename_track)):
            cell_dict_along[cell]=iris.load(filename_along)
            cell_dict_across[cell]=iris.load(filename_across)
            track_dict[cell]=pd.read_hdf(filename_track,'table')

    Processes_sum_Composite_along,Processes_mean_Composite_along=composite_cells_slice(cell_dict_along,track_dict,timing,shiftby=shiftby)
    Processes_sum_Composite_across,Processes_mean_Composite_across=composite_cells_slice(cell_dict_across,track_dict,timing,shiftby=shiftby)
    return Processes_sum_Composite_along,Processes_mean_Composite_along,Processes_sum_Composite_across,Processes_mean_Composite_across


def composite_Hydrometeors_mass_slice(savedir,timing,cells='all',shiftby='initiation'):
    '''
    Composite calculation for microphysical processes including loading data from files into dictionary
    '''
    from copy import deepcopy
    if cells=='all':
        cells=timing['cell'].values.astype(int)
    cell_dict_along={}
    cell_dict_across={}
    track_dict={}
    
    for cell in cells:
        savedir_cell=os.path.join(savedir,'cells',f'{cell}')    
        filename_along=os.path.join(savedir_cell,f'Hydrometeors_number_along.nc')
        filename_across=os.path.join(savedir_cell,f'Hydrometeors_number_across.nc')

        filename_track=os.path.join(savedir_cell,f'track_hydrometeors_integrated_TWC.h5')
        if (glob.glob(filename_along) and glob.glob(filename_across) and glob.glob(filename_track)):
            cell_dict_along[cell]=iris.load(filename_along)
            cell_dict_across[cell]=iris.load(filename_across)
            track_dict[cell]=pd.read_hdf(filename_track,'table')

    Hydrometeors_sum_Composite_along,Hydrometeors_mean_Composite_along=composite_cells_slice(cell_dict_along,track_dict,timing,shiftby=shiftby)
    Hydrometeors_sum_Composite_across,Hydrometeors_mean_Composite_across=composite_cells_slice(cell_dict_across,track_dict,timing,shiftby=shiftby)
    return Hydrometeors_sum_Composite_along,Hydrometeors_mean_Composite_along,Hydrometeors_sum_Composite_across,Hydrometeors_mean_Composite_across
