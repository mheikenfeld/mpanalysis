import iris
import numpy as np
# import pandas as pd
# import os
# import matplotlib.pyplot as plt
# import iris.plot as iplt
# import iris.quickplot as qplt
import logging


# def my_interpolate_3D2D_altitude(cube_in,cube_grid,coordinates_in=None,coordinates_out=None,method='linear'):
#     import numpy as np
#     from scipy.interpolate import RegularGridInterpolator,interp1d

#     #Calculate input grid for griddate method
#     source_points=[cube_in.coord(c).points for c in coordinates_in]
#     source_points_interp=tuple(source_points)

#     #Calculate output grid for griddate method
#     mesh_coord_0,mesh_coord_1=np.meshgrid(cube_in.coord(coordinates_in[0]).points,cube_grid.coord(coordinates_in[1]).points,indexing='ij')
#     mesh_coord_0,mesh_coord_2=np.meshgrid(cube_in.coord(coordinates_in[0]).points,cube_grid.coord(coordinates_in[2]).points,indexing='ij')
#     target_points_interp=np.stack(tuple([mesh_coord_0.flatten(),mesh_coord_1.flatten(),mesh_coord_2.flatten()])).T
#     #Calcululate interpolate output and bring it into the right shape:

#     Interpolator=RegularGridInterpolator(source_points_interp,cube_in.data,method=method,bounds_error=False,fill_value=0)

#     target_data=Interpolator(target_points_interp).reshape((cube_in.shape[0],cube_grid.coord(coordinates_out[1]).shape[0]))

#     Interpolator_altitude=RegularGridInterpolator(source_points_interp,cube_in.coord(coordinates_out[0]).points,method=method,bounds_error=False,fill_value=0)
#     altitude_data=Interpolator_altitude(target_points_interp).reshape((cube_in.shape[0],cube_grid.coord(coordinates_out[1]).shape[0]))

#     target_data_altitude=np.nan*np.ones(cube_grid.shape)
#     for i in range(target_data_altitude.shape[-1]):
#         target_data_altitude_interpolator=interp1d(altitude_data[:,i],target_data[:,i],bounds_error=False,fill_value=0)
#         target_data_altitude[...,i]=target_data_altitude_interpolator(cube_grid.coord(coordinates_out[0]).points)

#     #Store output in cube and adjust metadata:
#     cube_out=cube_grid.copy()
#     cube_out.data=target_data_altitude
#     cube_out.rename(cube_in.name())
#     cube_out.units=cube_in.units
#     return(cube_out)

# def my_interpolate_3D2D(cube_in,cube_grid,coordinates,method='linear'):
#     import numpy as np
#     from scipy.interpolate import RegularGridInterpolator
   
#     #Calculate input grid for griddate method
#     source_points=[cube_in.coord(c).points for c in coordinates]
#     source_points_interp=tuple(source_points)

#     #Calculate output grid for griddate method
#     mesh_coord_0,mesh_coord_1=np.meshgrid(cube_grid.coord(coordinates[0]).points,cube_grid.coord(coordinates[1]).points,indexing='ij')
#     mesh_coord_0,mesh_coord_2=np.meshgrid(cube_grid.coord(coordinates[0]).points,cube_grid.coord(coordinates[2]).points,indexing='ij')
#     target_points_interp=np.stack(tuple([mesh_coord_0.flatten(),mesh_coord_1.flatten(),mesh_coord_2.flatten()])).T

#     #Calcululate interpolate output and bring it into the right shape:
#     Interpolator=RegularGridInterpolator(source_points_interp,cube_in.data,method=method,bounds_error=False,fill_value=0)
#     target_data=Interpolator(target_points_interp).reshape(cube_grid.shape)

#     #Store output in cube and adjust metadata:
#     cube_out=cube_grid.copy()
#     cube_out.data=target_data
#     cube_out.rename(cube_in.name())
#     cube_out.units=cube_in.units

#     return(cube_out)

# def interpolate_to_cubes(cubelist_in,grid_along,grid_across,z_coord='model_level_number'):
#     from iris.cube import CubeList
#     cubelist_out_along= CubeList()
#     cubelist_out_across= CubeList()

#     for variable in cubelist_in:
        
        
#         if z_coord=='model_level_number':
#             cube_along=my_interpolate_3D2D(variable,grid_along,coordinates=['model_level_number','projection_y_coordinate','projection_x_coordinate'],method='linear')
#             cube_across=my_interpolate_3D2D(variable,grid_across,coordinates=['model_level_number','projection_y_coordinate','projection_x_coordinate'],method='linear')
#         elif z_coord=='geopotential_height':
#             cube_along=my_interpolate_3D2D_altitude(variable,grid_along,coordinates_in=['model_level_number','projection_y_coordinate','projection_x_coordinate'],coordinates_out=['geopotential_height','x_dx'],method='linear')
#             cube_across=my_interpolate_3D2D_altitude(variable,grid_across,coordinates_in=['model_level_number','projection_y_coordinate','projection_x_coordinate'],coordinates_out=['geopotential_height','x_dx'],method='linear')
#         else:
#             raise ValueError('z_coord must be model_level_number or geopotential_height ')
#         cubelist_out_along.append(cube_along)
#         cubelist_out_across.append(cube_across)
        
#         del cube_along
#         del cube_across
#     return cubelist_out_along,cubelist_out_across

# def make_cubes_alongacross(x,y,alpha,cube_sample,dx,width,height_levels=None,z_coord='model_level_number',interpolate_time=False):
#     """
#     Create cubes along and across the cell at the current timestep to use for interpolation to them from other variables 
#     Cube contains no data (zeros)
#         Parameters
#     ----------
#     x : float
#         x coordinate of the track
#     y : float
#         y coordinate of the track
#     alpha: float
#         angle betweeen track and x axis (range: -pi tp pi)
#     cube_sample: iris.cube.Cube
#                   sample cube of original data (for time and vertical coordinate)
#     dx     : float
#              horizontal spacing interpolated to (in m)
#     width  : int 
#              half width of the interpolated cube in (numbers of points)
#     height_levels: numpy ndarray
#                    vertical levels (geopotential_height) to interpolate to (in m)
            
                 
#     Returns
#     -------
    
#     grid_along   iris.cube.Cube
    
#     grid_across  iris.cube.Cube

#     """
#     time_coord=cube_sample.coord('time')
#     if z_coord=='model_level_number':
#         model_level_number_coord=iris.coords.DimCoord(cube_sample.coord('model_level_number')[::5].points,
#                                                       long_name='model_level_number',
#                                                       var_name='model_level_number',units=' ')
#     elif z_coord=='geopotential_height':
#         geopotential_height_coord=iris.coords.DimCoord(height_levels,long_name='geopotential_height',var_name='geopotential_height',units='m')
#     else:
#         raise ValueError('z_coord must be model_level_number or geopotential_height ')
       
      

#     # create new x coordinate (new cube)
#     x_dx=dx*(np.arange(-1.*width,1.*width+1))
#     x_dx_coord=iris.coords.DimCoord(x_dx,long_name='x_dx',var_name='x_dx',units='m')
        
#     #create coordinates x and y for new cube:  
#     x_along=x+np.expand_dims(x_dx,0)*np.cos(alpha)
#     y_along=y+np.expand_dims(x_dx,0)*np.sin(alpha)  
    
#     x_along_coord=iris.coords.AuxCoord(x_along,standard_name='projection_x_coordinate',units='m')
#     y_along_coord=iris.coords.AuxCoord(y_along,standard_name='projection_y_coordinate',units='m')
    
#     x_across=x+np.expand_dims(x_dx,0)*np.sin(alpha)
#     y_across=y-np.expand_dims(x_dx,0)*np.cos(alpha)    
#     x_across_coord=iris.coords.AuxCoord(x_across,standard_name='projection_x_coordinate',long_name='x',var_name='x',units='m')
#     y_across_coord=iris.coords.AuxCoord(y_across,standard_name='projection_y_coordinate',long_name='y',var_name='y',units='m')
    
#     if z_coord=='model_level_number':
#        grid_along=iris.cube.Cube(np.zeros((time_coord.shape[0],model_level_number_coord.shape[0],x_dx_coord.shape[0])),
#                                  long_name='cube_along',var_name='cube_along',
#                                  dim_coords_and_dims=[(time_coord,0),(model_level_number_coord, 1), (x_dx_coord, 2)],
#                                  aux_coords_and_dims=[(x_along_coord, [0,2]), (y_along_coord, [0,2])])
#        grid_across=iris.cube.Cube(np.zeros((time_coord.shape[0],model_level_number_coord.shape[0],x_dx_coord.shape[0])),
#                                  long_name='cube_across',var_name='cube_across',
#                                  dim_coords_and_dims=[(time_coord,0),(model_level_number_coord, 1), (x_dx_coord, 2)],
#                                  aux_coords_and_dims=[(x_across_coord, [0,2]), (y_across_coord, [0,2])])
#     elif z_coord=='geopotential_height':
#        grid_along=iris.cube.Cube(np.zeros((time_coord.shape[0],height_levels.shape[0],x_dx.shape[0])),
#                                  long_name='cube_along',var_name='cube_along',
#                                  dim_coords_and_dims=[(time_coord,0),(geopotential_height_coord, 1), (x_dx_coord, 2)],
#                                  aux_coords_and_dims=[(x_along_coord, [0,2]), (y_along_coord, [0,2])])
    
#        grid_across=iris.cube.Cube(np.ones((time_coord.shape[0],height_levels.shape[0],x_dx.shape[0])),
#                                  long_name='cube_across',var_name='cube_across',
#                                  dim_coords_and_dims=[(time_coord,0),(geopotential_height_coord, 1), (x_dx_coord, 2)],
#                                  aux_coords_and_dims=[(x_across_coord, [0,2]), (y_across_coord, [0,2])])
#     else:
#         raise ValueError('z_coord must be model_level_number or geopotential_height ')
#     return grid_along,grid_across

def extract_cell_cubes_subset(cubelist_in,mask,track,cell,
                              z_coord='model_level_number', height_levels=None):
    
    from iris.analysis import SUM
    from iris import Constraint
    from iris.cube import CubeList
    from iris.coords import AuxCoord
    import numpy as np
    from tobac import mask_cell,mask_cell_surface,get_bounding_box
    from copy import deepcopy
    
    track_i=track[track['cell']==cell]
    
    cubelist_cell_integrated_out=CubeList()
    cubelist_cell_sum=CubeList()

    for time_i in track_i['time'].values:    
        
        logging.debug('start extracting cubes for cell '+str(cell) + ' and time '+str(time_i))

        constraint_time = Constraint(time=time_i)
        mask_i=mask.extract(constraint_time)
        mask_cell_i=mask_cell(mask_i,cell,track_i,masked=False)
        mask_cell_surface_i=mask_cell_surface(mask_i,cell,track_i,masked=False,z_coord=z_coord)

        x_dim=mask_cell_surface_i.coord_dims('projection_x_coordinate')[0]
        y_dim=mask_cell_surface_i.coord_dims('projection_y_coordinate')[0]
        x_coord=mask_cell_surface_i.coord('projection_x_coordinate')
        y_coord=mask_cell_surface_i.coord('projection_y_coordinate')
    
        if (mask_cell_surface_i.core_data()>0).any():
            box_mask_i=get_bounding_box(mask_cell_surface_i.core_data(),buffer=1)
    
            box_mask=[[x_coord.points[box_mask_i[x_dim][0]],x_coord.points[box_mask_i[x_dim][1]]],
                     [y_coord.points[box_mask_i[y_dim][0]],y_coord.points[box_mask_i[y_dim][1]]]]
        else:
            box_mask=[[np.nan,np.nan],[np.nan,np.nan]]
    
        width=20
        dx=500
        x=track_i[track_i['time'].values==time_i]['projection_x_coordinate'].values[0]
        y=track_i[track_i['time'].values==time_i]['projection_y_coordinate'].values[0]

        n_add_width=2
        
        box_slice=[[x-(width+n_add_width)*dx,x+(width+n_add_width)*dx],[y-(width+n_add_width)*dx,y+(width+n_add_width)*dx]]
               
        x_min=np.nanmin([box_mask[0][0],box_slice[0][0]])
        x_max=np.nanmax([box_mask[0][1],box_slice[0][1]])
        y_min=np.nanmin([box_mask[1][0],box_slice[1][0]])
        y_max=np.nanmax([box_mask[1][1],box_slice[1][1]])

        constraint_x=Constraint(projection_x_coordinate=lambda cell: int(x_min) < cell < int(x_max))
        constraint_y=Constraint(projection_y_coordinate=lambda cell: int(y_min) < cell < int(y_max))

        constraint=constraint_time & constraint_x & constraint_y
    
        mask_cell_i=mask_cell_i.extract(constraint_x & constraint_y)
        mask_cell_surface_i=mask_cell_surface_i.extract(constraint_x & constraint_y)
           
        cubelist_i=cubelist_in.extract(constraint)
       
        cubelist_cell_sum.extend(sum_profile_mask(cubelist_i,height_levels,mask_cell_i))
    cubelist_cell_sum_out=cubelist_cell_sum.merge()        
    for cube in cubelist_cell_sum_out:
        cell_time_coord=AuxCoord(track_i['time_cell'].dt.total_seconds().values, units='s',long_name='time_cell')
        cube.add_aux_coord(cell_time_coord,cube.coord_dims('time')[0])

    for cube in cubelist_cell_sum_out:
        cubelist_cell_integrated_out.append(cube.collapsed(('geopotential_height'),SUM))
            
    track_cell_integrated=deepcopy(track_i)
#
    for cube in cubelist_cell_integrated_out:
        track_cell_integrated[cube.name()]=cube.core_data()
        
    return cubelist_cell_sum_out, cubelist_cell_integrated_out,track_cell_integrated

def extract_cell_cubes_subset_2D(cubelist_in,mask,track,cell,z_coord='model_level_number'):
    import iris
    from iris import Constraint
    from iris.cube import CubeList
    import numpy as np
    from tobac import mask_cell_surface,get_bounding_box
    from copy import deepcopy
    
    track_i=track[track['cell']==cell]
    
    cubelist_cell_sum=CubeList()

    for time_i in track_i['time'].values:    
        
        logging.debug('start extracting cubes for cell '+str(cell) + ' and time '+str(time_i))

        constraint_time = Constraint(time=time_i)
        mask_i=mask.extract(constraint_time)
        mask_cell_surface_i=mask_cell_surface(mask_i,cell,track_i,masked=False,z_coord=z_coord)

        x_dim=mask_cell_surface_i.coord_dims('projection_x_coordinate')[0]
        y_dim=mask_cell_surface_i.coord_dims('projection_y_coordinate')[0]
        x_coord=mask_cell_surface_i.coord('projection_x_coordinate')
        y_coord=mask_cell_surface_i.coord('projection_y_coordinate')
    
        if (mask_cell_surface_i.core_data()>0).any():
            box_mask_i=get_bounding_box(mask_cell_surface_i.core_data(),buffer=1)
    
            box_mask=[[x_coord.points[box_mask_i[x_dim][0]],x_coord.points[box_mask_i[x_dim][1]]],
                     [y_coord.points[box_mask_i[y_dim][0]],y_coord.points[box_mask_i[y_dim][1]]]]
        else:
            box_mask=[[np.nan,np.nan],[np.nan,np.nan]]
    
        width=20
        dx=500
        x=track_i[track_i['time'].values==time_i]['projection_x_coordinate'].values[0]
        y=track_i[track_i['time'].values==time_i]['projection_y_coordinate'].values[0]

        n_add_width=2
        
        box_slice=[[x-(width+n_add_width)*dx,x+(width+n_add_width)*dx],[y-(width+n_add_width)*dx,y+(width+n_add_width)*dx]]
               
        x_min=np.nanmin([box_mask[0][0],box_slice[0][0]])
        x_max=np.nanmax([box_mask[0][1],box_slice[0][1]])
        y_min=np.nanmin([box_mask[1][0],box_slice[1][0]])
        y_max=np.nanmax([box_mask[1][1],box_slice[1][1]])

        constraint_x=Constraint(projection_x_coordinate=lambda cell: int(x_min) < cell < int(x_max))
        constraint_y=Constraint(projection_y_coordinate=lambda cell: int(y_min) < cell < int(y_max))

        constraint=constraint_time & constraint_x & constraint_y
    
        mask_cell_surface_i=mask_cell_surface_i.extract(constraint_x & constraint_y)
        cubelist_i=cubelist_in.extract(constraint)

        cubelist_cell_sum.extend(sum_mask_surface(cubelist_i,mask_cell_surface_i))
    logging.debug(str(cubelist_cell_sum))
    cubelist_cell_sum_out=cubelist_cell_sum.merge()
    logging.debug(str(cubelist_cell_sum_out))
    for cube in cubelist_cell_sum_out:
        logging.debug(str(cube))
        logging.debug(str(cube.attributes))
        logging.debug(str(cube.coords()))

    if len(cubelist_cell_sum_out)==6:
        logging.debug(str(iris.util.describe_diff(cubelist_cell_sum_out[0],cubelist_cell_sum_out[3])))
        logging.debug(str(iris.util.describe_diff(cubelist_cell_sum_out[1],cubelist_cell_sum_out[4])))
        logging.debug(str(iris.util.describe_diff(cubelist_cell_sum_out[2],cubelist_cell_sum_out[5])))

    track_cell=deepcopy(track_i)
    for cube in cubelist_cell_sum_out:
        logging.debug(f'cube.shape: {cube.shape}')
        logging.debug(f'len(track_cell): {len(track_cell)}')
        logging.debug(f'cube.coord("time"): {cube.coord("time")}')
        logging.debug(f'track_cell[time]: {track_cell["time"]}')

        track_cell[cube.name()]=cube.core_data()
        
    return cubelist_cell_sum_out,track_cell




def collapse_profile_mask(variable_cube,height_levels_borders,mask_cell,coordinate='geopotential_height',method=None,fillnan=None):
    from iris.coords import DimCoord
    from iris.cube import Cube
    from iris import Constraint
    from numpy import array,transpose
    import numpy as np
    from tobac import mask_cube

    # Set midpoints of height bounds as coordinate values
    height_levels=height_levels_borders[0:-1]+0.5*np.diff(height_levels_borders)
    # Create array of height bounds
    bounds=transpose(array([height_levels_borders[0:-1],height_levels_borders[1:]]))
    #Create DimCoord of height levels and add coordinates
    dim_coordinate=DimCoord(height_levels,standard_name=coordinate,units=variable_cube.coord(coordinate).units,bounds=bounds)
    variable_cube_out=Cube(height_levels,units=variable_cube.units,long_name=variable_cube.name(),dim_coords_and_dims=[(dim_coordinate,0)])
    variable_cube_out.add_aux_coord(variable_cube.coord('time'))    
    for i_height,height in enumerate(height_levels_borders[:-1]):
        constraint_height=Constraint(geopotential_height= lambda cell: height_levels_borders[i_height]<= cell <height_levels_borders[i_height+1])
        mask=mask_cell
        variable_cube_i=mask_cube(variable_cube,mask.core_data()).extract(constraint_height)        
        variable_cube_i=variable_cube_i.collapsed(('model_level_number','x','y'),method)
        variable_cube_out_i=variable_cube_i.core_data()
        if fillnan is not None:
            variable_cube_out_i=fillnan
        variable_cube_out.data[i_height]=variable_cube_out_i
    return variable_cube_out
    
def sum_profile_mask(cubelist_in,height_levels_borders,mask_cell):    
    from iris.cube import CubeList
    from iris.analysis import SUM
    from dask.array.ma import masked_invalid
    cubelist_out = CubeList()
    for variable in cubelist_in:
        #Sum up values for height slices and mask for all cubes in cubelist_in:
        cube_cell_profile=collapse_profile_mask(variable_cube=variable,height_levels_borders=height_levels_borders,mask_cell=mask_cell,
                                             coordinate='geopotential_height',method=SUM)
        cube_cell_profile.rename(variable.name())
        cube_cell_profile.data=masked_invalid(cube_cell_profile.core_data())
        coord_names=[coord.name() for coord in cube_cell_profile.coords()]
        coord_names.remove('time')
        coord_names.remove('geopotential_height')
        for coord in coord_names:
            cube_cell_profile.remove_coord(coord)
        cubelist_out.append(cube_cell_profile)
    return cubelist_out


def sum_mask_surface(cubelist_in,mask_cell_surface):    
    from iris.cube import CubeList
    from iris.analysis import SUM
    from tobac import mask_cube
    cubelist_out = CubeList()
    for variable in cubelist_in:
        #Sum up values for height slices and mask for all cubes in cubelist_in:
        variable_cube_i=mask_cube(variable,mask_cell_surface.core_data())        
        variable_sum=variable_cube_i.collapsed(('x','y'),SUM)
        coord_names=[coord.name() for coord in variable_sum.coords()]
        coord_names.remove('time')
        for coord in coord_names:
            variable_sum.remove_coord(coord)
        cubelist_out.append(variable_sum)
    return cubelist_out

def calculate_alpha(x,y,i):
    """
    Calculate the angle of the track (from x axis) at each timestep based on the point at that 
    timestep and the two adjacent ones. 
    The first and last point are treated special with a calculation based on two points
    
        Parameters
    ----------
    x : numpy ndarray
        x coordinates of the track
    y : numpy ndarray
        y coordinates of the track
    i : integer
        current position on the track

    Returns
    -------
    alpha: float
    angle betweeen track and x axis (range: -pi tp pi)
    """
    # Special case: first point in track:
    if i==0:
        alpha=np.arctan2(y[i+1]-y[i],x[i+1]-x[i])
        
    # Special case: last point in track:
    elif i==len(x)-1:
        alpha=np.arctan2(y[i]-y[i-1],x[i]-x[i-1])
    # General case: all other points in track:
    else:
        alpha_1=np.arctan2(y[i+1]-y[i],x[i+1]-x[i])
        alpha_2=np.arctan2(y[i]-y[i-1],x[i]-x[i-1])
        alpha=np.arctan2(1./2*(np.sin(alpha_1)+np.sin(alpha_2)),1./2*(np.cos(alpha_1)+np.cos(alpha_2)))
    return alpha

def calculate_alpha_all(x,y):
    #calculate angle alpha between cell direction and x-axis for all points in the the track
    from scipy.ndimage.filters import gaussian_filter1d
    alpha=np.nan*np.ones(np.array(x.shape))
    for i in range(len(x)):
        alpha[i]=calculate_alpha(x,y,i)
        
    alpha=gaussian_filter1d(alpha,2.)
    return alpha

def calculate_xy_alphafrom_time(x,y,Time,time_i):
    """
    Calculate the position and the angle of the track (from x axis) at a specific timesteo based on the four adjacent timesteps. 
    The first two and the last two point are treated special with a calculation based on two points
    
        Parameters
    ----------
    x : numpy ndarray
        x coordinates of the track
    y : numpy ndarray
        y coordinates of the track
    Times: numpy ndarray
        list of dattime objexts of the track

    time_i : datetime_object
        current position on the track

    Returns
    -------
    alpha: float
    angle betweeen track and x axis (range: -pi tp pi)
    """
    
    i=1
    # Special case: first point in track:
    if i==0:
        alpha=np.arctan2(y[i+1]-y[i],x[i+1]-x[i])
        
    # Special case: last point in track:
    elif i==len(x)-1:
        alpha=np.arctan2(y[i]-y[i-1],x[i]-x[i-1])
    # General case: all other points in track:
    else:
        alpha_1=np.arctan2(y[i+1]-y[i],x[i+1]-x[i])
        alpha_2=np.arctan2(y[i]-y[i-1],x[i]-x[i-1])
        alpha=np.arctan2(1./2*(np.sin(alpha_1)+np.sin(alpha_2)),1./2*(np.cos(alpha_1)+np.cos(alpha_2)))
        
    return x,y,alpha

def interpolate_alongacross(Processes,Track,cell,dx,width,z_coord='model_level_number',height_levels=np.arange(2,20000,2000)):

    from iris import Constraint
    from iris.cube import CubeList
    Track_cell=Track[Track['cell']==cell]
    time=Track_cell['time'].values
    x=Track_cell['projection_x_coordinate'].values
    y=Track_cell['projection_y_coordinate'].values
    alpha=calculate_alpha_all(x,y)

    cubelist_Processes_along=CubeList()
    cubelist_Processes_across=CubeList()

    for i,time_i in enumerate(time):
        logging.debug(time_i)
        constraint_time=Constraint(time=time_i)
        grid_along,grid_across=make_cubes_alongacross(x=x[i],
                                                      y=y[i],
                                                      alpha=alpha[i],
                                                      cube_sample=Processes.extract(constraint_time)[0],
                                                      dx=dx,width=width,
                                                      z_coord=z_coord,
                                                      height_levels=height_levels)

        Processes_along_i,Processes_across_i=interpolate_to_cubes(Processes.extract(constraint_time),
                                                                  grid_along,grid_across,
                                                                  z_coord=z_coord)

        cubelist_Processes_along.extend(Processes_along_i)
        cubelist_Processes_across.extend(Processes_across_i)
    Processes_along=cubelist_Processes_along.concatenate()
    Processes_across=cubelist_Processes_across.concatenate()

    return Processes_along,Processes_across

def my_interpolate_3D2D_altitude(cube_in,cube_grid,coordinates_in=None,coordinates_out=None,method='linear'):
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator,interp1d

    #Calculate input grid for griddate method
    source_points=[cube_in.coord(c).points for c in coordinates_in]
    source_points_interp=tuple(source_points)

    #Calculate output grid for griddate method
    mesh_coord_0,mesh_coord_1=np.meshgrid(cube_in.coord(coordinates_in[0]).points,cube_grid.coord(coordinates_in[1]).points,indexing='ij')
    mesh_coord_0,mesh_coord_2=np.meshgrid(cube_in.coord(coordinates_in[0]).points,cube_grid.coord(coordinates_in[2]).points,indexing='ij')
    target_points_interp=np.stack(tuple([mesh_coord_0.flatten(),mesh_coord_1.flatten(),mesh_coord_2.flatten()])).T
    #Calcululate interpolate output and bring it into the right shape:

    Interpolator=RegularGridInterpolator(source_points_interp,cube_in.data,method=method,bounds_error=False,fill_value=0)

    target_data=Interpolator(target_points_interp).reshape((cube_in.shape[0],cube_grid.coord(coordinates_out[1]).shape[0]))

    Interpolator_altitude=RegularGridInterpolator(source_points_interp,cube_in.coord(coordinates_out[0]).points,method=method,bounds_error=False,fill_value=0)
    altitude_data=Interpolator_altitude(target_points_interp).reshape((cube_in.shape[0],cube_grid.coord(coordinates_out[1]).shape[0]))

    target_data_altitude=np.nan*np.ones(cube_grid.shape)
    for i in range(target_data_altitude.shape[-1]):
        target_data_altitude_interpolator=interp1d(altitude_data[:,i],target_data[:,i],bounds_error=False,fill_value=0)
        target_data_altitude[...,i]=target_data_altitude_interpolator(cube_grid.coord(coordinates_out[0]).points)

    #Store output in cube and adjust metadata:
    cube_out=cube_grid.copy()
    cube_out.data=target_data_altitude
    cube_out.rename(cube_in.name())
    cube_out.units=cube_in.units
    return(cube_out)

def my_interpolate_3D2D(cube_in,cube_grid,coordinates,method='linear'):
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
   
    #Calculate input grid for griddate method
    source_points=[cube_in.coord(c).points for c in coordinates]
    source_points_interp=tuple(source_points)

    #Calculate output grid for griddate method
    mesh_coord_0,mesh_coord_1=np.meshgrid(cube_grid.coord(coordinates[0]).points,cube_grid.coord(coordinates[1]).points,indexing='ij')
    mesh_coord_0,mesh_coord_2=np.meshgrid(cube_grid.coord(coordinates[0]).points,cube_grid.coord(coordinates[2]).points,indexing='ij')
    target_points_interp=np.stack(tuple([mesh_coord_0.flatten(),mesh_coord_1.flatten(),mesh_coord_2.flatten()])).T

    #Calcululate interpolate output and bring it into the right shape:
    Interpolator=RegularGridInterpolator(source_points_interp,cube_in.data,method=method,bounds_error=False,fill_value=0)
    target_data=Interpolator(target_points_interp).reshape(cube_grid.shape)

    #Store output in cube and adjust metadata:
    cube_out=cube_grid.copy()
    cube_out.data=target_data
    cube_out.rename(cube_in.name())
    cube_out.units=cube_in.units

    return(cube_out)

def interpolate_to_cubes(cubelist_in,grid_along,grid_across,z_coord='model_level_number'):
    from iris.cube import CubeList
    cubelist_out_along= CubeList()
    cubelist_out_across= CubeList()

    for variable in cubelist_in:
        
        
        if z_coord=='model_level_number':
            cube_along=my_interpolate_3D2D(variable,grid_along,coordinates=['model_level_number','projection_y_coordinate','projection_x_coordinate'],method='linear')
            cube_across=my_interpolate_3D2D(variable,grid_across,coordinates=['model_level_number','projection_y_coordinate','projection_x_coordinate'],method='linear')
        elif z_coord=='geopotential_height':        
            if variable.coord('geopotential_height').ndim==1:
                cube_along=my_interpolate_3D2D(variable,grid_along,coordinates=['geopotential_height','projection_y_coordinate','projection_x_coordinate'],method='linear')
                cube_across=my_interpolate_3D2D(variable,grid_across,coordinates=['geopotential_height','projection_y_coordinate','projection_x_coordinate'],method='linear')
            if variable.coord('geopotential_height').ndim>1:
                cube_along=my_interpolate_3D2D_altitude(variable,grid_along,coordinates_in=['model_level_number','projection_y_coordinate','projection_x_coordinate'],coordinates_out=['geopotential_height','x_dx'],method='linear')
                cube_across=my_interpolate_3D2D_altitude(variable,grid_across,coordinates_in=['model_level_number','projection_y_coordinate','projection_x_coordinate'],coordinates_out=['geopotential_height','x_dx'],method='linear')
            

        
        else:
            raise ValueError('z_coord must be model_level_number or geopotential_height ')
        cubelist_out_along.append(cube_along)
        cubelist_out_across.append(cube_across)
        
        del cube_along
        del cube_across
    return cubelist_out_along,cubelist_out_across

def make_cubes_alongacross(x,y,alpha,cube_sample,dx,width,height_levels=None,z_coord='model_level_number',interpolate_time=False):
    """
    Create cubes along and across the cell at the current timestep to use for interpolation to them from other variables 
    Cube contains no data (zeros)
        Parameters
    ----------
    x : float
        x coordinate of the track
    y : float
        y coordinate of the track
    alpha: float
        angle betweeen track and x axis (range: -pi tp pi)
    cube_sample: iris.cube.Cube
                  sample cube of original data (for time and vertical coordinate)
    dx     : float
             horizontal spacing interpolated to (in m)
    width  : int 
             half width of the interpolated cube in (numbers of points)
    height_levels: numpy ndarray
                   vertical levels (geopotential_height) to interpolate to (in m)
            
                 
    Returns
    -------
    
    grid_along   iris.cube.Cube
    
    grid_across  iris.cube.Cube

    """
    time_coord=cube_sample.coord('time')
    if z_coord=='model_level_number':
        model_level_number_coord=iris.coords.DimCoord(cube_sample.coord('model_level_number')[::5].points,
                                                      long_name='model_level_number',
                                                      var_name='model_level_number',units=' ')
    elif z_coord=='geopotential_height':
        geopotential_height_coord=iris.coords.DimCoord(height_levels,long_name='geopotential_height',var_name='geopotential_height',units='m')
    else:
        raise ValueError('z_coord must be model_level_number or geopotential_height ')
       
      

    # create new x coordinate (new cube)
    x_dx=dx*(np.arange(-1.*width,1.*width+1))
    x_dx_coord=iris.coords.DimCoord(x_dx,long_name='x_dx',var_name='x_dx',units='m')
        
    #create coordinates x and y for new cube:  
    x_along=x+np.expand_dims(x_dx,0)*np.cos(alpha)
    y_along=y+np.expand_dims(x_dx,0)*np.sin(alpha)  
    
    x_along_coord=iris.coords.AuxCoord(x_along,standard_name='projection_x_coordinate',units='m')
    y_along_coord=iris.coords.AuxCoord(y_along,standard_name='projection_y_coordinate',units='m')
    
    x_across=x+np.expand_dims(x_dx,0)*np.sin(alpha)
    y_across=y-np.expand_dims(x_dx,0)*np.cos(alpha)    
    x_across_coord=iris.coords.AuxCoord(x_across,standard_name='projection_x_coordinate',long_name='x',var_name='x',units='m')
    y_across_coord=iris.coords.AuxCoord(y_across,standard_name='projection_y_coordinate',long_name='y',var_name='y',units='m')
    
    if z_coord=='model_level_number':
       
       grid_along=iris.cube.Cube(np.zeros((time_coord.shape[0],model_level_number_coord.shape[0],x_dx_coord.shape[0])),
                                 long_name='cube_along',var_name='cube_along',
                                 dim_coords_and_dims=[(time_coord,0),(model_level_number_coord, 1), (x_dx_coord, 2)],
                                 aux_coords_and_dims=[(x_along_coord, [0,2]), (y_along_coord, [0,2])])
    
       grid_across=iris.cube.Cube(np.zeros((time_coord.shape[0],model_level_number_coord.shape[0],x_dx_coord.shape[0])),
                                 long_name='cube_across',var_name='cube_across',
                                 dim_coords_and_dims=[(time_coord,0),(model_level_number_coord, 1), (x_dx_coord, 2)],
                                 aux_coords_and_dims=[(x_across_coord, [0,2]), (y_across_coord, [0,2])])
    elif z_coord=='geopotential_height':
    
       grid_along=iris.cube.Cube(np.zeros((time_coord.shape[0],height_levels.shape[0],x_dx.shape[0])),
                                 long_name='cube_along',var_name='cube_along',
                                 dim_coords_and_dims=[(time_coord,0),(geopotential_height_coord, 1), (x_dx_coord, 2)],
                                 aux_coords_and_dims=[(x_along_coord, [0,2]), (y_along_coord, [0,2])])
    
       grid_across=iris.cube.Cube(np.ones((time_coord.shape[0],height_levels.shape[0],x_dx.shape[0])),
                                 long_name='cube_across',var_name='cube_across',
                                 dim_coords_and_dims=[(time_coord,0),(geopotential_height_coord, 1), (x_dx_coord, 2)],
                                 aux_coords_and_dims=[(x_across_coord, [0,2]), (y_across_coord, [0,2])])
    else:
        raise ValueError('z_coord must be model_level_number or geopotential_height ')
    return grid_along,grid_across


def interpolate_alongacross_mean(Processes,Track,cell,dx,width,z_coord='model_level_number',height_level_borders=np.arange(0,20000,2000)):

    from iris import Constraint
    from iris.cube import CubeList
    Track_cell=Track[Track['cell']==cell]
    time=Track_cell['time'].values
    x=Track_cell['projection_x_coordinate'].values
    y=Track_cell['projection_y_coordinate'].values
    alpha=calculate_alpha_all(x,y)

    cubelist_Processes_along=CubeList()
    cubelist_Processes_across=CubeList()

               



    for i,time_i in enumerate(time):
        logging.debug(time_i)
        
        constraint_time=Constraint(time=time_i)
        
        n_add_width=2
        box_slice=[[x[i]-(width+n_add_width)*dx,x[i]+(width+n_add_width)*dx],[y[i]-(width+n_add_width)*dx,y[i]+(width+n_add_width)*dx]]

        x_min=box_slice[0][0]
        x_max=box_slice[0][1]
        y_min=box_slice[1][0]
        y_max=box_slice[1][1]

        constraint_x=Constraint(projection_x_coordinate=lambda cell: int(x_min) < cell < int(x_max))
        constraint_y=Constraint(projection_y_coordinate=lambda cell: int(y_min) < cell < int(y_max))

        constraint=constraint_time & constraint_x & constraint_y

        grid_along,grid_across=make_cubes_alongacross_mean(x=x[i],
                                                      y=y[i],
                                                      alpha=alpha[i],
                                                      cube_sample=Processes.extract(constraint)[0],
                                                      dx=dx,width=width,
                                                      height_level_borders=height_level_borders)

        Processes_along_i,Processes_across_i=interpolate_to_cubes_mean(Processes.extract(constraint),
                                                                  grid_along,grid_across,
                                                                  height_level_borders=height_level_borders)

        cubelist_Processes_along.extend(Processes_along_i)
        cubelist_Processes_across.extend(Processes_across_i)
    Processes_along=cubelist_Processes_along.concatenate()
    Processes_across=cubelist_Processes_across.concatenate()

    return Processes_along,Processes_across

def interpolate_to_cubes_mean(cubelist_in,grid_along,grid_across,height_level_borders=None):
    from iris.cube import CubeList
    cubelist_out_along= CubeList()
    cubelist_out_across= CubeList()

    for variable in cubelist_in:

        if variable.coord('geopotential_height').ndim==1:
            cube_along=my_interpolate_3D2D_mean(variable,grid_along,coordinates=['geopotential_height','projection_y_coordinate','projection_x_coordinate'],method='linear',height_level_borders=height_level_borders)
            cube_across=my_interpolate_3D2D_mean(variable,grid_across,coordinates=['geopotential_height','projection_y_coordinate','projection_x_coordinate'],method='linear',height_level_borders=height_level_borders)
        # if variable.coord('geopotential_height').ndim>1:
        #     cube_along=my_interpolate_3D2D_altitude(variable,grid_along,coordinates_in=['model_level_number','projection_y_coordinate','projection_x_coordinate'],coordinates_out=['geopotential_height','x_dx'],method='linear')
        #     cube_across=my_interpolate_3D2D_altitude(variable,grid_across,coordinates_in=['model_level_number','projection_y_coordinate','projection_x_coordinate'],coordinates_out=['geopotential_height','x_dx'],method='linear')

        cubelist_out_along.append(cube_along)
        cubelist_out_across.append(cube_across)
        
        del cube_along
        del cube_across
    return cubelist_out_along,cubelist_out_across



def make_cubes_alongacross_mean(x,y,alpha,cube_sample,dx,width,height_level_borders=None,interpolate_time=False):
    """
    Create cubes along and across the cell at the current timestep to use for interpolation to them from other variables 
    Cube contains no data (zeros)
        Parameters
    ----------
    x : float
        x coordinate of the track
    y : float
        y coordinate of the track
    alpha: float
        angle betweeen track and x axis (range: -pi tp pi)
    cube_sample: iris.cube.Cube
                  sample cube of original data (for time and vertical coordinate)
    dx     : float
             horizontal spacing interpolated to (in m)
    width  : int 
             half width of the interpolated cube in (numbers of points)
    height_level_borders: numpy ndarray
                   vertical level borders (geopotential_height) to average between (in m)
            
                 
    Returns
    -------
    
    grid_along   iris.cube.Cube
    
    grid_across  iris.cube.Cube

    """
    time_coord=cube_sample.coord('time')

    height_levels=height_level_borders[0:-1]+0.5*np.diff(height_level_borders)
    geopotential_height_coord=iris.coords.DimCoord(height_levels,long_name='geopotential_height',var_name='geopotential_height',units='m')
       
      

    # create new x coordinate (new cube)
    x_dx=dx*(np.arange(-1.*width,1.*width+1))
    x_dx_coord=iris.coords.DimCoord(x_dx,long_name='x_dx',var_name='x_dx',units='m')
        
    #create coordinates x and y for new cube:  
    x_along=x+np.expand_dims(x_dx,0)*np.cos(alpha)
    y_along=y+np.expand_dims(x_dx,0)*np.sin(alpha)  
    
    x_along_coord=iris.coords.AuxCoord(x_along,standard_name='projection_x_coordinate',units='m')
    y_along_coord=iris.coords.AuxCoord(y_along,standard_name='projection_y_coordinate',units='m')
    
    x_across=x+np.expand_dims(x_dx,0)*np.sin(alpha)
    y_across=y-np.expand_dims(x_dx,0)*np.cos(alpha)    
    x_across_coord=iris.coords.AuxCoord(x_across,standard_name='projection_x_coordinate',long_name='x',var_name='x',units='m')
    y_across_coord=iris.coords.AuxCoord(y_across,standard_name='projection_y_coordinate',long_name='y',var_name='y',units='m')

    grid_along=iris.cube.Cube(np.zeros((time_coord.shape[0],height_levels.shape[0],x_dx.shape[0])),
                              long_name='cube_along',var_name='cube_along',
                              dim_coords_and_dims=[(time_coord,0),(geopotential_height_coord, 1), (x_dx_coord, 2)],
                              aux_coords_and_dims=[(x_along_coord, [0,2]), (y_along_coord, [0,2])])

    grid_across=iris.cube.Cube(np.ones((time_coord.shape[0],height_levels.shape[0],x_dx.shape[0])),
                              long_name='cube_across',var_name='cube_across',
                              dim_coords_and_dims=[(time_coord,0),(geopotential_height_coord, 1), (x_dx_coord, 2)],
                              aux_coords_and_dims=[(x_across_coord, [0,2]), (y_across_coord, [0,2])])
    
    return grid_along,grid_across

def my_interpolate_3D2D_mean(cube_in,cube_grid,coordinates,method='linear',height_level_borders=None):
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator

    # Create array of height bounds
    bounds=np.transpose(np.array([height_level_borders[0:-1],height_level_borders[1:]]))

    #Calculate input grid for griddate method
    source_points=[cube_in.coord(c).points for c in coordinates]
    source_points_interp=tuple(source_points)
    coord_0=cube_in.coord(coordinates[0]).points
    coord_1=cube_grid.coord(coordinates[1]).points
    coord_2=cube_grid.coord(coordinates[2]).points
    #Calculate output grid for griddate method
    mesh_coord_0,mesh_coord_1=np.meshgrid(coord_0,coord_1,indexing='ij')
    mesh_coord_0,mesh_coord_2=np.meshgrid(coord_0,coord_2,indexing='ij')
    target_points_interp=np.stack(tuple([mesh_coord_0.flatten(),mesh_coord_1.flatten(),mesh_coord_2.flatten()])).T
    
    #Calcululate interpolate output and bring it into the right shape:
    Interpolator=RegularGridInterpolator(source_points_interp,cube_in.data,method=method,bounds_error=False,fill_value=0)
    target_data_intermediate=Interpolator(target_points_interp).reshape((cube_in.shape[0],cube_grid.shape[2]))
    
    target_data_final=np.zeros(cube_grid.shape)
    for i in range(bounds.shape[0]):
        target_data_final[0,i]=np.mean(target_data_intermediate[(coord_0>=bounds[i,0]) & (coord_0<bounds[i,1])],axis=0)

    #Store output in cube and adjust metadata:
    cube_out=cube_grid.copy()
    cube_out.data=target_data_final
    cube_out.rename(cube_in.name())
    cube_out.units=cube_in.units

    return cube_out


