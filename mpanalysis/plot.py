import numpy as np
import matplotlib.pyplot as plt

def plot_piecharts_color_time(cubelist_in, Aux,
                              colors_in=None,names_in=None,
                              axes=None,
                              scaling='linear', minvalue=None, maxvalue=None,piecharts_rasterized=True,
                              legend_piecharts=False,legend_overlay=False, loc_legend='upper left',
                              legend_piecharts_pos=(1,1),legend_overlay_pos=(1,0.5),
                              fontsize_legend=plt.rcParams['legend.fontsize'],
                              xlabel=True, ylabel=True,xlim=None,ylim=None,
                              overlay=False, contour_ice=None, contour_liquid=None,
                              scale=True,vscale=None, loc_scale='upper left', unit_scale='',
                              fontsize_scale=None,x_shift=None):
    
    from matplotlib import patches
#    from scipy.ndimage import zoom
    from matplotlib import lines
    from piecharts import piecharts

    cube_1 = cubelist_in[0]
    ratios = np.zeros([len(cubelist_in), cube_1.coord('time').shape[0], cube_1.coord('geopotential_height').shape[0]])
    colors = []
    for i, cube in enumerate(cubelist_in):
        ratios[i, :, :] = np.abs(cube.data)/np.diff(cube_1.coord('geopotential_height').points)[0]
        if cube.name() in colors_in:
            colors.append(colors_in[cube.name()])
        else:
            colors.append('grey')
            
    ratios[ratios>1e30]=0

    coord_names=[coord.name() for coord in cube_1.coords()]
    if 'time_cell' in coord_names:
        time_coord=cube_1.coord('time_cell')
    elif 'time' in coord_names:
        time_coord=cube_1.coord('time')
    else:
        raise ValueError('missing time coordinate, must have coordinate time or time_cell')

    if time_coord.units=='minute':
        x = time_coord.points
    elif time_coord.units=='s':
        x = time_coord.points/60
    else:
        x = [dt.total_seconds()/60 for dt in time_coord.units.num2date(time_coord.points)-time_coord.units.num2date(time_coord.points[0])]

    if x_shift:
        x=x+x_shift

    y = cube_1.coord('geopotential_height').points / 1000
    if xlim==None:
         xlim=[x[0] - np.diff(x)[0], x[-1] + np.diff(x)[0]]
    if ylim==None: 
         xlim=[y[0] - np.diff(y)[0], y[-1] + np.diff(y)[0]]
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    
    
    if overlay:
        handles_overlay = []
        #n_zoom = 5
        # LWC = Aux.extract_strict('liquid water content').data
        # if np.unique(LWC).size > 1:
        #     axes.contour(zoom(x[1:-1], n_zoom), zoom(y[1:-1], n_zoom), zoom(LWC[1:-1, 1:-1], n_zoom),
        #                  levels=[contour_liquid],
        #                  linewidth=0.5, linestyle='-', colors='darkblue', zorder=2, label='Liquid water content')
        # handles_overlay.append(lines.Line2D([], [], color='darkblue', label='Liquid water cont.'))
        # IWC = Aux.extract_strict('ice water content').data
        # if np.unique(IWC).size > 1:
        #     axes.contour(zoom(x[1:-1], n_zoom), zoom(y[1:-1], n_zoom), zoom(IWC[1:-1, 1:-1], n_zoom), levels=[contour_ice],
        #                  linewidth=0.5, linestyle='-', colors='darkgrey', zorder=2, label='Ice water content')
        # handles_overlay.append(lines.Line2D([], [], color='darkgrey', label='Ice water cont.'))
        
        xx,yy=np.meshgrid(x,y,indexing='xy')

        T = Aux.extract_strict('air temperature').data
        if np.unique(T).size > 1:
            axes.contour(x, y, T.transpose(),
                         levels=[273.15],
                         linewidth=0.5, linestyle='-', colors='darkred', label='Melting level')

            # axes.contour(zoom(x[1:-1], n_zoom), zoom(y[1:-1], n_zoom), zoom(yy[1:-1, 1:-1], n_zoom),
            #              levels=[1,2,4,6,10,13,17],
            #              linewidth=1, linestyle='-', color='darkred', label='Melting level fake')
  
        handles_overlay.append(lines.Line2D([], [], color='darkred', label='Melting level'))
        
        if legend_overlay:
            legend_contours = axes.legend(handles=handles_overlay, loc=loc_legend, bbox_to_anchor=legend_overlay_pos, frameon=False,fontsize=fontsize_legend)
            axes.add_artist(legend_contours)


    piecharts(ratios, x, y, colors, axes=axes, scaling=scaling, vmin=0, vmax=maxvalue,
              scale=True,vscale=vscale, loc_scale='upper left', unit_scale=unit_scale, fontsize_scale=fontsize_scale,
              legend=False, loc_legend='upper right', labels=None, rasterized=piecharts_rasterized, zorder=4)
    if xlabel:
        xlabel_str = 'Simulation time [min]'
        axes.set_xlabel(xlabel_str)
    if ylabel:
        ylabel_str = 'altitude [km]'
        axes.set_ylabel(ylabel_str)
    if legend_piecharts:
        patches_legend = []
        for cube in cubelist_in:
            patches_legend.append(patches.Patch(color=colors_in[cube.name()], label=names_in[cube.name()]))

        legend_pie = axes.legend(handles=patches_legend, loc='upper left',
                                 bbox_to_anchor=legend_piecharts_pos, frameon=False,
                                 fontsize=fontsize_legend)
        axes.add_artist(legend_pie)
    return axes

def plot_piecharts_slice(cubelist_in, Aux, colors_in=None,names_in=None, axes=None,
                         scaling='linear', fraction=1, minvalue=0, maxvalue=None,piecharts_rasterized=True,
                         xlabel=True, ylabel=True,
                         scale=True,vscale=None, loc_scale='upper left', unit_scale='', fontsize_scale=plt.rcParams['legend.fontsize'],
                         legend_piecharts=False,legend_overlay=False, loc_legend='upper left',legend_piecharts_pos=(1,0.3),legend_overlay_pos=(1,0.5), fontsize_legend=plt.rcParams['legend.fontsize'],
                         labels=None,
                         xlim=None,ylim=None,z_coord=None,
                         overlay=False, contour_ice=None, contour_liquid=None,
                         wind_vector=False,wind_vector_pos=(0.75, 0.95)):
    from scipy.ndimage import zoom
    from matplotlib import patches
    from matplotlib import lines
    from matplotlib import font_manager
    from piecharts import piecharts 

    cube_1 = cubelist_in[0]
    if z_coord == None:
        coord_names=[coord.name() for coord in cube_1.coords()]
        if 'geopotential_height' in coord_names:
            z_coord='geopotential_height'
        elif 'model_level_number' in coord_names:
            z_coord='model_level_number'

    
    ratios = np.zeros([len(cubelist_in), cube_1.coord('x_dx').shape[0], cube_1.coord(z_coord).shape[0]])
    colors = []
        
    if names_in is None:
        names_in={}
        for name in colors_in:
            names_in[name]=name
    for i, cube in enumerate(cubelist_in):
        ratios[i, :, :] = np.abs(cube.data.transpose())
        if cube.name() in colors_in:
            colors.append(colors_in[cube.name()])
        else:
            colors.append('grey')

    x = cube_1.coord('x_dx').points / 1000
    
    if z_coord == 'geopotential_height':
        y = cube_1.coord('geopotential_height').points / 1000
    elif z_coord == 'model_level_number':
        y = cube_1.coord('model_level_number').points
    if not xlim:
        xlim=[x[0] - np.diff(x)[0], x[-1] + np.diff(x)[0]]
    axes.set_xlim(xlim)
    if not ylim:
        ylim=[y[0],y[-1]]
    axes.set_ylim(ylim)
    if fontsize_scale == None:
        fontsize_main = plt.rcParams['font.size']
        fontsize_scale = font_manager.font_scalings[fontsize_legend] * fontsize_main
    if overlay:
        u = Aux.extract_strict('u_x').data
        w = Aux.extract_strict('upward_air_velocity').data
        if wind_scale is 'auto':
            scale=wind_scale
        else:
            dx=(diff(x)[0]+diff(y)[0])/2
            scale=dx/wind_scale
            
        plot_quiver = axes.quiver(x, y, u, w, w, color='grey', cmap=plt.cm.coolwarm, clim=[-5.0, 5.0], scale=scale,scale_units='xy', linewidth=0.7, zorder=3, label='Wind vectors',rasterized=False)
        if wind_vector:
            axes.quiverkey(plot_quiver, wind_vector_pos[0], wind_vector_pos[1], 30, '30 m/s', color='k', labelpos='E', coordinates='axes', fontproperties={'size': fontsize_scale})
        n_zoom = 5
        handles_overlay = []
        handles_overlay.append(plt.quiverkey(plot_quiver, 100, 100, 2, 'arrow 1', coordinates='data'))
        LWC = Aux.extract_strict('liquid water content').data
        if np.unique(LWC).size > 1:
            axes.contour(zoom(x[1:-1], n_zoom), zoom(y[1:-1], n_zoom), zoom(LWC[1:-1, 1:-1], n_zoom),
                         levels=[contour_liquid],
                         linewidth=0.5, linestyle='-', colors='darkblue', zorder=2, label='Liquid water content')
        handles_overlay.append(lines.Line2D([], [], color='darkblue', label='Liquid water cont.'))
        IWC = Aux.extract_strict('ice water content').data
        if np.unique(IWC).size > 1:
            axes.contour(zoom(x[1:-1], n_zoom), zoom(y[1:-1], n_zoom), zoom(IWC[1:-1, 1:-1], n_zoom), levels=[contour_ice],
                         linewidth=0.5, linestyle='-', colors='darkgrey', zorder=2, label='Ice water content')
        handles_overlay.append(lines.Line2D([], [], color='darkgrey', label='Ice water cont.'))
        T = Aux.extract_strict('air temperature').data
        if np.unique(T).size > 1:
            axes.contour(zoom(x[1:-1], n_zoom), zoom(y[1:-1], n_zoom), zoom(T[1:-1, 1:-1], n_zoom),
                         levels=[273.15], linewidth=0.5, linestyle=':', colors='darkred', zorder=2, label='Melting level')
        handles_overlay.append(lines.Line2D([], [], color='darkred', label='Melting level'))
        if legend_overlay:
            legend_contours = axes.legend(handles=handles_overlay, loc=loc_legend, bbox_to_anchor=legend_overlay_pos, frameon=False,fontsize=fontsize_legend)
            axes.add_artist(legend_contours)
    piecharts(ratios, x, y, colors, axes=axes, scaling=scaling, vmin=0, vmax=maxvalue,
              scale=True,vscale=vscale, loc_scale='upper left', unit_scale=unit_scale, fontsize_scale=fontsize_scale,
              legend=False, loc_legend='upper right', labels=None, rasterized=piecharts_rasterized, zorder=4)
    
    if legend_piecharts:
        patches_legend = []
        # names_in['Other']='Other'
        for cube in cubelist_in:
            patches_legend.append(patches.Patch(color=colors_in[cube.name()], label=names_in[cube.name()]))

        legend_pie = axes.legend(handles=patches_legend,
                                       loc=loc_legend,
                                       bbox_to_anchor=legend_piecharts_pos, frameon=False,fontsize=fontsize_legend)
        axes.add_artist(legend_pie)
    if z_coord == 'geopotential_height':
        ylabel_str = 'altitude [km]'
    elif z_coord == 'model_level_number':
            ylabel_str = 'model level'

    xlabel_str = 'x [km]'
    if xlabel:
        axes.set_xlabel(xlabel_str)
    if ylabel:
        axes.set_ylabel(ylabel_str)
    #axes.set_xticks([-10, -10, 0, 10, 20])
    return axes

def plot_processes_color_time(Processes,Aux,
                              microphysics_scheme=None,
                              colors_processes='all',
                              axes=None, **kargs):
    from mpdiag import processes_colors
    processes_colors, processes_names = processes_colors(microphysics_scheme=microphysics_scheme,
                                                         colors_processes=colors_processes)
    axes = plot_piecharts_color_time(Processes,Aux, colors_in=processes_colors,names_in=processes_names, axes=axes, **kargs)
    return axes


def plot_hydrometeors_color_time(Hydrometeors,Aux,
                                 axes=None,
                                 microphysics_scheme=None,
                                 **kargs):
    from mpdiag import hydrometeors_colors
    dict_hydrometeors_colors, dict_hydrometeors_names = hydrometeors_colors(microphysics_scheme=microphysics_scheme)
    axes = plot_piecharts_color_time(Hydrometeors,Aux,
                                     colors_in=dict_hydrometeors_colors,
                                     names_in=dict_hydrometeors_names, 
                                     axes=axes, **kargs)
    return axes

def plot_processes_color_slice(Processes, Aux,
                               microphysics_scheme=None, colors_processes='all',
                               axes=None, **kwargs):
    from mpdiag import processes_colors
    dict_processes_colors, dict_processes_names = processes_colors(microphysics_scheme=microphysics_scheme, colors_processes=colors_processes)
    axes = plot_piecharts_slice(Processes, Aux,
                                axes=axes,
                                colors_in=dict_processes_colors,names_in=dict_processes_names,
                                **kwargs)
    return axes


def plot_hydrometeors_color_slice(Hydrometeors, Aux,
                                  microphysics_scheme=None, colors_processes='all',
                                  axes=None,
                                  **kwargs):
    from mpdiag import hydrometeors_colors
    dict_hydrometeors_colors, dict_hydrometeors_names = hydrometeors_colors(microphysics_scheme=microphysics_scheme)
    axes = plot_piecharts_slice(Hydrometeors, Aux,
                                axes=axes,
                                colors_in=dict_hydrometeors_colors,names_in=dict_hydrometeors_names,
                                **kwargs)
    return axes


# Plots for Composites:

def plot_processes_composite_integrated(Processes_Composite,
                                        maxvalue=None,mp=None,aggregate_min=None,
                                        xlim=None,ylim=None,xlim_profile=None,ylim_integrated=None,title=None,
                                        figsize=(20/2.54,10/2.54),height_ratios=[1.8, 1],width_ratios=[4, 1]):
    from matplotlib.axes import Axes
    from matplotlib.gridspec import GridSpec
    from mpdiag import processes_colors
    from iris.analysis import MEAN,SUM
    from copy import deepcopy
    from iris.coord_categorisation import add_categorised_coord
    from iris.cube import CubeList
    #from mpanalysis.plot import plot_processes_color_time

    dict_processes_colors, dict_processes_names = processes_colors(microphysics_scheme=mp,colors_processes='lumped')

    if xlim is None:
        xlim=[Processes_Composite[0].coord('time').points[0],
              Processes_Composite[0].coord('time').points[-1]]
    if ylim is None:
        ylim=[Processes_Composite[0].coord('geopotential_height').points[0]/1000,
              Processes_Composite[0].coord('geopotential_height').points[-1]/1000]

    fig, ax = plt.subplots(nrows=2, ncols=2,
                           #sharex='col', sharey='row',
                            gridspec_kw={'height_ratios': height_ratios, 'width_ratios': width_ratios},
                            figsize=figsize)
    
    fig.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.93,wspace=0.1,hspace=0.2)

    Processes_Composite_copy=deepcopy(Processes_Composite)
    
    if aggregate_min is not None:

        def get_min(coord,value):
                minutes = value
                return np.floor(minutes/aggregate_min)*aggregate_min+aggregate_min/2

        processes_aggregated=CubeList()
        for cube in Processes_Composite_copy:
            if aggregate_min==5:
                add_categorised_coord(cube, 'time_aggregated', 'time', get_min)
                processes_aggregated.append(cube.aggregated_by(['time_aggregated'], MEAN))

        processes_piecharts=processes_aggregated
    else:
        processes_piecharts=Processes_Composite

    plot_processes_color_time(processes_piecharts,
                                      None,
                                      axes=ax[0,0], 
                                      microphysics_scheme=mp, colors_processes='lumped', 
                                      scaling='linear', minvalue=0, maxvalue=maxvalue,vscale=maxvalue,
                                      piecharts_rasterized=False,
                                      legend_piecharts=True,fontsize_legend=6,legend_piecharts_pos=(1.05,-0.25),
                                      legend_overlay=False,legend_overlay_pos=(1,-0.5),
                                      overlay=False,
                                      xlabel=False, ylabel=False,
                                      xlim=xlim,
                                      ylim=ylim,
                                      scale=True, unit_scale='kg s$^{-1}$ m$^{-1}$', fontsize_scale=6,
                                      x_shift=0)

    ax[0,0].plot([0,0],[0,20000],color='grey',ls='-')
    for cube in Processes_Composite:
        color=dict_processes_colors[cube.name()]
        ax[0,1].plot(abs(cube.collapsed(('time'),SUM).data),cube.coord('geopotential_height').points/1000,color=color)
        ax[1,0].plot(cube.coord('time').points,abs(cube.collapsed(('geopotential_height'),SUM).data),color=color)
        #ax1[1,0].set_ylim(0,1000)
#     ax[1,0].plot([0,0],[0,2e10],color='grey',ls='-')
    ax[1,1].axis('off')
    ax[1,0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax[0,1].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax[1,0].set_xlim(xlim)    
    ax[1,0].set_ylim(ylim_integrated)
    ax[0,1].set_ylim(ylim)
    ax[0,1].set_xlim(xlim_profile)
    ax[0,0].set_ylabel('altitude (km)')
    ax[1,0].set_xlabel('time (min)')
    ax[1,0].set_ylabel('integrated (kg s$^{-1}$)')
    ax[0,0].xaxis.set_tick_params(labelbottom=False)
    ax[0,1].yaxis.set_tick_params(labelleft=False)
    ax[0,1].set_xlabel('integrated (kg m$^{-1}$)',labelpad=10)
    if title:  
        ax[0,0].set_title(title,loc='left')
    
    return fig

def plot_hydrometeors_composite_integrated(Hydrometeors_Composite,
                                           maxvalue=None,aggregate_min=None,mp=None, 
                                           xlim=None,ylim=None,xlim_profile=None,ylim_integrated=None,
                                           title=None,
                                           figsize=(20/2.54,10/2.54),height_ratios=[1.8, 1],width_ratios=[4, 1]):

    from mpdiag import hydrometeors_colors
    from iris.analysis import MEAN,SUM
    from iris.coord_categorisation import add_categorised_coord
    from iris.cube import CubeList

    from copy import deepcopy
    dict_hydrometeors_colors, dict_hydrometeors_names = hydrometeors_colors(microphysics_scheme=mp)

    if xlim is None:
        xlim=[Hydrometeors_Composite[0].coord('time').points[0],Hydrometeors_Composite[0].coord('time').points[-1]]
    if ylim is None:
        ylim=[Hydrometeors_Composite[0].coord('geopotential_height').points[0]/1000,Hydrometeors_Composite[0].coord('geopotential_height').points[-1]/1000]

    fig, ax = plt.subplots(nrows=2, ncols=2, 
                           #sharex='col', sharey='row',
                               gridspec_kw={'height_ratios': height_ratios, 'width_ratios': width_ratios},
                               figsize=figsize)
    
    fig.subplots_adjust(left=0.1,right=0.95,bottom=0.1,top=0.93,wspace=0.1,hspace=0.2)

    
    Hydrometeors_Composite_copy=deepcopy(Hydrometeors_Composite)
    
    if aggregate_min is not None:

        def get_min(coord,value):
                minutes = value
                return np.floor(minutes/aggregate_min)*aggregate_min+aggregate_min/2

        hydrometeors_aggregated=CubeList()
        for cube in Hydrometeors_Composite_copy:
            if aggregate_min==5:
                add_categorised_coord(cube, 'time_aggregated', 'time', get_min)
                hydrometeors_aggregated.append(cube.aggregated_by(['time_aggregated'], MEAN))

        hydrometeors_piecharts=hydrometeors_aggregated
    else:
        hydrometeors_piecharts=hydrometeors_Composite

    
    plot_hydrometeors_color_time(hydrometeors_piecharts,
                                  Aux=None,
                                  axes=ax[0,0], 
                                  microphysics_scheme=mp, 
                                  scaling='linear', minvalue=0, maxvalue=maxvalue,vscale=maxvalue,
                                  piecharts_rasterized=False,
                                  legend_piecharts=True,fontsize_legend=6,legend_piecharts_pos=(1.05,-0.25),
                                  legend_overlay=False,legend_overlay_pos=(1,-0.6),
                                  overlay=False,
                                  xlabel=False, ylabel=False,
                                  xlim=xlim,
                                  ylim=ylim,
                                  scale=True, unit_scale='kg m$^{-1}$ s$^{-1}$', fontsize_scale=6,
                                  x_shift=0)

    ax[0,0].plot([0,0],[0,20000],color='grey',ls='-')
    for cube in Hydrometeors_Composite:
        color=dict_hydrometeors_colors[cube.name()]
        ax[0,1].plot(cube.collapsed(('time'),MEAN).data,cube.coord('geopotential_height').points/1000,color=color)
        ax[1,0].plot(cube.coord('time').points,cube.collapsed(('geopotential_height'),SUM).data,color=color)
        #ax1[1,0].set_ylim(0,1000)
#     ax[1,0].plot([0,0],[0,2e10],color='grey',ls='-')
    ax[1,1].axis('off')
    
    ax[0,0].set_ylabel('altitude (km)')
    ax[1,0].set_xlabel('time (min)')
    ax[1,0].set_xlim(xlim)    
    ax[1,0].set_ylim(ylim_integrated)
    ax[0,1].set_ylim(ylim)
    ax[0,1].set_xlim(xlim_profile)

    ax[1,0].set_ylabel('integrated (kg $^{-1}$)')
    ax[1,0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax[0,1].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax[0,0].xaxis.set_tick_params(labelbottom=False)
    ax[0,1].yaxis.set_tick_params(labelleft=False)
    ax[0,1].set_xlabel('integrated (kg m$^{-1}$)',labelpad=10)
    
    if title:  
        ax[0,0].set_title(title,loc='left')
    
    return fig