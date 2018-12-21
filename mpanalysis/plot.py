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
        x = cube_1.coord('time_cell').points/60
    else:
        x = cube_1.coord('time').points
     
    if x_shift:
        x=x+x_shift

    #x = np.array([cube_1.coord('time').units.num2date(cube_1.coord('time').points[0]).hour * 60 + cube_1.coord('time').units.num2date(cube_1.coord('time').points[i]).minute for i in range(len(cubelist_in[0].coord('time').units.num2date(cube_1.coord('time').points)))])
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
        plot_quiver = axes.quiver(x, y, u, w, w, color='grey', cmap=plt.cm.coolwarm, clim=[-5.0, 5.0], scale=1000, linewidth=0.7, zorder=3, label='Wind vectors',rasterized=False)
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
        names_in['Other']='Other'
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
