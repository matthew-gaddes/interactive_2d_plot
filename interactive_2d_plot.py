#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 17:15:23 2020

@author: matthew
"""
import sys
sys.path.append('/home/matthew/university_work/python_stuff/python_scripts/')
from small_plot_functions import matrix_show



#%% function attempt 


def exploratory_2d_plot(xy, c, spatial_data = None, temporal_data = None,
                        inset_axes_side = 0.1, arrow_length = 0.1, figsize = (10,6), labels = None):
    """ Data are plotted in a 2D space, and when hovering over a point, further information about it (e.g. what image it is)  appears in an inset axes.  
    Inputs:
        xy | rank 2 array | e.g. 2x100, the x and y positions of each data
        c | rank 1 array | e.g. 100, value used to set the colour of each data point
        spatial_data | dict or None | contains 'images_r2' in which the images are stored as row vectors and 'mask', which converts a row vector back to a masked array
        temporal_data | dict or None | contains 'timecourses' as time signals as row vectors and 'xvals' which are the times for each item in the timecourse.   
        inset_axes_side | float | inset axes side length as a fraction of the full figure.  
        arrow_length | float | lenth of arrow from data point to inset axes, as a fraction of the full figure.  
        figsize | tuple |  standard Matplotlib figsize tuple, in inches.  
        labels | dict or None | title for title, xlabel for x axis label, and ylabel for y axis label

    Returns:
        Interactive figure
    History:
        2020/09/09 | MEG | Modified from a sript in the ICASAR package.  
    
    """
    def remove_axes2_and_arrow(fig):
        """ Given a figure that has a second axes and an annotation arrow due to a 
        point having been hovered on, remove this axes and annotation arrow.  
        Inputs:
            fig | matplotlib figure 
        Returns:
        History:
            2020/09/08 | MEG | Written
        """
        # 1: try and remove any axes except the primary one
        try:
            fig.axes[1].remove()                
        except:
            pass
        
        # 2: try and remove any annotation arrows
        for art in axes1.get_children():
            if isinstance(art, matplotlib.patches.FancyArrow):
                try:
                    art.remove()        
                except:
                    continue
            else:
                continue
        fig.canvas.draw_idle()                                          # update the figure
    
    
    def axes_data_to_fig_percent(axes_lims, fig_lims, point):
        """ Given a data point, find where on the figure it plots (ie convert from axes coordinates to figure coordinates) 
        Inputs:
            axes_xlims | tuple | usually just the return of something like: axes1.get_ylim()
            fig_lims | tuple | the limits of the axes in the figure.  usuall (0.1, 0.9)  for an axes made with something like this: axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])                  # main axes
            point |float | point in data coordinates
        Returns:
            fig_position | float | where the data point is in the figure.  (0,0) would be the lower left corner.  
        History:
            2020/09/08 | MEG | Written
            
        """
        gradient = (fig_lims[1] - fig_lims[0])/(axes_lims[1] - axes_lims[0])
        y_intercept = fig_lims[0] - (gradient * axes_lims[0])
        fig_position = (gradient * point) + y_intercept
        return fig_position
    
    def calculate_insetaxes_offset(lims, points, offset_length):
        """
        The offsets between the inset axes and the point are different depending on which quadrant of the graph the point is in.  
        Inputs:
            lims | list | length is equal to the number of dimensions.  Filled with tuples of the axes limits.  
            point | list | length is equal to the number of diemsions. Filled with points.  
            offset_length | float | length of the arrow.  
        Returns:
            offsets | list | length is equal to the number of dimensions.  Length of offset for inset axes in each dimension.  
        History:
            2020/09/08 | MEG | Written
        """
        
        offsets = []
        for dim_n in range(len(lims)):                                        # loop through each dimension.  
            dim_centre = np.mean(lims[dim_n])
            if points[dim_n] < dim_centre:
                offsets.append(-offset_length)
            else:
                offsets.append(offset_length)
        return offsets
    
    
    def col_to_ma(col, pixel_mask):
        """ A function to take a column vector and a 2d pixel mask and reshape the column into a masked array.  
        Useful when converting between vectors used by BSS methods results that are to be plotted
        Inputs:
            col | rank 1 array | 
            pixel_mask | array mask (rank 2)
        Outputs:
            source | rank 2 masked array | colun as a masked 2d array
        """
        import numpy.ma as ma 
        import numpy as np
        
        source = ma.array(np.zeros(pixel_mask.shape), mask = pixel_mask )
        source.unshare_mask()
        source[~source.mask] = col.ravel()   
        return source

    def hover(event):
        if event.inaxes == axes1:                                                       # determine if the mouse is in the axes
            cont, ind = sc.contains(event)                                              # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
            if cont:                                                                    # if on point
                remove_axes2_and_arrow(fig)                                             # remove the axes and arrow created when hovering on the point (incase cursor moves from one point to next without going off a point)
                point_n = ind['ind'][0]                                                 # get the index of which data point we're hovering on in a simpler form.      
                
                # 1: Add the annotation arrow (from inset axes to data point)
                arrow_lengths = calculate_insetaxes_offset([axes1.get_xlim(), axes1.get_ylim()], 
                                                          [xy[0,point_n], xy[1,point_n]], arrow_length)                               # calculate the length of the arrow, which depends which quadrant we're in (as the arrow always go away from the plot)
                axes1.arrow(xy[0,point_n] + arrow_lengths[0], xy[1,point_n] + arrow_lengths[1],                                       # add the arrow.  Notation is all a bit backward as head is fixed at end, so it has to be drawn backwards.  
                            -arrow_lengths[0], -arrow_lengths[1], clip_on = False, zorder = 999)                                # clip_on makes sure it's visible, even if it goes off the edge of the axes.  

                # 2: Add the inset axes                
                fig_x = axes_data_to_fig_percent(axes1.get_xlim(), (0.1, 0.9), xy[0,point_n] + arrow_lengths[0])                   # convert position on axes to position in figure, ready to add the inset axes
                fig_y = axes_data_to_fig_percent(axes1.get_ylim(), (0.1, 0.9), xy[1,point_n] + arrow_lengths[1])                   # ditto for y dimension
                if arrow_lengths[0] > 0 and arrow_lengths[1] > 0:                                                                    # top right quadrant
                    inset_axes = fig.add_axes([fig_x, fig_y, inset_axes_side, inset_axes_side], anchor = 'SW')                       # create the inset axes, simple case, anochored to lower left forner
                elif arrow_lengths[0] < 0 and arrow_lengths[1] > 0:                                                                  # top left quadrant
                    inset_axes = fig.add_axes([fig_x - inset_axes_side, fig_y, inset_axes_side, inset_axes_side], anchor = 'SE')     # create the inset axes, nudged in x direction, anchored to lower right corner
                elif arrow_lengths[0] > 0 and arrow_lengths[1] < 0:                                                                  # lower right quadrant
                    inset_axes = fig.add_axes([fig_x, fig_y - inset_axes_side, inset_axes_side, inset_axes_side], anchor = 'NW')     # create the inset axes, nudged in y direction
                else:
                    inset_axes = fig.add_axes([fig_x - inset_axes_side, fig_y - inset_axes_side, inset_axes_side, inset_axes_side], anchor = 'NE')  # create the inset axes, nudged in both x and y
                
                # 3: Plot on the inset axes
                if temporal_data is not None:
                    inset_axes.plot(temporal_data['xvals'], temporal_data['time_courses'][point_n,])                            # draw the inset axes time course graph
                if spatial_data is not None:
                    inset_axes.imshow(col_to_ma(spatial_data['images_r2'][point_n,], spatial_data['mask']))                  # or draw the inset axes image
                inset_axes.set_xticks([])                                                                                       # and remove ticks (and so labels too) from x
                inset_axes.set_yticks([])                                                                                       # and from y
                fig.canvas.draw_idle()                                                                                          # update the figure.  
            else:                                                                       # else not on a point
                remove_axes2_and_arrow(fig)                                             # remove the axes and arrow created when hovering on the point                       
        else:                                                                           # else not in the axes
            remove_axes2_and_arrow(fig)                                                 # remove the axes and arrow created when hovering on the point (incase cursor moves from one point to next without going off a point)
    
    import matplotlib.pyplot as plt
    import matplotlib

    # norm = plt.Normalize(1,4)
    # cmap = plt.cm.RdYlGn
    if temporal_data is None and spatial_data is None:                                                                  # check inputs
        raise Exception("One of either spatial or temporal data must be supplied.  Exiting.  ")
    if temporal_data is not None and spatial_data is not None:
        raise Exception("Only either spatial or temporal data can be supplied, but not both.  Exiting.  ")
    
    fig = plt.figure(figsize = figsize)                                             # create the figure, size set in function args.  
    axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])                                      # main axes
    sc = axes1.scatter(xy[0,],xy[1,],c=c, s=100) #, cmap=cmap, norm=norm)

    # try:
    #     fig.canvas.set_window_title(labels['title'])
    #     fig.suptitle(labels['title'])
    # except:
    #     pass
    # try:
    #     axes1.set_xlabel(labels['xlabel'])
    # except:
    #     pass
    # try:
    #     axes1.set_ylabel(labels['ylabel'])
    # except:
    #     pass
           
    fig.canvas.mpl_connect("motion_notify_event", hover)                                # connect the figure and the function.  


#%% temporal data example

import numpy as np
import numpy.ma as ma

xy = np.random.rand(2,15)
time_courses = np.cumsum(np.random.randn(15,40), axis = 0)
xvals = np.arange(0,40)
c = np.random.randint(1,5,size=15)

temporal_data = {'time_courses' : time_courses,
                 'xvals'        : xvals}

labels = {'title' : 'Temporal Example',
          'xlabel' : 'x',
          'ylabel' : 'y'}

exploratory_2d_plot(xy, c, temporal_data = temporal_data, inset_axes_side = 0.1, arrow_length = 0.05, figsize = (10,6), labels = labels)


#%% spatial data example

spatial_maps = np.random.rand(15,100,100)

mask = np.where(np.random.randint(0,2, (100,100)) == 1, np.ones((100,100)), np.zeros((100,100))).astype(bool)

spatial_maps_ma = ma.array(spatial_maps, mask = np.repeat(mask[np.newaxis,], 15, axis = 0))

n_pixs = len(ma.compressed(spatial_maps_ma[0,]))
spatial_maps_r2 = np.zeros((15, n_pixs))

for i in range(15):
    spatial_maps_r2[i,] = ma.compressed(spatial_maps_ma[i,])
    

spatial_data = {'images_r2' : spatial_maps_r2,
                'mask'   : mask}
labels['title'] = 'Spatial Example'

exploratory_2d_plot(xy, c, spatial_data = spatial_data, inset_axes_side = 0.1, arrow_length = 0.05, figsize = (10,6), labels = labels)

#%% Old version as a script.  
#Version where the axes are drawn each time, using a more object orientated approach.  

# def remove_axes2_and_arrow(fig):
#     """ Given a figure that has a second axes and an annotation arrow due to a 
#     point having been hovered on, remove this axes and annotation arrow.  
#     Inputs:
#         fig | matplotlib figure 
#     Returns:
#     History:
#         2020/09/08 | MEG | Written
#     """
#     # 1: try and remove any axes except the primary one
#     try:
#         fig.axes[1].remove()                
#     except:
#         pass
    
#     # 2: try and remove any annotation arrows
#     for art in axes1.get_children():
#         if isinstance(art, matplotlib.patches.FancyArrow):
#             try:
#                 art.remove()        
#             except:
#                 continue
#         else:
#             continue
#     fig.canvas.draw_idle()                                          # update the figure


# def axes_data_to_fig_percent(axes_lims, fig_lims, point):
#     """ Given a data point, find where on the figure it plots (ie convert from axes coordinates to figure coordinates) 
#     Inputs:
#         axes_xlims | tuple | usually just the return of something like: axes1.get_ylim()
#         fig_lims | tuple | the limits of the axes in the figure.  usuall (0.1, 0.9)  for an axes made with something like this: axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])                  # main axes
#         point |float | point in data coordinates
#     Returns:
#         fig_position | float | where the data point is in the figure.  (0,0) would be the lower left corner.  
#     History:
#         2020/09/08 | MEG | Written
        
#     """
#     gradient = (fig_lims[1] - fig_lims[0])/(axes_lims[1] - axes_lims[0])
#     y_intercept = fig_lims[0] - (gradient * axes_lims[0])
#     fig_position = (gradient * point) + y_intercept
#     return fig_position

# def calculate_insetaxes_offset(lims, points, offset_length):
#     """
#     The offsets between the inset axes and the point are different depending on which quadrant of the graph the point is in.  
#     Inputs:
#         lims | list | length is equal to the number of dimensions.  Filled with tuples of the axes limits.  
#         point | list | length is equal to the number of diemsions. Filled with points.  
#         offset_length | float | length of the arrow.  
#     Returns:
#         offsets | list | length is equal to the number of dimensions.  Length of offset for inset axes in each dimension.  
#     History:
#         2020/09/08 | MEG | Written
#     """
    
#     offsets = []
#     for dim_n in range(len(lims)):                                        # loop through each dimension.  
#         dim_centre = np.mean(lims[dim_n])
#         if points[dim_n] < dim_centre:
#             offsets.append(-offset_length)
#         else:
#             offsets.append(offset_length)
#     return offsets


# import matplotlib.pyplot as plt
# import matplotlib
# import numpy as np; np.random.seed(1)

# inset_axes_side = 0.1
# arrow_length = 0.1

# x = np.random.rand(15)
# y = np.random.rand(15)
# time_courses = np.cumsum(np.random.randn(15,40), axis = 0)
# xvals = np.arange(0,40)
# c = np.random.randint(1,5,size=15)
# norm = plt.Normalize(1,4)
# cmap = plt.cm.RdYlGn

# fig = plt.figure()
# axes1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])                  # main axes
# sc = axes1.scatter(x,y,c=c, s=100, cmap=cmap, norm=norm)

# def hover(event):
#     if event.inaxes == axes1:                                              # determine if the mouse is in the axes
#         print('in axes', end = '')
#         cont, ind = sc.contains(event)                                              # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
#         if cont:
#             print('on point')
#             remove_axes2_and_arrow(fig)                                                 # remove the axes and arrow created when hovering on the point (incase cursor moves from one point to next without going off a point)
#             point_n = ind['ind'][0]          
#             arrow_lengths = calculate_insetaxes_offset([axes1.get_xlim(), axes1.get_ylim()], 
#                                                       [x[point_n], y[point_n]], arrow_length)
#             axes1.arrow(x[point_n] + arrow_lengths[0], y[point_n] + arrow_lengths[1],                                       # add the arrow.  Notation is all a bit backward as head is fixed at end, so it has to be drawn backwards.  
#                         -arrow_lengths[0], -arrow_lengths[1], clip_on = False, zorder = 999)                                # clip_on makes sure it's visible, even if it goes off the edge of the axes.  
            
#             fig_x = axes_data_to_fig_percent(axes1.get_xlim(), (0.1, 0.9), x[point_n] + arrow_lengths[0])                   # convert position on axes to position in figure, ready to add the inset axes
#             fig_y = axes_data_to_fig_percent(axes1.get_ylim(), (0.1, 0.9), y[point_n] + arrow_lengths[1])                   # ditto for y dimension
            
#             print(arrow_lengths)
#             if arrow_lengths[0] > 0 and arrow_lengths[1] > 0:                                                               # top right quadrant
#                 inset_axes = fig.add_axes([fig_x, fig_y, inset_axes_side, inset_axes_side])                                         # create the inset axes, simple case
#             elif arrow_lengths[0] < 0 and arrow_lengths[1] > 0:                                                             # top left quadrant
#                 inset_axes = fig.add_axes([fig_x - inset_axes_side, fig_y, inset_axes_side, inset_axes_side])                                         # create the inset axes, simple case
#             elif arrow_lengths[0] > 0 and arrow_lengths[1] < 0:                                                             # lower right quadrant
#                 inset_axes = fig.add_axes([fig_x, fig_y - inset_axes_side, inset_axes_side, inset_axes_side])                                         # create the inset axes, simple case
#             else:
#                 inset_axes = fig.add_axes([fig_x - inset_axes_side, fig_y - inset_axes_side, inset_axes_side, inset_axes_side])                                         # create the inset axes, simple case
                
#             inset_axes.plot(xvals, time_courses[point_n,])                                                                  # draw the inset axes figure
#             inset_axes.set_xticks([])                                                                                       # and remove ticks (and so labels too) from x
#             inset_axes.set_yticks([])                                                                                       # and from y
#             fig.canvas.draw_idle()                                                                                          # update the figure.  
#         else:
#             print('  off point')
#             remove_axes2_and_arrow(fig)                                                 # remove the axes and arrow created when hovering on the point

                    
#     else:
#         print('out of axes')
#         remove_axes2_and_arrow(fig)                                                 # remove the axes and arrow created when hovering on the point (incase cursor moves from one point to next without going off a point)
        

# fig.canvas.mpl_connect("motion_notify_event", hover)

