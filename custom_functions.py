
def nd2_metadata(imported_nd2, trackmate_z_spacing=0.5, trackmate_micron_per_px=0.1405):
    min_z = np.min(raw.metadata['z_coordinates'])
    max_z = np.max(raw.metadata['z_coordinates'])
    z_range = max_z - min_z
    z_slices = np.max(raw.metadata['z_levels'])+1
    
    metadata_dict = {
    'channel':imported_nd2.metadata['channels'],
#    'frames':imported_nd2.metadata['frames'],
    'num_frames':imported_nd2.metadata['num_frames'],
    's_per_frame':imported_nd2.metadata['experiment']['loops'][0]['sampling_interval']/1000,
    's_total':imported_nd2.timesteps[-1]/1000,
    'micron_per_px':imported_nd2.metadata['pixel_microns'],
    'min_z':min_z,
    'max_z':max_z,
    'z_range':max_z - min_z,
    'z_slices':z_slices,
    'z_spacing':z_range/z_slices}
    
    #metadata = pd.DataFrame(metadata_dict).set_index('channel')
    
    return metadata_dict

def movie_3d(stack, channels_list, vis_params, points=None, traj=None,
             z_low=None, z_high=None, y_low=None, y_high=None, x_low=None, x_high=None,
             script=None, movie_name='default_movie', 
             merge=True, separate=True, close_window=True):
    if merge == True:
        viewer = napari.Viewer(ndisplay=3)

        for channel in channels_list:
            vis_channel = vis_params.loc[channel]
            stack_channel = stack[...,vis_channel['stack_index']]
            if z_low and z_high:
                new_stack_channel = np.zeros_like(stack_channel)
                new_stack_channel[:,z_low:z_high,...] = stack_channel[:,z_low:z_high,...]
                stack_channel = new_stack_channel
            scale_channel = [1,vis_channel['z_spacing'],
                             vis_channel['micron_per_px'],
                             vis_channel['micron_per_px']]
            scale_trackmate = [1,vis_channel['micron_per_px'],
                              vis_channel['micron_per_px'],
                              vis_channel['micron_per_px']]

            viewer.add_image(stack_channel, colormap=vis_channel['color'],
                         name=vis_channel['label'], scale=scale_channel,
                         contrast_limits=vis_channel['contrast'], blending='additive')
            
        if type(points) != type(None):
            viewer.add_points(points, face_color='#ffffff00', edge_color='white', 
                  size=5, edge_width=3, opacity=0.7, scale=scale_trackmate,
                  symbol='disc', name='vesicle')

        if type(traj) != type(None):
            viewer.add_shapes(traj, shape_type='path', edge_color='white',
                              name='track', scale=scale_trackmate[1:],
                  edge_width=0.5, opacity = 0.7)

        if script:
            movie = Movie(myviewer=viewer)
            movie.create_state_dict_from_script(script)
            movie.make_movie(name = 'movies/'+movie_name+'_merge.mp4',
                             resolution = 320, fps = 20)
        if close_window == True:
            viewer.close()
    
    if separate == True:
        for channel in channels_list:
            viewer = napari.Viewer(ndisplay=3)
            vis_channel = vis_params.loc[channel]
            stack_channel = stack[...,vis_channel['stack_index']]
            if z_low and z_high:
                new_stack_channel = np.zeros_like(stack_channel)
                new_stack_channel[:,z_low:z_high,...] = stack_channel[:,z_low:z_high,...]
                stack_channel = new_stack_channel
            scale_channel = [1,vis_channel['z_spacing'],
                             vis_channel['micron_per_px'],
                             vis_channel['micron_per_px']]
            scale_trackmate = [1,vis_channel['micron_per_px'],
                              vis_channel['micron_per_px'],
                              vis_channel['micron_per_px']]

            viewer.add_image(stack_channel, 
                         name=vis_channel['label'], scale=scale_channel,
                         contrast_limits=vis_channel['contrast'], blending='additive')

            if type(points) != type(None):
                viewer.add_points(points, face_color='#ffffff00', edge_color='red', 
                      size=5, edge_width=3, opacity=0.7, scale=scale_trackmate,
                      symbol='disc', name='vesicle')
                
            if type(traj) != type(None):
                viewer.add_shapes(traj, shape_type='path', edge_color='yellow',
                                  name='track', scale=scale_trackmate[1:],
                      edge_width=0.5, opacity = 0.7)

            if script:
                movie = Movie(myviewer=viewer)
                movie.create_state_dict_from_script(script)
                movie.make_movie(name = 'movies/'+movie_name+'_'+vis_channel['label']+'.mp4',
                                 resolution = 320, fps = 20)
            if close_window == True:
                viewer.close()
                
                
def findInflectionPoint(x,y,debug=False):
    # this function calculates the point where the function has a kink and returns the index of that point
    x = x.loc[:np.argmax(y)]
    y = y.loc[:np.argmax(y)]
    slope, intercept, r, prob2, see = stats.linregress(x, y)
    mx = x.mean()
    sx2 = ((x-mx)**2).sum()
    sd_intercept = see * np.sqrt(1./len(x) + mx*mx/sx2)
    sd_slope = see * np.sqrt(1./sx2)

    inter_slope_intercept = [slope,sd_slope,intercept,sd_intercept]

    if debug == True:
        plt.plot(x, y, 'o', label='original data',color='b')
        plt.plot(x, intercept + slope*x, 'r', label='fitted line')
        print("slope={}±{}, intercept={}±{}, r={}, prob2={}, see={}".format(slope,sd_slope, intercept,sd_intercept, r, prob2, see))
    

# find the inflection point

    newy = y-(intercept+slope*x)
    posdiff = newy.loc[y > np.max(y)*0.5]
    #newy = np.array(newy)
    #inflectionPoint = np.min(newy)
    peak = posdiff.idxmax()
    inflectionPointIndex = np.argmin(newy.loc[:peak])
    inflectionPointIndex = newy.loc[:peak].idxmin()
    
    if debug == True:
        plt.vlines(x[inflectionPointIndex], min(y), max(y), label='inflection')
        plt.legend()
        plt.show()
    return(inflectionPointIndex)