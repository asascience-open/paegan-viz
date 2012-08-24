import numpy as np
import netCDF4, sys, os
import matplotlib
import matplotlib.pyplot
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import subprocess
import multiprocessing
    
class CFTrajectory(object):
    def __init__(self, filepath):
        self.nc = netCDF4.Dataset(filepath)
        self.path = filepath
    def plot_summary(self, view=(25, -75), bathy="/media/sf_Python/paegan/paegan/resources/bathymetry/ETOPO1_Bed_g_gmt4.grd"):
        fig = matplotlib.pyplot.figure()#figsize=(20,16)) # call a blank figure
        ax = fig.gca(projection='3d') # line with points
        
        for i in range(self.nc.variables['lat'].shape[1]):
            p_proj_lons = self.nc.variables['lon'][:,i]
            p_proj_lats = self.nc.variables['lat'][:,i]
            ax.plot(p_proj_lons,
                    p_proj_lats,
                    self.nc.variables['depth'][:,i],
                    marker='o', c='red') # each particle

        visual_bbox = (p_proj_lons.min()-.8, p_proj_lats.min()-.8,
                       p_proj_lons.max()+.8, p_proj_lats.max()+.8)#tracks.buffer(1).bounds
        
        #coast_line = Shoreline(point=midpoint, spatialbuffer=1.5).linestring
        #c_lons, c_lats = coast_line.xy
        #c_lons = np.array(c_lons)
        #c_lats = np.array(c_lats)
        #c_lons = np.where((c_lons >= visual_bbox[0]) & (c_lons <= visual_bbox[2]), c_lons, np.nan)
        #c_lats = np.where((c_lats >= visual_bbox[1]) & (c_lats <= visual_bbox[3]), c_lats, np.nan)

        #add bathymetry
        nc1 = netCDF4.Dataset(os.path.normpath(bathy))
        x = nc1.variables['x']
        y = nc1.variables['y']

        x_indexes = np.where((x[:] >= visual_bbox[0]) & (x[:] <= visual_bbox[2]))[0]
        y_indexes = np.where((y[:] >= visual_bbox[1]) & (y[:] <= visual_bbox[3]))[0]

        x_min = x_indexes[0] 
        x_max = x_indexes[-1]
        y_min = y_indexes[0]
        y_max = y_indexes[-1]

        lons = x[x_min:x_max]
        lats = y[y_min:y_max]
        bath = nc1.variables['z'][y_min:y_max,x_min:x_max]

        x_grid, y_grid = np.meshgrid(lons, lats)

        mpl_extent = matplotlib.transforms.Bbox.from_extents(visual_bbox[0],visual_bbox[1],visual_bbox[2],visual_bbox[3])
        
        CNorm = matplotlib.colors.Normalize(vmin=-200,
                                            vmax=500,
                                            )
        s = ax.plot_surface(x_grid,y_grid,bath, rstride=1, cstride=1,
            cmap="gist_earth", shade=True, linewidth=0, antialiased=False,
            norm=CNorm, edgecolors=None)
            #edgecolors=None,) # bathymetry
        fig.colorbar(s)
        #ax.plot(c_lons, c_lats, clip_box=mpl_extent, clip_on=True, color='c') # shoreline
        ax.set_xlim3d(visual_bbox[0],visual_bbox[2])
        ax.set_ylim3d(visual_bbox[1],visual_bbox[3])
        ax.set_zmargin(0.1)
        ax.view_init(*view)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_zlabel('Depth (m)')
        #matplotlib.pyplot.show()
        fig.savefig('trajectory.png')
    
    def plot_animate(self, output, view=(25, -75), bathy="/media/sf_Python/paegan/paegan/resources/bathymetry/ETOPO1_Bed_g_gmt4.grd",
                     frame_prefix='_paegan'):
        fig2 = matplotlib.pyplot.figure(figsize=(6,6)) # call a blank figure
        #ax2 = fig2.gca(projection='3d') # line with points
        #ax2 = fig2.add_axes(rect=(0,1,0,1), projection='3d')
        ax2 = fig2.add_subplot(111, projection='3d')
        visual_bbox = (self.nc.variables['lon'][:,0].min()-.15, self.nc.variables['lat'][:,0].min()-.15,
                       self.nc.variables['lon'][:,0].max()+.15, self.nc.variables['lat'][:,0].max()+.15)#tracks.buffer(1).bounds
                           
        #add bathymetry
        nc1 = netCDF4.Dataset(os.path.normpath(bathy))
        x = nc1.variables['x']
        y = nc1.variables['y']

        x_indexes = np.where((x[:] >= visual_bbox[0]) & (x[:] <= visual_bbox[2]))[0]
        y_indexes = np.where((y[:] >= visual_bbox[1]) & (y[:] <= visual_bbox[3]))[0]

        x_min = x_indexes[0] 
        x_max = x_indexes[-1]
        y_min = y_indexes[0]
        y_max = y_indexes[-1]

        lons = x[x_min:x_max]
        lats = y[y_min:y_max]
        bath = nc1.variables['z'][y_min:y_max,x_min:x_max]

        x_grid, y_grid = np.meshgrid(lons, lats)

        mpl_extent = matplotlib.transforms.Bbox.from_extents(visual_bbox[0],visual_bbox[1],visual_bbox[2],visual_bbox[3])
            
        CNorm = matplotlib.colors.Normalize(vmin=-200,
                                            vmax=300,
                                            )
                                            
        #bath[bath>0] = np.log(bath[bath>0]) * 10                                  
        s = ax2.plot_surface(x_grid, y_grid, bath, rstride=1, cstride=1,
                cmap="gist_earth",  linewidth=0, antialiased=False,
                norm=CNorm, shade=True, edgecolor='.2')
                
        ax2.set_xlim3d(visual_bbox[0],visual_bbox[2])
        ax2.set_ylim3d(visual_bbox[1],visual_bbox[3])
        ax2.view_init(*view)
        
        ax2.set_zmargin(1)
        #ax2.set_xlabel('Longitude')
        #ax2.set_ylabel('Latitude')
        #ax2.set_zlabel('Depth (m)')
        #ax2.set_frame_on(False)
        #ax2.set_position([0,0,1,1])
        fname = []
        p_proj_lons = []
        p_proj_lats = []
        p_proj_depth = []
        
        lat = self.nc.variables['lat'][:,:]
        lon = self.nc.variables['lon'][:,:]
        depth = self.nc.variables['depth'][:,:]
        def create_image(ax2, i, lat, lon, depth, prefix, c):
            for q in [0,1]:
                fname.append('%s%04d.png' % (frame_prefix, c))
                
                for j in range(self.nc.variables['particle'].shape[0]):
                    line, = ax2.plot(lon[:i+1,j],
                             lat[:i+1,j],
                             depth[:i+1,j], ':', c='gray',
                             linewidth=.5, markersize=5, markerfacecolor='r',) # each particle)
                    line.set_markevery((i-3,1))
                #ax2.set_zlim3d(-200,0)
                fig2.savefig(fname[c], dpi=100, bbox_inches='tight')
                c += 1
                
            return c
            
        p = []
        c = 0
        for i in range(self.nc.variables['time'].shape[0]):
            c = create_image(ax2, i, lat, lon, depth, frame_prefix, c)
            
            #p.append(multiprocessing.Process(target=create_image, args=(ax2, i, fname[i], lat, lon, depth)))
            #p[i].start()
        #[proc.join() for proc in p]
            
        save_animation(output, fname, frame_prefix='_paegan')
        
def save_animation(filename, fnames, fps=10, codec='mpeg4', clear_temp=True,
    frame_prefix='_tmp'):
    '''
    Saves a movie file by drawing every frame.

    *filename* is the output filename, eg :file:`mymovie.mp4`

    *fps* is the frames per second in the movie

    *codec* is the codec to be used,if it is supported by the output method.

    *clear_temp* specifies whether the temporary image files should be
    deleted.

    *frame_prefix* gives the prefix that should be used for individual
    image files.  This prefix will have a frame number (i.e. 0001) appended
    when saving individual frames.
    '''
    from subprocess import Popen, PIPE
    #command = ['ffmpeg', '-y', '-r', str(fps), '-b', '1800k', '-i',
    #    '%s%%04d.png' % frame_prefix, filename]
    
    command = ('mencoder',
               'mf://*.png',
               '-mf',
               'type=png:w=400:h=400:fps='+str(fps),
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=mpeg4',
               '-oac',
               'copy',
               '-o',
               filename)
    proc = Popen(command, shell=False,
        stdout=PIPE, stderr=PIPE)
    proc.wait()

    
    #Delete temporary files
    if clear_temp:
        for fname in fnames:
            os.remove(fname)
    
