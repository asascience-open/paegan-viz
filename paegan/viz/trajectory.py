import numpy as np
import netCDF4, sys, os
import matplotlib
import matplotlib.pyplot
from matplotlib import cm, animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
import subprocess
import multiprocessing
from paegan.transport.shoreline import Shoreline
from shapely.geometry import Point
from cv2 import VideoWriter, imread
from cv import CV_FOURCC
    
semo = multiprocessing.Semaphore(multiprocessing.cpu_count() - 1)
    
class CFTrajectory(object):
    def __init__(self, filepath):
        self.nc = netCDF4.Dataset(filepath)
        self.path = filepath
        
    def plot_summary(self, view=(25, -75), bathy=os.path.join(__file__,"../../resources/bathymetry/ETOPO1_Bed_g_gmt4.grd")):
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
    
    def plot_animate(self, output, temp_folder=None, view=(45, -75), bathy=os.path.join(__file__,"../../resources/bathymetry/ETOPO1_Bed_g_gmt4.grd"),
                     frame_prefix='_paegan', extent=None, stride=None):
       
        if temp_folder == None:
            temp_folder = os.path.dirname(output)

        if extent == None:
            visual_bbox = (self.nc.variables['lon'][:,0].min()-.6, self.nc.variables['lat'][:,0].min()-.75,
                           self.nc.variables['lon'][:,0].max()+.6, self.nc.variables['lat'][:,0].max()+.75)#tracks.buffer(1).bounds
        else:
            visual_bbox = extent
            
        pt = Point(((visual_bbox[2]-visual_bbox[0])/2)+visual_bbox[0],((visual_bbox[3]-visual_bbox[1])/2)+visual_bbox[1])
        coast_line = Shoreline(point=pt, spatialbuffer=1.5).linestring
        c_lons, c_lats = coast_line.xy
        c_lons = np.array(c_lons)
        c_lats = np.array(c_lats)
        c_lons = np.where((c_lons >= visual_bbox[0]) & (c_lons <= visual_bbox[2]), c_lons, np.nan)
        c_lats = np.where((c_lats >= visual_bbox[1]) & (c_lats <= visual_bbox[3]), c_lats, np.nan)               
        #add bathymetry
        if stride == None:
            if visual_bbox[2] - visual_bbox[0] < 1.5:
                stride = 1
            else:
                stride = 2
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
        bath[bath>0] = 0
        #bath = bath.astype(np.float32)
        bath[bath<-800] = -800#np.nan
        x_grid, y_grid = np.meshgrid(lons, lats)

        mpl_extent = matplotlib.transforms.Bbox.from_extents(visual_bbox[0],visual_bbox[1],visual_bbox[2],visual_bbox[3])
            
        CNorm = matplotlib.colors.Normalize(vmin=-400,
                                            vmax=300,
                                            )
        '''                                   
        #bath[bath>0] = np.log(bath[bath>0]) * 10                                  
        s = ax2.plot_surface(x_grid, y_grid, bath, rstride=stride, cstride=stride,
                cmap="Blues_r",  linewidth=.01, antialiased=False,
                norm=CNorm, shade=True, edgecolor='#6183A6')
                
        ax2.set_xlim3d(visual_bbox[0],visual_bbox[2])
        ax2.set_ylim3d(visual_bbox[1],visual_bbox[3])
        ax2.view_init(*view)
        
        #ax2.set_zmargin(50)
        ax2.set_xlabel('Longitude')
        ax2.set_ylabel('Latitude')
        ax2.set_zlabel('Depth (m)')
        #ax2.set_frame_on(False)
        #ax2.set_position([0,0,1,1])
        ax2.xaxis.set_ticklabels([])
        ax2.yaxis.set_ticklabels([])
        ax2.zaxis.set_ticklabels(['Surface'])
        ax2.zaxis.set_ticks([0])
        ax2.grid(False)
        #ax2.set_zlim(-200, 100)
        '''
        mgr = multiprocessing.Manager()
        fname = mgr.list()

        p_proj_lons = []
        p_proj_lats = []
        p_proj_depth = []
        
        lat = self.nc.variables['lat'][:,:]
        lon = self.nc.variables['lon'][:,:]
        depth = self.nc.variables['depth'][:,:]
        time = netCDF4.num2date(self.nc.variables['time'][:], self.nc.variables['time'].units)
            
        datetimeformat = '%Y-%m-%d %H:%M'

        p = []
        c = 0
        length = self.nc.variables['particle'].shape[0]
            
        def render_frame(visual_bbox, c_lons, c_lat, mpl_extent, 
                         x_grid, y_grid, bath, stride, view,
                         length, lat, lon, depth, temp_folder,
                         frame_prefix, c, semo):
            import numpy as np
            import netCDF4, sys, os
            import matplotlib
            import matplotlib.pyplot
            from matplotlib import cm, animation
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib.ticker import MultipleLocator
            with semo:
                fig2 = matplotlib.pyplot.figure(figsize=(12,6))
                ax2 = fig2.add_subplot(111, projection='3d')
                ax3 = fig2.add_axes([.75, .1, .15, .3])
                ax4 = fig2.add_axes([.2, .1, .15, .3])
                subbox = visual_bbox#(self.nc.variables['lon'][:,0].min(), self.nc.variables['lat'][:,0].min(),
                         #self.nc.variables['lon'][:,0].max(), self.nc.variables['lat'][:,0].max())
                ax3.plot(c_lons, c_lats, clip_box=mpl_extent, clip_on=True, color='c') # shoreline
               
                ax3.set_xlim(subbox[0],subbox[2])
                ax3.set_ylim(subbox[1],subbox[3])
                #ax3.pcolor(x_grid, y_grid, bath, cmap="Blues_r", norm=CNorm)

                ax2.plot_surface(x_grid, y_grid, bath, rstride=stride, cstride=stride,
                    cmap="Blues_r",  linewidth=0.01, antialiased=False,
                    norm=CNorm, shade=True, edgecolor='#6183A6')
                
                ax2.plot(c_lons, c_lats, np.zeros_like(c_lons))
                
                ax2.set_xlim3d(visual_bbox[0],visual_bbox[2])
                ax2.set_ylim3d(visual_bbox[1],visual_bbox[3])
                ax2.view_init(*view)
                
                ax2.set_title(time[i].strftime(datetimeformat) + " - " + time[i+2].strftime(datetimeformat))
                #ax2.set_zmargin(50)
                ax3.set_xlabel('Longitude')
                ax3.set_ylabel('Latitude')
                ax3.tick_params(axis='both', which='major', labelsize=10)
                ax3.yaxis.set_ticks_position('right')
                ax3.ticklabel_format(axis='x', style='plain')
                ax3.xaxis.set_major_locator(MultipleLocator(.5))
                ax3.grid(True)
                ax4.set_ylabel('Depth (m)')
                ax4.set_ylim(-200, 0)
                ax4.tick_params(axis='x', which='major', labelsize=10)
                #ax4.set_xlim(1,3)
                ax4.xaxis.set_ticklabels([])
                #ax3.xaxis.set_ticklabels(np.unique(c_lons.astype(int)))
                ax2.set_zlabel('Depth (m)')
                #ax2.set_frame_on(False)
                #ax2.set_position([0,0,1,1])
                ax2.xaxis.set_ticklabels([])
                ax2.yaxis.set_ticklabels([])
                #ax2.zaxis.set_ticklabels(['Surface'])
                ax2.zaxis.set_ticks(range(-800,100,200))
                ax2.grid(False)
                #ax2.set_zlim(-200, 100)
                
                for j in range(length):
                    # Each particle
                    ax2.plot(lon[i:i+3,j], lat[i:i+3,j], depth[i:i+3,j], ':', c='r', linewidth=2, markersize=5, markerfacecolor='r')
                    ax3.plot(lon[:i+3,j], lat[:i+3,j], c='.2',
                         linewidth=.5, markersize=5, markerfacecolor='r',)
                    ax3.scatter(lon[i+2,j], lat[i+2,j], c='r')
                    ax4.plot(range(i+3), depth[:i+3,j], c='r', linewidth=.5, aa=True)
                    ax4.scatter(np.ones_like(depth[i+2,j])*(i+2), depth[i+2,j], c='r')
                    if i == 2:
                        ax4.set_xlim(i-2,i+2.25)
                    elif i >= 3:
                        ax4.set_xlim(i-3,i+2.25)
                    else:
                        ax4.set_xlim(i,i+2.25)
                #ax2.scatter(lon[i,:], lat[i,:], depth[i,:], zdir='z', c='r')
                ax2.set_zlim3d(-800,25)
                
                image_path = os.path.join(temp_folder,'%s%04d.png' % (frame_prefix, c))
                fig2.savefig(image_path, dpi=350, bbox_inches='tight')
                fname.append(image_path)

                del ax2, ax3, ax4, subbox, fig2
        
        jobs = []
        for i in range(self.nc.variables['time'].shape[0])[:-4:2]:
            p = multiprocessing.Process(target=render_frame, args=(visual_bbox, c_lons, c_lats, mpl_extent, 
                                                                   x_grid, y_grid, bath, stride, view,
                                                                   length, lat, lon, depth, temp_folder,
                                                                   frame_prefix, c, semo)
                                       )
            p.start()
            jobs.append(p)
            c += 1
        
        for j in jobs:
            j.join(120)

        return save_animation(output, sorted(fname), frame_prefix=frame_prefix)
        
def save_animation(filename, files, fps=10, codec=None, clear_temp=True, frame_prefix='_tmp'):

    '''
    Saves a movie file by drawing every frame.

    *filename* is the absolute path to the output filename, eg :file:`mymovie.mp4`

    *files* is a list of files to add to the movie

    *fps* is the frames per second in the movie

    *codec* is the codec to be used,if it is supported by the output method.

    *clear_temp* specifies whether the temporary image files should be
    deleted.
    '''
    if len(files) > 0:

        fps = max(fps,10)

        if codec is None:
            #codec = CV_FOURCC('D','I','B',' ')
            #codec = CV_FOURCC('D','I','V','X')
            codec = CV_FOURCC('X','V','I','D')
            #codec = CV_FOURCC('X','2','6','4')

        # To get correct width and height for video
        height,width,bands = imread(files[0]).shape
        vw = VideoWriter(filename, codec, fps, (width, height), True)

        if vw is None:
            print "Error creating video writer"

        for fname in files:

            # 2.0
            ig = imread(fname)
            vw.write(ig)
            vw.write(ig)

            if clear_temp:
                os.remove(fname)

        del vw
        return True

    return False
