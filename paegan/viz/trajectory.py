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
    
    def plot_animate(self, output, view=(45, -75), bathy=os.path.join(__file__,"../../resources/bathymetry/ETOPO1_Bed_g_gmt4.grd"),
                     frame_prefix='_paegan', extent=None, stride=None):
       
        
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
        fname = []
        p_proj_lons = []
        p_proj_lats = []
        p_proj_depth = []
        
        lat = self.nc.variables['lat'][:,:]
        lon = self.nc.variables['lon'][:,:]
        depth = self.nc.variables['depth'][:,:]
        time = netCDF4.num2date(self.nc.variables['time'][:], self.nc.variables['time'].units)
        def create_image(ax2, ax3, ax4, i, lat, lon, depth, prefix, c):
            for q in [0,1]:
                fname.append('%s%04d.png' % (frame_prefix, c))
                for j in range(self.nc.variables['particle'].shape[0]):
                    line, = ax2.plot(lon[i:i+3,j],
                                     lat[i:i+3,j],
                                     depth[i:i+3,j], ':', c='r',
                                     linewidth=2, markersize=5, markerfacecolor='r',) # each particle)
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
                    line.set_markevery((i-3,1))
                #ax2.scatter(lon[i,:], lat[i,:], depth[i,:], zdir='z', c='r')
                ax2.set_zlim3d(-800,25)
                
                fig2.savefig(fname[c], dpi=350, bbox_inches='tight')
                c += 1
            return c
            
        p = []
        c = 0
        for i in range(self.nc.variables['time'].shape[0])[:-4:2]:
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
            s = ax2.plot_surface(x_grid, y_grid, bath, rstride=stride, cstride=stride,
                cmap="Blues_r",  linewidth=0.01, antialiased=False,
                norm=CNorm, shade=True, edgecolor='#6183A6')
            ax2.plot(c_lons, c_lats, np.zeros_like(c_lons))
            
            
            ax2.set_xlim3d(visual_bbox[0],visual_bbox[2])
            ax2.set_ylim3d(visual_bbox[1],visual_bbox[3])
            ax2.view_init(*view)
            datetimeformat = '%y/%m/%d %H:%M'
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
            c = create_image(ax2, ax3, ax4,  i, lat, lon, depth, frame_prefix, c)

        save_animation(output, fname, frame_prefix=frame_prefix)
        
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
               'mf://%s*.png'%frame_prefix,
               '-mf',
               'type=png:w=600:h=300:fps='+str(fps),
               '-ovc',
               'lavc',
               '-lavcopts',
               #'vcodec=mpeg4',
               'vcodec=msmpeg4v2:vbitrate=1000', 
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
    
