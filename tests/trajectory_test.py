import os
import unittest

from paegan.viz.trajectory import CFTrajectory, save_animation

class TrajectoryTests(unittest.TestCase):
        
    def test_full_animation(self):
        traj = CFTrajectory(os.path.join(os.path.dirname(__file__), "resources/trajectories.nc"))
        output_movie = os.path.join(os.path.dirname(__file__), "output/from_netcdf.avi")
        status = traj.plot_animate(output_movie, bathy="/home/dev/Development/paegan/paegan/resources/bathymetry/ETOPO1_Bed_g_gmt4.grd")
        assert status is True
        
    def test_save_animation(self):
        root_dir = os.path.join(os.path.dirname(__file__), "resources/trajectory_images")
        files = sorted(os.listdir(root_dir))
        files = map(lambda x: os.path.join(root_dir,x),files)
        #print files
        output_movie = os.path.join(os.path.dirname(__file__), "output/from_files.avi")
        status = save_animation(output_movie, files, clear_temp=False)
        assert status is True