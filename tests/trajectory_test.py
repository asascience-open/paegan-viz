import sys

sys.path.append("/media/sf_Python/paegan")


from paegan.viz.trajectory import CFTrajectory

traj = CFTrajectory("/media/sf_Python/paegan/model_run_output.nc")

traj.plot_animate('my_test_animation.mp4', view=(0,-90))
