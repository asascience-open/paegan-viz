from paegan.viz.trajectory import CFTrajectory

traj = CFTrajectory("resources/trajectories.nc")
traj.plot_animate('my_test_animation.mp4', view=(0,-90))
