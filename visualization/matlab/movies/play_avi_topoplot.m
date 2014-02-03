file = '/mdkm1/2/kmdev/projects/jsherfey/sleep/movies/s4_18_SlowOscillations_500frametarget.avi';
data = aviread(file);

fps  = 10; % frames per second (dt=100ms ==> fps=1/dt=10)
movie(data,1,fps);

