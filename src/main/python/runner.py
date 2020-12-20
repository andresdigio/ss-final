import os

jar_path = 'C:/Users/Andres/ss-final/target/ss-final-1.0-jar-with-dependencies.jar'

def run_simulation(N, T, dt, o, path):
    cmd = 'java -jar {:s} -data -n {:d} -T {:g} -dt {:g} -o {:g} -out {:s}'.format(jar_path, N, T, dt, o, path)
    print(cmd)
    os.system(cmd)


for n in [10, 15]:
    dir = 'n{:d}'.format(n)
    if not os.path.exists('data/' + dir):
        os.mkdir('data/' + dir)
    for o in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        odir = dir + '/{:d}'.format(int(o*100))
        if not os.path.exists('data/' + odir):
            os.mkdir('data/' + odir)
        for i in [1, 2, 3, 4, 5]:
            run_simulation(n, 1500, 5e-5, o, odir + '/{:d}.data'.format(i))
