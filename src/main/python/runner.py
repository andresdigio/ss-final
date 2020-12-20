import os

jar_path = '/Users/andres/itba/ss-final/target/ss-final-1.0-jar-with-dependencies.jar'

def run_simulation(T, dt, o, path):
    cmd = 'java -jar {:s} -data -T {:g} -dt {:g} -o {:g} -out {:s}'.format(jar_path, T, dt, o, path)
    print(cmd)
    os.system(cmd)

for o in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    for i in [1, 2, 3, 4, 5]:
        run_simulation(1500, 5e-5, o, '{:d}/{:d}.data'.format(int(o*100), i))
