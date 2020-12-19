import os

jar_path = 'C:/Users/Andres/ss-final/target/ss-final-1.0-jar-with-dependencies.jar'

def run_simulation(ST, dt, o, path):
    cmd = 'java -jar {:s} -data -ST {:g} -dt {:g} -o {:g} -out {:s}'.format(jar_path, ST, dt, o, path)
    print(cmd)
    os.system(cmd)

for o in [0.5, 0.7, 0.9]:
    for i in [1]:
    #for i in [1, 2, 3, 4, 5]:
        run_simulation(1000, 5e-5, 0.5, '{:d}/{:d}.data'.format(int(o*100), i))
