package simulation;

import analysis.DataExporter;
import model.FastAtan;
import model.Particle;
import ovito.Ovito;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import java.util.*;
import java.util.stream.Collectors;
import org.apache.commons.cli.*;


public class Orbit {
    private static final double SUN_MASS = 5000;
    private static final double SUN_RADIUS = 50;
    private static final double PARTICLE_MASS = 1;
    private static final double PARTICLE_RADIUS = 2;
    private static final double SPAWN_DISTANCE = 300;
    private static final double VT0 = 13;
    private static final double VN0 = 2;
    public static final double G = 10;     // Gravitational constant! Not to be confused with gravitational acceleration (g)

    private static double STABILITY_TIME = 0;
    private static double MAX_TIME = 0;
    public static double dt;

    public static int timeIdx = 0;
    public static double lastDeathTime = 0;
    public static double time = 0;

    public static int N = 50;

    private static double kn = 10e5;
    private static double kt = 2*kn;

    private static double slidingWindow = 0.8;     // Probando a ojimetro
//    private static List<Double> collidedParticles = new ArrayList<>();
    private static double collidedParticles = 0;

    private static NumberFormat defaultFormat = NumberFormat.getPercentInstance();

    private static List<Particle> particles = new ArrayList<>();
    public static Particle sun;
    public static double o = 0.55;

    private static final int FRAME_COUNT = 100;
    private static int OVITO_DT = (int) (1 / dt);
    private static final double WIDTH = 300;
    private static final double HEIGHT = 300;

    private static boolean exportFrames = false, exportData = false;
    private static String fileOut = null;

    private static Map<DataType, List<Double>> data = new HashMap<>();
    private static final double forceLoss = 0.5;

    private static List<Integer> exportIdxs = new ArrayList<>();
    private static int EXPORT_DT = (int) (15 / dt);

    public enum Orientation {
        CLOCK,
        COUNTER
    }

    public enum DataType {
        TIME("t", "0.###"),
        PARTICLE_COUNT("n", "0"),
        CLOCKWISE_PARTICLES("clock", "0"),
        COLLISIONS("collisions", "0");

        private String name;
        private DecimalFormat fmt;
        DataType(String name, String fmt) {
            this.name = name;
            this.fmt = new DecimalFormat(fmt);
        }

        public DecimalFormat getFmt() {
            return fmt;
        }

        @Override
        public String toString() {
            return name;
        }
    }


    public static void main(String[] args) throws Exception {
        parseArguments(args);

        for (int i = 2; i <= 2; i++) {
            initializeDataArrays();
            timeIdx = 0;
            simulate();

            DataExporter dataExporter = new DataExporter();
            if (fileOut == null)
             fileOut = dt + "," + N + ".data";

            System.out.println("Saving data file to: " + fileOut);
            dataExporter.export(data, fileOut);
        }
    }

    private static void initializeDataArrays() {
        for (DataType dataType : getVariableDataHeader())
            data.put(dataType, new ArrayList<>());
    }

    private static boolean isDone() {
        if (MAX_TIME == 0) {
            return time - lastDeathTime >= STABILITY_TIME;
        }
        return time >= MAX_TIME;
    }

    private static void simulate() throws IOException {
        Ovito ovito = new Ovito();
        List<Particle> ovitoParticles = new ArrayList<>();
        int i = 0;
        sun = Particle.builder().mass(SUN_MASS).radius(SUN_RADIUS).x(0).y(0).id(0).build();
        initializeParticles(N, SPAWN_DISTANCE, o);
        System.out.println(particles.size());

        time = 0;
        lastDeathTime = 0;

        while (!isDone()) {
            particles.forEach(Particle::clearForces);
            computeCollisions();
            computeGravityForces();
            updateParticlePositions();
            time += dt;

            Collection<Particle> deadParticles = overlappingSun(particles);
            if (deadParticles.size() > 0) {
                lastDeathTime = time;
            }
            particles.removeAll(deadParticles);

            time += dt;

            if (exportFrames && i % OVITO_DT == 0) {
                ovitoParticles.clear();
                ovitoParticles.add(sun);
                ovitoParticles.addAll(particles);
                ovito.createFile(i/OVITO_DT, ovitoParticles, WIDTH, HEIGHT);
            }


            if (i % EXPORT_DT == 0) {
                if (exportData) {
                    data.get(DataType.TIME).add(time);
                    data.get(DataType.PARTICLE_COUNT).add((double) particles.size());
                    data.get(DataType.CLOCKWISE_PARTICLES).add((double) getParticlesMoving(Orientation.CLOCK));
                    data.get(DataType.COLLISIONS).add(collidedParticles);
                    exportIdxs.add(timeIdx);
                }

                collidedParticles = 0;
            }
            i++;
            timeIdx++;
        }
        System.out.println(particles.size());


    }

    private static Long getParticlesMoving(Orientation orientation) {
        return particles.stream().filter(p -> p.getOrientation() == orientation).count();
    }

    private static Collection<Particle> overlappingSun(Collection<Particle> particles) {
        return particles.stream().filter(p -> p.isOverlapping(sun)).collect(Collectors.toCollection(ArrayList::new));
    }

    private static void computeCollisions() {
        for(int i = 0; i < particles.size(); i++) {
            Particle pi = particles.get(i);
            for (int j = i + 1; j < particles.size(); j++) {
                Particle pj = particles.get(j);
                if (pi.isOverlapping(pj)) {
                    pi.setColliding(true);
                    pj.setColliding(true);
                    computeCollision(pi, pj);
                    collidedParticles++;
                }
            }
        }
    }


    private static void computeCollision(Particle pi, Particle pj){
        ArrayList<Double> v_rel = pi.v_rel(pj);
        double enx = pi.enx(pj);
        double eny = pi.eny(pj);

        double etx = -eny;
        double ety = enx;

        double overlap = pi.overlap(pj);

        double Fn = -kn * overlap;
        double Ft = -kt * overlap * (v_rel.get(0) * etx + v_rel.get(1) * ety);
        double Fx = Fn * enx + Ft * etx;
        double Fy = Fn * eny + Ft * ety;
        Fx *= forceLoss;
        Fy *= forceLoss;

        pi.addNormalForce(Fn);
        pi.addTangentialForce(Math.abs(Ft));
        pj.addTangentialForce(Math.abs(Ft));
        pi.applyForce(Fx, Fy);
        pj.applyForce(-Fx, -Fy);
    }

    private static void updateParticlePositions(){
        particles.parallelStream().forEach(p -> p.move(dt));
    }

    private static void computeGravityForces() {
        particles.forEach(Orbit::computeGravityForces);
    }

    private static void computeGravityForces(Particle p) {
        if(p.isColliding())
            return;

        double theta = FastAtan.atan2(p.getY(), p.getX());
        double fMag = G * p.getMass() * sun.getMass() / p.distance_sq(sun);
        double fx = - fMag * Math.cos(theta);
        double fy = - fMag * Math.sin(theta);

        p.applyForce(fx, fy);
    }

    private static int initializeParticles(int N, double distance, double orientationProportion) {
        Particle p;
        particles.clear();

        int n = 0;
        int id = 1;
        particles = new ArrayList<>();

        int positiveCount = (int) Math.floor(orientationProportion * N);
        for (; n < N; n++) {
            do {
                double theta = Math.random() * 2*Math.PI;
                double x = distance * Math.cos(theta);
                double y = distance * Math.sin(theta);

                long positive = particles.stream().filter(particle -> particle.getOrientation() == Orientation.CLOCK).count();
                double sign = positive >= positiveCount ? 1 : -1;

                double vt = VT0 * Math.signum(sign);

                double vn = (Math.random() - 0.5) * VN0;

                double vx = - vt * Math.sin(theta) + vn * Math.cos(theta);
                double vy = vt * Math.cos(theta) + vn * Math.sin(theta);
                p = Particle.builder().x(x).y(y).mass(PARTICLE_MASS).radius(PARTICLE_RADIUS).vx(vx).vy(vy).fn(0).ft(0).id(id).build();
            } while (!p.isValid(particles));
            id++;
            particles.add(p);
        }

        return n;
    }

    private static void parseArguments(String[] args) throws ParseException {
        CommandLine commandLine;
        Option option_dt = Option.builder("dt")
                .required(false)
                .desc("dt")
                .hasArg()
                .build();
        Option export_frames = Option.builder("frames")
                .required(false)
                .desc("export_frames")
                .build();
        Option export_data = Option.builder("data")
                .required(false)
                .desc("export_data")
                .build();
        Option max_time = Option.builder("T")
                .required(false)
                .desc("max_time")
                .hasArg()
                .build();
        Option stability_time = Option.builder("ST")
                .required(false)
                .desc("stability_time")
                .hasArg()
                .build();
        Option orientation = Option.builder("o")
                .required(false)
                .desc("orientation_rate")
                .hasArg()
                .build();
        Option output = Option.builder("out")
                .required(false)
                .desc("output")
                .hasArg()
                .build();
        Option n = Option.builder("n")
                .required(false)
                .desc("n")
                .hasArg()
                .build();

        Options options = new Options();
        options.addOption(option_dt);
        options.addOption(export_frames);
        options.addOption(export_data);
        options.addOption(max_time);
        options.addOption(stability_time);
        options.addOption(orientation);
        options.addOption(output);
        options.addOption(n);
        CommandLineParser parser = new DefaultParser();

        commandLine = parser.parse(options, args);

        if(commandLine.hasOption("dt")) {
            dt = Double.parseDouble(commandLine.getOptionValue("dt"));
            OVITO_DT = (int) (1 / dt);
            EXPORT_DT = (int) (1 / dt);
        }
        if(commandLine.hasOption("ST")) {
            STABILITY_TIME = Double.parseDouble(commandLine.getOptionValue("ST"));
        }
        if(commandLine.hasOption("T")) {
            MAX_TIME = Double.parseDouble(commandLine.getOptionValue("T"));
        }
        if (commandLine.hasOption("o")) {
            o = Double.parseDouble(commandLine.getOptionValue("o"));
        }
        if (commandLine.hasOption("out")) {
            fileOut = commandLine.getOptionValue("out");
        }
        if (commandLine.hasOption("n")) {
            N = Integer.parseInt(commandLine.getOptionValue("n"));
        }

        exportFrames  = commandLine.hasOption("frames");
        exportData = commandLine.hasOption("data");
    }

    public static List<Orbit.DataType> getVariableDataHeader() {
        return Arrays.asList(
                DataType.TIME,
                DataType.PARTICLE_COUNT,
                DataType.CLOCKWISE_PARTICLES,
                DataType.COLLISIONS);
    }

}
