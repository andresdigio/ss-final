package simulation;

import model.FastAtan;
import model.Particle;
import ovito.Ovito;

import java.text.NumberFormat;

import java.util.*;
import java.util.stream.Collectors;

public class Orbit {
    private static final double SUN_MASS = 5000;
    private static final double SUN_RADIUS = 50;
    private static final double PARTICLE_MASS = 1;
    private static final double PARTICLE_RADIUS = 5;
    private static final double SPAWN_DISTANCE = 300;
    private static final double VT0 = 13;
    private static final double VN0 = 2;
    public static final double GRAVITY = 9.8;

    private static final double MAX_TIME = 300;

    private static final double dt = 1e-5;
    private static final int N = 50;

    private static double kn = 10e5;
    private static double kt = 2*kn;

    private static List<Double> energyList;
    private static NumberFormat defaultFormat = NumberFormat.getPercentInstance();

    private static List<Particle> particles = new ArrayList<>();
    private static Particle sun;

    private static final int FRAME_COUNT = 500;
    private static final int OVITO_DT = (int) (MAX_TIME / (dt * FRAME_COUNT));
    private static final double WIDTH = 300;
    private static final double HEIGHT = 300;

    public enum Orientation {
        CLOCK,
        COUNTER
    }

    public static void main(String[] args) {
        energyList = new ArrayList<>();
        simulate();
    }

    private static void simulate() {
        Ovito ovito = null;
        List<Particle> ovitoParticles = new ArrayList<>();
        int i=0;
        try {
            ovito = new Ovito();
        }
        catch (Exception e){
            System.out.println(e);
        }
        int particleCount = initializeParticles(N, SPAWN_DISTANCE, 0.8);
        sun = Particle.builder().mass(SUN_MASS).radius(SUN_RADIUS).x(0).y(0).id(0).build();

        double time = 0;

        while(time < MAX_TIME) {
            energyList.add(getSystemEnergy(particles));
            particles.forEach(Particle::clearForces);
            computeCollisions();
            computeGravityForces();
            particles.removeAll(overlappingSun(particles));

            updateParticlePositions();

            time += dt;

            if(i % OVITO_DT == 0 && ovito != null){
                System.out.println("Progress: " + defaultFormat.format(time / MAX_TIME));
                System.out.println(computeOrientationCount(particles));
                ovitoParticles.clear();
                ovitoParticles.add(sun);
                ovitoParticles.addAll(particles);
                ovito.createFile(i/OVITO_DT, ovitoParticles, WIDTH, HEIGHT);
            }
            i++;
        }

        System.out.println("Min system energy: " + energyList.stream().min(Double::compareTo).get());
        System.out.println("Max system energy: " + energyList.stream().max(Double::compareTo).get());

    }

    private static Map<Orientation, Long> computeOrientationCount(Collection<Particle> particles) {
        Map<Orientation, Long> map = new HashMap<>();
        map.put(Orientation.CLOCK, particles.stream().filter(p -> p.orientation() == Orientation.CLOCK).count());
        map.put(Orientation.COUNTER, particles.stream().filter(p -> p.orientation() == Orientation.COUNTER).count());
        return map;
    }

    private static double getSystemEnergy(Collection<Particle> particles) {
        return particles.stream().mapToDouble(p -> p.computeEnergy(sun)).sum();
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
        double fMag = GRAVITY * p.getMass() * sun.getMass() / p.distance_sq(sun);
        double fx = - fMag * Math.cos(theta);
        double fy = - fMag * Math.sin(theta);

        p.applyForce(fx, fy);
    }

    private static int initializeParticles(int N, double distance, double orientationProportion) {
        Particle p;

        long tiempoRazonableNanos = 2 * (long) 1e9;
        long t0 = System.nanoTime();
        long t = t0;
        int n = 0;
        int id = 1;
        while((t-t0) < tiempoRazonableNanos && n < N) {
            boolean add = true;
            do {
                double theta = Math.random() * 2*Math.PI;
                double x = distance * Math.cos(theta);
                double y = distance * Math.sin(theta);

                // Clockwise if < 0, counterclockwise if > 0
                double orientationRand = Math.random() - orientationProportion;

                double vt = VT0 * Math.signum(orientationRand);

                double vn = (Math.random() - 0.5) * VN0;

                double vx = - vt * Math.sin(theta) + vn * Math.cos(theta);
                double vy = vt * Math.cos(theta) + vn * Math.sin(theta);
                p = Particle.builder().id(n).x(x).y(y).mass(PARTICLE_MASS).radius(PARTICLE_RADIUS).vx(vx).vy(vy).fn(0).ft(0).id(id).build();
                t = System.nanoTime();
                if(t-t0 > tiempoRazonableNanos) {
                    add = false;
                    break;
                }
            } while (!p.isValid(particles));

            if(add) {
                id++;
                particles.add(p);
                n++;
            }
            t = System.nanoTime();
        }

        return n;
    }


}
