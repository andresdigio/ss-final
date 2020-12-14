package model;

import lombok.Builder;
import lombok.Data;
import simulation.Orbit;

import java.util.ArrayList;
import java.util.List;
import java.util.function.ToDoubleBiFunction;

@Data
@Builder(toBuilder =  true)
public class Particle {
    private final long id;
    private final double mass;
    private final double radius;
    private final double potentialRadius;
    private double x, y;
    private Particle prevParticle;
    private double vx, vy;
    private double ax, ay;
    private double fn;
    private double ft;
    private ToDoubleBiFunction<Double, Double> force;
    private boolean isColliding = false;

    public void applyForce(double fx, double fy) {
        ax += fx/mass;
        ay += fy/mass;
    }

    public void clearForces() {
        ax = 0;
        ay = 0;
        isColliding = false;
    }

    public boolean isValid(List<Particle> particles) {
        for(Particle p : particles){
            if(p != this) {
                double dist_sq = distance_sq(p);
                if (dist_sq <= (p.radius + radius) * (p.radius + radius))
                    return false;
            }
        }

        return true;
    }

    public double distance_sq(Particle other) {
        double dx = other.x - x;
        double dy = other.y - y;

        return dx*dx + dy*dy;
    }

    public double calculateAngle(Particle other) {
        double dx = other.x - x;
        double dy = other.y - y;

        return Math.atan2(dy, dx);
    }

    public double getVelocity() {
        return Math.sqrt(vx*vx + vy*vy);
    }

    public double getVelocity(double dt) {
        if (prevParticle != null) {
            vx = (x - prevParticle.getX()) / dt;
            vy = (y - prevParticle.getY()) / dt;
        }

        return Math.sqrt(vx*vx + vy*vy);
    }

    public void move(double dt) {
        if (prevParticle == null) {
            prevParticle = Particle.builder().x(x-vx*dt).y(y-vy*dt).build();
        }

        double dt2  = dt*dt;
        double nextX = 2*x - prevParticle.x + ax * dt2;
        double nextY = 2*y - prevParticle.y + ay * dt2;

        prevParticle = Particle.builder().x(x).y(y).build();
        x = nextX;
        y = nextY;
        vx = (x - prevParticle.x) / dt;
        vy = (y - prevParticle.y) / dt;
    }

    public void checkPrevParticle(double dt) {
        if (prevParticle == null) {
            prevParticle = Particle.builder().x(x - vx * dt).y(y - vy * dt).build();
        }
    }

    public Orbit.Orientation orientation() {
        if (x > 0) {
            return vy > 0 ? Orbit.Orientation.COUNTER : Orbit.Orientation.CLOCK;
        } else if (x < 0) {
            return vy < 0 ? Orbit.Orientation.COUNTER : Orbit.Orientation.CLOCK;
        } else {
            return vx > 0 ? Orbit.Orientation.CLOCK : Orbit.Orientation.COUNTER;
        }
    }

    private double getFx(){
        return mass * ax;
    }

    private double getFy(){
        return mass * ay;
    }

    public boolean isOverlapping(Particle p){
        return overlap(p) > 0;
    }

    public double overlap(Particle p){
        return radius + p.getRadius() - getVectorDistance(p);
    }

    public double getVectorDistance(Particle p){
        return Math.sqrt(distance_sq(p));
    }

    public double enx(Particle p){
        return (p.getX() - x) / getVectorDistance(p);
    }

    public double eny(Particle p){
        return (p.getY() - y) / getVectorDistance(p);
    }

    public ArrayList<Double> v_rel(Particle p){
        ArrayList<Double> v_rel = new ArrayList<>();
        v_rel.add(vx - p.vx);
        v_rel.add(vy - p.vy);
        return v_rel;
    }

    public double kineticEnergy() {
        return 0.5*mass*Math.pow(getVelocity(),2);
    }

    public double potentialEnergy(Particle sun) {
        return Orbit.GRAVITY * getVectorDistance(sun) * mass;
    }

    public double computeEnergy(Particle sun) {
        return potentialEnergy(sun) + kineticEnergy();
    }

    public double velocity(){
        return Math.sqrt(vx*vx + vy*vy);
    }

    public void addNormalForce(double Fn) { this.fn += Fn; }

    public void addTangentialForce(double Ft) { this.ft += Ft; }
}
