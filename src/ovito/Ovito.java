package ovito;

import model.Particle;
import simulation.Orbit;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Locale;

public class Ovito {
    private static final Path results = Paths.get("results");

    public Ovito() throws IOException {
        try {
            Files.createDirectory(results);
        } catch (FileAlreadyExistsException ignored) {}
    }

    public void createFile(int i, Collection<Particle> particles, double MAX_X, double MAX_Y){
        try{
            Files.deleteIfExists(results.resolve(getFileName(i)));
            FileWriter file = new FileWriter(results.resolve(getFileName(i)).toString(), true);
            if(!setFileHeader(file, particles.size(), MAX_X, MAX_Y)) {
                return;
            }

            for(Particle p : particles) {
                if(p.getId() == 0){
                    file.write(p.getId() + " " + p.getX() + " " + p.getY() + " " + p.getRadius() + " " + "0" + "\n");
                }
                else {
                    file.write(p.getId() + " " + p.getX() + " " + p.getY() + " " + p.getRadius() + " " + (p.getOrientation() == Orbit.Orientation.CLOCK ? "1" : "-1") + "\n");
                }
            }
            file.close();
        }
        catch (IOException ignored) {
        }
    }

    private boolean setFileHeader(FileWriter file, int particleCount, double MAX_X, double MAX_Y){
        try {
            file.write(particleCount+2 + "\n");
            file.write("\n");
            file.write("-2 0.0 0.0 0.0 0\n");
            file.write(String.format(Locale.ENGLISH, "-1 %f %f 0.0 0\n", MAX_X, MAX_Y));
        }
        catch (IOException e){
            return false;
        }
        return true;
    }

    private String getFileName(int i) {
        return "ovito-" + i + ".xyz";
    }
}
