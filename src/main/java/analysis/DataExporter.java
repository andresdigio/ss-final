package analysis;

import simulation.Orbit;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

public class DataExporter {
    private static final Path cwd = Paths.get(System.getProperty("user.dir"));
    private static final Path projectDir = cwd.getParent().getParent().getParent();
    private static final Path results = projectDir.resolve("results");
    private static final Path dataDir = projectDir.resolve("data");

    public DataExporter() throws IOException {
        try {
            Files.createDirectory(results);
            Files.createDirectory(dataDir);
        } catch (FileAlreadyExistsException ignored) {

        }
    }

    public void export(Map<Orbit.DataType, List<Double>> data, String fileName) throws IOException {
        FileWriter file = createFile(fileName);
        int rows = data.get(Orbit.DataType.TIME).size();    // Amount of time steps

        for (int i = 0; i < rows; i++) {
            StringBuilder sb = new StringBuilder();

            sb.append(Orbit.dt).append(",").append(Orbit.N).append(",").append(Orbit.o).append(",");

            List<Orbit.DataType> dataHeader = Orbit.getVariableDataHeader();

            for (Orbit.DataType dataType : dataHeader) {
                Double columnValue = data.get(dataType).get(i);
                String formattedValue = dataType.getFmt().format(columnValue);

                sb.append(formattedValue).append(",");
            }

            sb.deleteCharAt(sb.lastIndexOf(","));
            sb.append("\n");
            file.write(sb.toString());
        }

        file.close();
    }

    public FileWriter createFile(String fileName) throws IOException {
        Files.deleteIfExists(dataDir.resolve(fileName));
        FileWriter file = new FileWriter(dataDir.resolve(fileName).toString(), true);
        file.write(getFileHeader() + "\n");

        return file;
    }

    private String getFileHeader() {
        List<Orbit.DataType> dataHeader = Orbit.getVariableDataHeader();

        return "dt,N,o," + dataHeader.stream().map(Objects::toString).collect(Collectors.joining(","));
    }
}
