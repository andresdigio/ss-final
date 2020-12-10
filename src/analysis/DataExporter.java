package analysis;

import simulation.Orbit;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.text.NumberFormat;
import java.util.Objects;
import java.util.stream.Collectors;

public class DataExporter {
    private static final Path results = Paths.get("results");
    private static final Path dataDir = Paths.get("data");

    public DataExporter() throws IOException {
        try {
            Files.createDirectory(results);
            Files.createDirectory(dataDir);
        } catch (FileAlreadyExistsException ignored) {

        }
    }

    public void createDataFile(String fileName, String constantHeader, String constantsValues, List<Orbit.DataType> dataHeader, List<Double> time, Map<Orbit.DataType, List<Double>> data) throws IOException {
        DecimalFormat timePattern = new DecimalFormat("0.###");
        Files.deleteIfExists(dataDir.resolve(fileName));
        FileWriter file = new FileWriter(dataDir.resolve(fileName).toString(), true);
        file.write(constantHeader + ",t," + dataHeader.stream().map(Objects::toString).collect(Collectors.joining(",")) + "\n");
        for (int i = 0; i < time.size(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(constantsValues).append(",").append(timePattern.format(time.get(i))).append(",");
            for (Orbit.DataType type : dataHeader) {
                String datum = type.getFmt().format(data.get(type).get(i));
                sb.append(datum).append(",");
            }
            sb.deleteCharAt(sb.lastIndexOf(","));
            sb.append("\n");
            file.write(sb.toString());
        }
        file.close();
    }

}
