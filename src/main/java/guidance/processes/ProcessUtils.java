package guidance.processes;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;


public class ProcessUtils {

    public static int exec(String cmd, String outputFile, String errorFile) throws IOException {
        // Create the process
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));

        // Remove unnecessary environment
        pb.environment().remove("LD_PRELOAD");

        // Start the process
        Process p = pb.start();

        // Handle the output and the error
        readOutputAndError(p.getInputStream(), outputFile, p.getErrorStream(), errorFile);

        // Retrieve the exit value
        int exitValue = -1;
        try {
            exitValue = p.waitFor();
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
        }

        // Return the exit value
        return exitValue;
    }

    private static void readOutputAndError(InputStream in, String outputPath, InputStream error, String errorPath) throws IOException {
        try (BufferedInputStream bisInp = new BufferedInputStream(in);
                BufferedOutputStream bosInp = new BufferedOutputStream(new FileOutputStream(outputPath))) {

            byte[] b = new byte[1024];
            int read;
            while ((read = bisInp.read(b)) >= 0) {
                bosInp.write(b, 0, read);
            }
        }

        try (BufferedInputStream bisErr = new BufferedInputStream(error);
                BufferedOutputStream bosErr = new BufferedOutputStream(new FileOutputStream(errorPath))) {

            byte[] b = new byte[1024];
            int read;
            while ((read = bisErr.read(b)) >= 0) {
                bosErr.write(b, 0, read);
            }
        }
    }

}
