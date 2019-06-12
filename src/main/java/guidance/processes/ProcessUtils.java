package guidance.processes;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;


public class ProcessUtils {

    /**
     * Executes the given cmd and retrieves the STDOUT and STDERR to the given files. At the end of the command
     * execution returns the exitValue
     * 
     * @param cmd
     * @param outputFile
     * @param errorFile
     * @return
     * @throws IOException
     */
    public static int execute(String cmd, String outputFile, String errorFile) throws IOException {
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
    
    public static int execute(String cmd, String outputFile, String errorFile, String variableToErase) throws IOException {
        // Create the process
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));

        // Remove unnecessary environment
        pb.environment().remove("LD_PRELOAD");
        pb.environment().remove(variableToErase);

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
    
    /**
     * Executes the given cmd. At the end of the command execution returns the exitValue
     * 
     * @param cmd
     * @return
     * @throws IOException
     */
    public static int executeWithoutOutputsSNP(String cmd) throws IOException {

        //String[] commandArray = { "bash", "-c", cmd };

        // Create the process
        //ProcessBuilder pb = new ProcessBuilder(commandArray);
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));

        // Remove unnecessary environment
        pb.environment().remove("LD_PRELOAD");

        // Start the process
        Process p = pb.start();

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

    /**
     * Executes the given cmd. At the end of the command execution returns the exitValue
     * 
     * @param cmd
     * @return
     * @throws IOException
     */
    public static int executeWithoutOutputs(String cmd) throws IOException {

        String[] commandArray = { "bash", "-c", cmd };

        // Create the process
        ProcessBuilder pb = new ProcessBuilder(commandArray);
        //ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));

        // Remove unnecessary environment
        pb.environment().remove("LD_PRELOAD");

        // Start the process
        Process p = pb.start();

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

    /**
     * Execute a bash command. We can handle complex bash commands including multiple executions (; | && ||), quotes,
     * expansions ($), escapes (\), e.g.: "cd /abc/def; mv ghi 'older ghi '$(whoami)"
     * 
     * @param command
     * @return true if bash got started, but your command may have failed.
     */
    public static int executeBashCommand(String command, String outputFile, String errorFile) throws IOException {

        Runtime r = Runtime.getRuntime();
        // Use bash -c so we can handle things like multi commands separated by ; and
        // things like quotes, $, |, and \. My tests show that command comes as
        // one argument to bash, so we do not need to quote it to make it one thing.
        // Also, exec may object if it does not have an executable file as the first thing,
        // so having bash here makes it happy provided bash is installed and in path.
        String[] commands = { "bash", "-c", command };

        Process p = r.exec(commands);

        int exitValue = -1;
        try {
            exitValue = p.waitFor();
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
        }

        readOutputAndError(p.getInputStream(), outputFile, p.getErrorStream(), errorFile);

        return exitValue;

    }

    /**
     * Execute a bash command. We can handle complex bash commands including multiple executions (; | && ||), quotes,
     * expansions ($), escapes (\), e.g.: "cd /abc/def; mv ghi 'older ghi '$(whoami)"
     * 
     * @param command
     * @return true if bash got started, but your command may have failed.
     */
    public static int executeBashCommandWithoutOutput(String command) throws IOException, InterruptedException {

        Runtime r = Runtime.getRuntime();
        // Use bash -c so we can handle things like multi commands separated by ; and
        // things like quotes, $, |, and \. My tests show that command comes as
        // one argument to bash, so we do not need to quote it to make it one thing.
        // Also, exec may object if it does not have an executable file as the first thing,
        // so having bash here makes it happy provided bash is installed and in path.
        String[] commands = { "bash", "-c", command };

        System.out.println("[DEBUG] Executed command: " + command);

        Process p = r.exec(commands);

        int exitValue = -1;
        try {
            exitValue = p.waitFor();
        } catch (InterruptedException ie) {
            Thread.currentThread().interrupt();
        }

        return exitValue;

    }

    /**
     * Processes the input streams to their given output and error files
     * 
     * @param in
     * @param outputPath
     * @param error
     * @param errorPath
     * @throws IOException
     */
    private static void readOutputAndError(InputStream in, String outputPath, InputStream error, String errorPath) throws IOException {
    	
        try (BufferedInputStream bisInp = new BufferedInputStream(in);
             BufferedOutputStream bosInp = new BufferedOutputStream(new FileOutputStream(outputPath));
				BufferedInputStream bisErr = new BufferedInputStream(error);
				BufferedOutputStream bosErr = new BufferedOutputStream(new FileOutputStream(errorPath))
        		) {

            byte[] bOut = new byte[1024];
            byte[] bErr = new byte[1024];
            int readOut = -1;
            int readErr = -1;
            while ((readOut = bisInp.read(bOut)) >= 0 || (readErr = bisErr.read(bErr)) >= 0) {
                if(readOut >= 0) {
                	bosInp.write(bOut, 0, readOut);
                	readOut = -1;
                }
                if(readErr >= 0) {
                	bosErr.write(bErr, 0, readErr);
                	readErr = -1;
                }
            }
        }
    }

}
