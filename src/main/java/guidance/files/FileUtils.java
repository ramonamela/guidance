package guidance.files;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import guidance.exceptions.GuidanceTaskException;
import guidance.processes.ProcessUtils;
import guidance.utils.ParseCmdLine;

public class FileUtils {

	// Logger
	private static final Logger LOGGER = LogManager.getLogger("Console");

	// Messages
	private static final String GUIDANCE_MSG_HEADER = "[Guidance]";
	private static final String GUIDANCE_MSG_CREATE_FILE = GUIDANCE_MSG_HEADER + " Creating ";

	// TODO: this constants are associated to a function that should not be here
	private static final String SAMTOOLSBINARY = "SAMTOOLSBINARY";

	private static final String STDOUT_EXTENSION = ".stdout";
	private static final String STDERR_EXTENSION = ".stderr";

	// File paths
	private static final String COMMON_REL_PATH = "common";
	private static final String MIXED_REL_PATH = "mixed";
	private static final String CHR_REL_PATH = "Chr_";
	private static final String ASSOCIATIONS_REL_PATH = "associations";

	private FileUtils() {
		// Private constructor to avoid instantiation
	}

	/**
	 * Method that creates all the output directory structure for the results.
	 * 
	 * @param parsingArgs
	 * @param myOutDir
	 * @param refPanels
	 * @param startChr
	 * @param endChr
	 */
	public static void createDirStructure(ParseCmdLine parsingArgs, String myOutDir, List<String> refPanels,
			int startChr, int endChr) {
		String mixedCohort = null;
		// We create the first directory: the cohort directory.
		mixedCohort = parsingArgs.getCohort();

		String tmpOutDir = "";
		tmpOutDir = myOutDir + File.separator + mixedCohort;

		LOGGER.info("[Guidance] Creating " + tmpOutDir);
		FileUtils.createDir(tmpOutDir);

		// We create the second directory: the REFPANEL directory.
		// It can be => common: for files from the beginning to phasing
		// => PANEL : for files for imputation for each panel

		// Now I create the directories for common and combined
		String tmpPanelDir = tmpOutDir + File.separator + COMMON_REL_PATH;
		FileUtils.createDir(tmpPanelDir);

		// Next level: Create mixed directories.
		String mixedOutDir = tmpPanelDir + File.separator + MIXED_REL_PATH;
		LOGGER.info(GUIDANCE_MSG_CREATE_FILE + mixedOutDir);
		FileUtils.createDir(mixedOutDir);
		for (int i = startChr; i <= endChr; i++) {
			String tmpChrDir = mixedOutDir + File.separator + CHR_REL_PATH + i;
			FileUtils.createDir(tmpChrDir);
		}

		// Now I create the directories for imputeOutDir
		for (int j = 0; j < refPanels.size(); j++) {
			String rPanel = refPanels.get(j);
			tmpPanelDir = tmpOutDir + File.separator + rPanel;
			FileUtils.createDir(tmpPanelDir);

			// Next level: Create mixed directories.
			mixedOutDir = tmpPanelDir + File.separator + MIXED_REL_PATH;
			LOGGER.info(GUIDANCE_MSG_CREATE_FILE + mixedOutDir);
			FileUtils.createDir(mixedOutDir);
			for (int i = startChr; i <= endChr; i++) {
				String tmpChrDir = mixedOutDir + File.separator + CHR_REL_PATH + i;
				FileUtils.createDir(tmpChrDir);
			}
		}

		// Then we create the directories for the Association files
		// This directory is created by default in myOutDir/association.
		String tmpMyOutDir1 = myOutDir + File.separator + ASSOCIATIONS_REL_PATH;
		FileUtils.createDir(tmpMyOutDir1);
		String testTypeName = null;
		int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
		for (int kk = 0; kk < numberOfTestTypes; kk++) {
			testTypeName = parsingArgs.getTestTypeName(kk);
			String tmpMyOutDir = myOutDir + File.separator + ASSOCIATIONS_REL_PATH + File.separator + testTypeName;
			FileUtils.createDir(tmpMyOutDir);

			for (int j = 0; j < refPanels.size(); j++) {
				String rPanel = refPanels.get(j);
				String assocDir = null;
				assocDir = tmpMyOutDir + File.separator + mixedCohort + "_for_" + rPanel;
				LOGGER.info(GUIDANCE_MSG_CREATE_FILE + assocDir);
				FileUtils.createDir(assocDir);

				for (int i = startChr; i <= endChr; i++) {
					String tmpChrDir = assocDir + "/Chr_" + i;
					FileUtils.createDir(tmpChrDir);
				}

				String summaryDir = assocDir + File.separator + "summary";
				LOGGER.info(GUIDANCE_MSG_CREATE_FILE + summaryDir);
				FileUtils.createDir(summaryDir);
			}

			// Finally, if there are several panel and if combinedRefPanel = 1
			// then we have to create the output directory for the combined ref_panel
			boolean refPanelCombine = parsingArgs.getRefPanelCombine();

			if (refPanelCombine) {
				// We take the first refPanel name
				String rPanel = refPanels.get(0);
				String combinedRefPanel = rPanel;
				for (int j = 1; j < refPanels.size(); j++) {
					rPanel = refPanels.get(j);
					combinedRefPanel = combinedRefPanel + "_" + rPanel;
				}

				String combinedAssocOutDir = null;
				combinedAssocOutDir = myOutDir + File.separator + ASSOCIATIONS_REL_PATH + File.separator + testTypeName
						+ File.separator + mixedCohort + "_combined_panels_" + combinedRefPanel;

				LOGGER.info(GUIDANCE_MSG_CREATE_FILE + combinedAssocOutDir);
				FileUtils.createDir(combinedAssocOutDir);

			} // End if
		} // End for kk < numberOfTestTyes

		// Now create the structure for the phenotype analysis results
		String phenomeAnalysisOutDir = myOutDir + File.separator + ASSOCIATIONS_REL_PATH + File.separator
				+ "pheno_analysis" + File.separator;
		LOGGER.info(GUIDANCE_MSG_CREATE_FILE + phenomeAnalysisOutDir);
		FileUtils.createDir(phenomeAnalysisOutDir);

		phenomeAnalysisOutDir = phenomeAnalysisOutDir + mixedCohort;
		// for(int tt=0; tt< numberOfTestTypes; tt++) {
		// testTypeName = parsingArgs.getTestTypeName(tt);
		// phenomeAnalysisOutDir = phenomeAnalysisOutDir + "_" + testTypeName;
		// }
		LOGGER.info(GUIDANCE_MSG_CREATE_FILE + phenomeAnalysisOutDir);
		FileUtils.createDir(phenomeAnalysisOutDir);
	}

	/**
	 * Method that creates a single directory
	 * 
	 * @param outputDir
	 */
	private static void createDir(String outputDir) {
		File tmpOutDir = new File(outputDir);
		boolean existOutDir = tmpOutDir.exists();
		if (!existOutDir) {
			boolean successExistOutDir = (new File(outputDir)).mkdir();
			if (!successExistOutDir) {
				LOGGER.fatal("[createDir] Error, cannot create " + tmpOutDir + " directory");
				System.exit(1);
			}
		}
	}

	/**
	 * Copy a file from @source to @dest. Returns true if success, false otherwise
	 * 
	 * @param source
	 * @param dest
	 * @throws IOException 
	 */
	public static void copy(String source, String dest) throws IOException {
		File sourceFile = new File(source);
		File destFile = new File(dest);
		Files.copy(sourceFile.toPath(), destFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
	}

	/**
	 * Moves a file or a directory from @source to @{dest}. Returns true if success,
	 * false otherwise
	 * 
	 * @param source
	 * @param dest
	 * @return
	 */
	public static boolean move(String source, String dest) {
		File srcFile = new File(source);
		File destFile = new File(dest);

		// Rename file (or directory)
		return srcFile.renameTo(destFile);
	}

	/**
	 * Deletes the file denoted by the filename @file
	 * 
	 * @param file
	 * @throws GuidanceTaskException
	 */
	public static void delete(String file) throws GuidanceTaskException {
		boolean success = new File(file).delete();
		if (!success) {
			throw new GuidanceTaskException("ERROR: Cannot erase file " + file);
		}
	}

	/**
	 * Creates an empty file if it does not exist
	 * 
	 * @param fileName
	 * @param moduleName
	 * @return true if the file was created or already existent, false otherwise
	 * @throws IOException
	 */
	public static boolean createEmptyFile(String fileName, String moduleName) throws IOException {
		File fa = new File(fileName);
		if (!fa.exists()) {
			System.out.println(moduleName + " The file " + fileName + " does not exist, then we create an empty file");
			return fa.createNewFile();
		}

		// The file already exists
		return true;
	}

	/**
	 * Method to zip a file
	 * 
	 * @param sourceFilePath
	 * @param destZipFilePath
	 */
	public static void gzipFile(String sourceFilePath, String destZipFilePath) {
		byte[] buffer = new byte[1024];
		File fInput = new File(sourceFilePath);
		String path = fInput.getParent();
		Path dir = Paths.get(path);

		System.out.println("Files in " + path + " :");
		try (DirectoryStream<Path> stream = Files.newDirectoryStream(dir, "*")) {
			for (Path file : stream) {
				System.out.println(file.toAbsolutePath().toString());
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("End files in sandbox");

		try (FileOutputStream fileOutputStream = new FileOutputStream(destZipFilePath);
				GZIPOutputStream gzipOuputStream = new GZIPOutputStream(fileOutputStream);
				FileInputStream fileInput = new FileInputStream(sourceFilePath)) {

			int bytesRead;
			while ((bytesRead = fileInput.read(buffer)) > 0) {
				gzipOuputStream.write(buffer, 0, bytesRead);
			}

			gzipOuputStream.finish();
		} catch (IOException ioe) {
			System.err.println("ERROR: Cannot zip file");
			/*
			 * try { Thread.sleep(600000); } catch (InterruptedException e) { // TODO
			 * Auto-generated catch block e.printStackTrace(); }
			 */
			ioe.printStackTrace();
		}
	}

	// TODO: remove when removing bgzipFile
	private static String loadFromEnvironment(String envVarName, String methodHeader) throws GuidanceTaskException {
		String envVar = System.getenv(envVarName);
		if (envVar == null) {
			throw new GuidanceTaskException(
					methodHeader + "Error, " + envVarName + " environment variable is not defined in .bashrc!!!");
		}

		return envVar;
	}

	// TODO: this function should not be here, just a "fast fix"
	public static void bgzipFile(String input, String output) throws GuidanceTaskException {

		String samToolsBinary = FileUtils.loadFromEnvironment(SAMTOOLSBINARY, "[samtoolsBgzipFile]");

		long startTime = System.currentTimeMillis();

		String cmd = samToolsBinary + "/bgzip -f " + input + " " + output;

		LOGGER.debug("\n[DEBUG] Command: " + cmd);
		LOGGER.debug(" ");

		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, input + STDOUT_EXTENSION, input + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		if (exitValue != 0) {
			System.err.println("[phasing] Warning executing samtoolsBgzip job, exit value is: " + exitValue);
			System.err.println("[phasing]                         (This warning is not fatal).");
		}

		File tmpOutput = new File(input + ".gz");
		if (tmpOutput.exists()) {
			boolean success = FileUtils.move(input + ".gz", output);
			if (!success) {
				throw new GuidanceTaskException(
						"[samtoolsBgzip] Error, the file " + input + ".gz was not succesfully renamed to " + output);
				// File was not successfully renamed
			}
		} else {
			throw new GuidanceTaskException(
					"[samtoolsBgzip] Error, the file " + input + ".gz was not succesfully created");
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
	}

	/**
	 * Method to uncompress a file
	 * 
	 * @param compressedFile
	 * @param decompressedFile
	 */
	public static void gunzipFile(String compressedFile, String decompressedFile) {
		byte[] buffer = new byte[1024];

		try (FileInputStream fileIn = new FileInputStream(compressedFile);
				GZIPInputStream gZIPInputStream = new GZIPInputStream(fileIn);
				FileOutputStream fileOutputStream = new FileOutputStream(decompressedFile)) {

			int bytesRead;
			while ((bytesRead = gZIPInputStream.read(buffer)) > 0) {
				fileOutputStream.write(buffer, 0, bytesRead);
			}
		} catch (IOException ioe) {
			System.err.println("ERROR: Cannot zip file");
			ioe.printStackTrace();
		}
	}
	
	public static void recursiveSearch(String filePath) {
		System.out.println("Looking into the folder of " + filePath);
		File file = new File(filePath);
		file = new File(file.getParent());
		System.out.println("The folder is " + file.getParent());
		File[] filesList = file.listFiles();
		for (File f : filesList) {
			if (f.isFile()) {
				System.out.println("File name is  -------------->" + f.getName());
			}
		}
	}
	
	public static void getFile(String filename) {
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filename)))){
			String line = reader.readLine();
		} catch (FileNotFoundException e) {
			System.err.println("[DEBUG] File " + filename + " does not exist");
		} catch (IOException e) {
			System.err.println("[DEBUG] Error when bringing back " + filename);
		}
	}

}
