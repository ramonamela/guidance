package guidance.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.logging.log4j.Logger;

import guidance.exceptions.GuidanceEnvironmentException;

public class GeneralUtils {

	public static void flushCommands(File listOfStages, ArrayList<String> listOfCommands, boolean option)
			throws IOException {

		if (option) {
			try (BufferedWriter writer = new BufferedWriter(new FileWriter(listOfStages, true))) {
				while (!listOfCommands.isEmpty()) {
					String str = listOfCommands.get(0);
					listOfCommands.remove(0);
					writer.write(str);
					writer.newLine();
					writer.newLine();
				}

				// Close the file with the list of commands... writer.flush();
				writer.close();
			}
		}

	}
	
	/**
	 * Method to print the current status of each stage (0: unactive, 1:active) of
	 * Guidance workflow.
	 * 
	 * @param parsingArgs
	 */
	public static void printStagesStatus(Logger LOGGER, ParseCmdLine parsingArgs) {
		// Verify the status of each stage:
		LOGGER.info("[Guidance] Current Status of each stage of the whole workflow:");
		LOGGER.info("[Guidance] convertFromBedToBed      " + parsingArgs.getStageStatus("convertFromBedToBed"));
		LOGGER.info("[Guidance] createRsIdList           " + parsingArgs.getStageStatus("createRsIdList"));
		LOGGER.info("[Guidance] phasingBed               " + parsingArgs.getStageStatus("phasingBed"));
		LOGGER.info("[Guidance] phasing                  " + parsingArgs.getStageStatus("phasing"));
		LOGGER.info("[Guidance] createListOfExcludedSnps " + parsingArgs.getStageStatus("createListOfExcludedSnps"));
		LOGGER.info("[Guidance] filterHaplotypes         " + parsingArgs.getStageStatus("filterHaplotypes"));
		LOGGER.info("[Guidance] imputeWithImpute         " + parsingArgs.getStageStatus("imputeWithImpute"));
		LOGGER.info("[Guidance] imputeWithMinimac        " + parsingArgs.getStageStatus("imputeWithMinimac"));
		LOGGER.info("[Guidance] filterByInfo             " + parsingArgs.getStageStatus("filterByInfo"));
		LOGGER.info("[Guidance] qctoolS                  " + parsingArgs.getStageStatus("qctoolS"));
		LOGGER.info("[Guidance] snptest                  " + parsingArgs.getStageStatus("snptest"));
		LOGGER.info("[Guidance] collectSummary           " + parsingArgs.getStageStatus("collectSummary"));
		LOGGER.info("[Guidance] mergeTwoChunks           " + parsingArgs.getStageStatus("mergeTwoChunks"));
		LOGGER.info("[Guidance] filterByAll              " + parsingArgs.getStageStatus("filterByAll"));
		LOGGER.info("[Guidance] jointCondensedFiles      " + parsingArgs.getStageStatus("jointCondensedFiles"));
		LOGGER.info("[Guidance] jointFilteredByAllFiles  " + parsingArgs.getStageStatus("jointFilteredByAllFiles"));
		LOGGER.info("[Guidance] generateTopHits          " + parsingArgs.getStageStatus("generateTopHits"));
		LOGGER.info("[Guidance] generateQQManhattanPlots " + parsingArgs.getStageStatus("generateQQManhattanPlots"));
		LOGGER.info("[Guidance] combinePanelsComplex     " + parsingArgs.getStageStatus("combinePanelsComplex"));
		LOGGER.info("[Guidance] combGenerateManhattanTop " + parsingArgs.getStageStatus("combGenerateManhattanTop"));
		LOGGER.info("[Guidance] phenoAnalysis            " + parsingArgs.getStageStatus("phenoAnalysis"));
		LOGGER.info("[Guidance] tasku                    " + parsingArgs.getStageStatus("tasku"));
		LOGGER.info("[Guidance] taskv                    " + parsingArgs.getStageStatus("taskv"));
		LOGGER.info("[Guidance] taskw                    " + parsingArgs.getStageStatus("taskw"));
		LOGGER.info("[Guidance] taskx                    " + parsingArgs.getStageStatus("taskx"));
		LOGGER.info("[Guidance] tasky                    " + parsingArgs.getStageStatus("tasky"));
		LOGGER.info("[Guidance] taskz                    " + parsingArgs.getStageStatus("taskz"));
	}
	
	public static File createListOfStages(Logger LOGGER, String listOfStagesFileName) throws IOException {
		// Verify whether the file exists or not.
		File listOfStages = new File(listOfStagesFileName);
		if (!listOfStages.exists()) {
			LOGGER.info("\n[Guidance] File to store the tasks list: " + listOfStagesFileName);
		} else {
			PrintWriter writer = new PrintWriter(listOfStages);
			writer.print("");
			writer.close();
			listOfStages.delete();
			LOGGER.info("\n[Guidance] File to store the tasks list (overwritten): " + listOfStagesFileName);
	
		}
		if (!listOfStages.createNewFile()) {
			LOGGER.error("[Guidance] Error on main cannot create list of stages file: " + listOfStages);
		}
		return listOfStages;
	}
	
	/**
	 * Method to print a all the environment variables of the system. It is only
	 * used for debug
	 */
	public static void printEnVariables(Logger LOGGER) {
		long freeMemory = Runtime.getRuntime().freeMemory() / 1_048_576;
		long totalMemory = Runtime.getRuntime().totalMemory() / 1_048_576;
		long maxMemory = Runtime.getRuntime().maxMemory() / 1_048_576;

		LOGGER.debug("JVM freeMemory: " + freeMemory);
		LOGGER.debug("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
		LOGGER.debug("JVM maxMemory also equals to maximum heap size of JVM   : " + maxMemory);

		Map<String, String> env = System.getenv();
		LOGGER.debug("--------------------------------------");
		LOGGER.debug("Environmental Variables in Master:");
		for (Entry<String, String> envVar : env.entrySet()) {
			LOGGER.debug(envVar.getKey() + " = " + envVar.getValue());
		}
		LOGGER.debug("--------------------------------------");
		ArrayList<String> objects = new ArrayList<>();
		for (int ii = 0; ii < 10000000; ii++) {
			objects.add((Integer.toString(10 * 2_710)));
		}

		freeMemory = Runtime.getRuntime().freeMemory() / 1_048_576;
		totalMemory = Runtime.getRuntime().totalMemory() / 1_048_576;
		maxMemory = Runtime.getRuntime().maxMemory() / 1_048_576;

		LOGGER.debug("Used Memory in JVM: " + (maxMemory - freeMemory));
		LOGGER.debug("freeMemory in JVM: " + freeMemory);
		LOGGER.debug("totalMemory in JVM shows current size of java heap : " + totalMemory);
		LOGGER.debug("maxMemory in JVM: " + maxMemory);
	}
	
	/**
	 * Verifies that the given environment variable is correctly defined
	 * 
	 * @param envVar
	 * @throws GuidanceEnvironmentException
	 */
	public static void verifyEnvVarDefined(String envVar) throws GuidanceEnvironmentException {
		if (envVar == null) {
			throw new GuidanceEnvironmentException("[Guidance] Error, " + envVar
					+ " environment variable in not present in .bashrc. You must define it properly");
		}
	}

	/**
	 * Verifies that the given environment variable is correctly defined and points
	 * to a valid directory
	 * 
	 * @param envVar
	 * @throws GuidanceEnvironmentException
	 */
	public static void verify(String envVar) throws GuidanceEnvironmentException {
		verifyEnvVarDefined(envVar);

		File f = new File(envVar);
		if (!f.exists() || f.isDirectory()) {
			throw new GuidanceEnvironmentException("[Guidance] Error, " + envVar
					+ " does not exist or it is not a binary file. Please check your .bashrc");
		}
	}
	

	/**
	 * Method to print the general information of Guidance
	 * 
	 */
	public static void printGuidancePackageVersion(Logger LOGGER, String GUIDANCE_VERSION, String GUIDANCE_BUILDNUMBER) {
		LOGGER.warn("[Guidance] *****************************************************************");
		LOGGER.warn("[Guidance] ** This is the Guidance framework to performing imputation,    **");
		LOGGER.warn("[Guidance] ** GWAS and Phenotype analysis of large scale GWAS datasets.   **");
		LOGGER.warn("[Guidance] ** Version: Guidance " + GUIDANCE_VERSION + " build " + GUIDANCE_BUILDNUMBER
				+ "           **");
		LOGGER.warn("[Guidance] ** Date release: 20-Jul-2016                                   **");
		LOGGER.warn("[Guidance] ** Contact: http://cg.bsc.es/guidance                          **");
		LOGGER.warn("[Guidance] *****************************************************************\n");
	}
}