/*
 *  Copyright 2002-2017 Barcelona Supercomputing Center (www.bsc.es)
 *  Life Science Department, 
 *  Computational Genomics Group (http://www.bsc.es/life-sciences/computational-genomics)
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *
 *
 *  Last update: $LastChangedDate: 2017-14-08 11:36:54 +0100 (Mon, 14 Ago 2017) $
 *  Revision Number: $Revision: 16 $
 *  Last revision  : $LastChangedRevision: 16 $
 *  Written by     : Friman Sanchez C.
 *                 : friman.sanchez@gmail.com
 *  Modified by    : COMPSs Support
 *                 : support-compss@bsc.es
 *                
 *  Guidance web page: http://cg.bsc.es/guidance/
 *
 */

package guidance;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Date;
import java.util.LinkedList;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.text.DateFormat;
import java.text.SimpleDateFormat;

import guidance.GuidanceImpl;
import guidance.exceptions.GuidanceEnvironmentException;
import guidance.exceptions.GuidanceTaskException;
import guidance.files.AssocFiles;
import guidance.files.CombinedPanelsFiles;
import guidance.files.CommonFiles;
import guidance.files.FileUtils;
import guidance.files.ImputationFiles;
import guidance.files.MergeFiles;
import guidance.files.PhenomeAnalysisFiles;
import guidance.files.ResultsFiles;
import guidance.utils.ChromoInfo;
import guidance.utils.Headers;
import guidance.utils.ParseCmdLine;

import es.bsc.compss.api.COMPSs;

//import es.bsc.compss.api.COMPSs;

/**
 * The higher-level class of Guidance
 */
public class Guidance {

	// Package version
	private static final String GUIDANCE_VERSION;
	private static final String GUIDANCE_BUILDNUMBER;

	// Logger
	private static final Logger LOGGER = LogManager.getLogger("Console");

	// Debug mode
	private static final boolean DEBUG = false;

	// Threshold
	private static final double PVA_THRESHOLD = 5e-8;
	private static final String PVA_THRESHOLD_STR = Double.toString(PVA_THRESHOLD);

	// Execution mode variables
	private static final String FILTERED = "filtered";
	private static final String CONDENSED = "condensed";

	// Environment variables that contain the information of the applications that
	// will be used in the Guidance
	// execution
	private static final String PLINK_BINARY = System.getenv("PLINKBINARY");
	private static final String GTOOL_BINARY = System.getenv("GTOOLBINARY");
	private static final String R_SCRIPT_BIN_DIR = System.getenv("RSCRIPTBINDIR");
	private static final String R_SCRIPT_DIR = System.getenv("RSCRIPTDIR");
	private static final String QCTOOL_BINARY = System.getenv("QCTOOLBINARY");
	private static final String EAGLEBINARY = System.getenv("EAGLEBINARY");
	private static final String PLINKBINARY = System.getenv("PLINKBINARY");
	private static final String SHAPEIT_BINARY = System.getenv("SHAPEITBINARY");
	private static final String IMPUTE2_BINARY = System.getenv("IMPUTE2BINARY");
	private final static String MINIMAC_BINARY = System.getenv("MINIMACBINARY");
	private static final String SNPTEST_BINARY = System.getenv("SNPTESTBINARY");
	private static final String JAVA_HOME = System.getenv("JAVA_HOME");
	private static final String SAMTOOLSBINARY = System.getenv("SAMTOOLSBINARY");

	private static final String SEX1 = GuidanceImpl.getSex1();
	private static final String SEX2 = GuidanceImpl.getSex2();
	private static final String NO_SEX = GuidanceImpl.getNoSex();

	private static File listOfStages;
	private static ArrayList<String> listOfCommands = new ArrayList<>();

	static {
		// Load Guidance version and build number
		String guidanceVersion = "";
		String guidanceBuildnumber = "";
		try {
			Properties props = new Properties();
			props.load(Guidance.class.getResourceAsStream("/version.properties"));
			guidanceVersion = props.getProperty("guidance.version");
			guidanceBuildnumber = props.getProperty("guidance.build");
		} catch (IOException ioe) {
			LOGGER.error("ERROR: Cannot load guidance version and buildnumber. Skip", ioe);
		} finally {
			GUIDANCE_VERSION = guidanceVersion;
			GUIDANCE_BUILDNUMBER = guidanceBuildnumber;
		}
	}

	/**
	 * GUIDANCE main code
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		// Verify that all environment variables have been defined correctly
		verifyEnvVar();

		// Print information of Guidance version
		printGuidancePackageVersion();
		if (DEBUG) {
			printEnVariables();
		}

		// Get the input arguments
		ParseCmdLine parsingArgs = new ParseCmdLine(args);

		// Verify and print the status of each stage
		printStagesStatus(parsingArgs);
		LOGGER.info("\n[Guidance] Verifyed stages status.");

		// Get the file name where the list of commands is going to be saved
		// (listOfStagesFile)
		String listOfStagesFileName = parsingArgs.getListOfStagesFile();

		// Verify whether the file exists or not.
		listOfStages = new File(listOfStagesFileName);
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

		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		String datestring = dateFormat.format(date);

		listOfCommands.add("####################################################################");
		listOfCommands.add("# List of tasks executed by Guidance workflow");
		listOfCommands.add("# Date: " + datestring);
		listOfCommands.add("# Parameters of the execution: ");
		listOfCommands.add("####################################################################");

		// Get the names of the reference panel to be used in the execution
		ArrayList<String> rpanelTypes = new ArrayList<>(parsingArgs.getRpanelTypes());
		String outDir = parsingArgs.getOutDir();

		LOGGER.info("[Guidance] We print testTypes information");
		int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
		for (int kk = 0; kk < numberOfTestTypes; kk++) {
			String tmpTestType = parsingArgs.getTestTypeName(kk);
			String tmpResponseVar = parsingArgs.getResponseVar(kk);
			String tmpCovariables = parsingArgs.getCovariables(kk);
			LOGGER.info("[Guidance] " + tmpTestType + " = " + tmpResponseVar + ":" + tmpCovariables);
		}

		// Main code of the work flow:
		doMixed(parsingArgs, outDir, rpanelTypes);

		// Finally, we print the commands in the output file defined for this.
		flushCommands();
		LOGGER.info("[Guidance] Everything is working with Guidance, just wait...");
	}

	private static void flushCommands() throws IOException {

		try (BufferedWriter writer = new BufferedWriter(new FileWriter(listOfStages, true))) {
			while (!listOfCommands.isEmpty()) {
				String str = listOfCommands.get(0);
				listOfCommands.remove(0);
				writer.write(str);
				writer.newLine();
				writer.newLine();
			}

			// Close the file with the list of commands...
			writer.flush();
			writer.close();
		}
	}

	/**
	 * Method to print a all the environment variables of the system. It is only
	 * used for debug
	 */
	private static void printEnVariables() {
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
	 * Method that performs the complete work flow when "mixed" type of GWAS is
	 * chosen.
	 * 
	 * @param parsingArgs
	 * @param outDir
	 * @param rpanelTypes
	 * @param generalChromoInfo
	 * @param listOfCommands
	 * @throws IOException
	 * @throws GuidanceTaskException
	 */
	private static void doMixed(ParseCmdLine parsingArgs, String outDir, List<String> rpanelTypes)
			throws IOException, GuidanceTaskException {
		// Create some general objects
		int startChr = parsingArgs.getStart();
		int endChr = parsingArgs.getEnd();
		// int endChrNormal = parsingArgs.getEndNormal();

		String exclCgatFlag = parsingArgs.getExclCgatSnp();
		// String exclSVFlag = parsingArgs.getExclSVSnp();

		// double mafThreshold = parsingArgs.getMafThreshold();
		// double infoThreshold = parsingArgs.getInfoThreshold();
		// double hweCohortThreshold = parsingArgs.getHweCohortThreshold();
		// double hweCasesThreshold = parsingArgs.getHweCasesThreshold();
		// double hweControlsThreshold = parsingArgs.getHweControlsThreshold();

		// String mafThresholdS = Double.toString(mafThreshold);
		// String infoThresholdS = Double.toString(infoThreshold);
		// String hweCohortThresholdS = Double.toString(hweCohortThreshold);
		// String hweCasesThresholdS = Double.toString(hweCasesThreshold);
		// String hweControlsThresholdS = Double.toString(hweControlsThreshold);

		int chunkSize = parsingArgs.getChunkSize();

		String inputFormat = null;

		// Create the names for the common files
		CommonFiles commonFilesInfo = new CommonFiles(parsingArgs, outDir);

		// Create the whole directory structure.
		LOGGER.info("[Guidance] Creating the directory structures for the outputs...");
		FileUtils.createDirStructure(parsingArgs, outDir, rpanelTypes, startChr, endChr);

		// Create the names for mixed files
		ImputationFiles imputationFilesInfo = new ImputationFiles(parsingArgs, outDir, rpanelTypes);

		// Create the names for Association files
		AssocFiles assocFilesInfo = new AssocFiles(parsingArgs, outDir, rpanelTypes);

		// Create the names for Merge files
		MergeFiles mergeFilesInfo = new MergeFiles(parsingArgs, outDir, rpanelTypes);

		// Create the names for Results Files. Take into account this class it to
		// generate file name for results
		// of each testType and each rPanelType. The results for combined panels or
		// phenome analysis are created
		// by other class.
		ResultsFiles resultsFilesInfo = new ResultsFiles(parsingArgs, outDir, rpanelTypes);

		CombinedPanelsFiles combinedPanelsFilesInfo = new CombinedPanelsFiles(parsingArgs, outDir, rpanelTypes);

		PhenomeAnalysisFiles phenomeAnalysisFilesInfo = new PhenomeAnalysisFiles(parsingArgs, outDir, rpanelTypes);

		LOGGER.info("[Guidance] Directory structures have been created...");

		// Create all file names used in the workflow.
		// Also, define which of them will be temporal or permanent.
		setFinalStatusForCommonFiles(parsingArgs, commonFilesInfo);

		LOGGER.info("[Guidance] All the common file names used during the workflow have been defined...");

		setFinalStatusForImputationFiles(parsingArgs, imputationFilesInfo, rpanelTypes);

		LOGGER.info("[Guidance] All the imputation file names used during the workflow have been defined...");

		setFinalStatusForAssocFiles(parsingArgs, assocFilesInfo, rpanelTypes);

		LOGGER.info("[Guidance] All the names used during the workflow have been defined...");

		// COMPSs.barrier();

		// The number of Chromos to process is endChr - startChr + 1;
		for (int chr = startChr; chr <= endChr; chr++) {
			// We get the output pairs file name for mixed
			String mixedPairsFile = commonFilesInfo.getPairsFile(chr);
			String mixedSampleFile = commonFilesInfo.getSampleFile(chr);

			String theChromo = Integer.toString(chr);
			String mixedGenFile = null;
			String bedFile = null;
			String bimFile = null;
			String famFile = null;
			String mixedBedFile = null;
			String mixedBimFile = null;
			String mixedFamFile = null;
			String mixedBedToBedLogFile = null;

			String gmapFile = "";
			String phasingTool = parsingArgs.getPhasingTool();
			// String imputationTool = parsingArgs.getImputationTool();

			if (phasingTool.equals("shapeit")) {
				gmapFile = parsingArgs.getGmapDir() + "/" + parsingArgs.getGmapFileName(chr);

			} else if (phasingTool.equals("eagle")) {
				gmapFile = parsingArgs.getGmapDir() + "/" + parsingArgs.getGmapFileNameEagle();

				String mixedGmapFile = commonFilesInfo.getGmapFile(chr);

				if (parsingArgs.getStageStatus("phasingBed") == 1) {
					doCreateSplitedFiles(parsingArgs, gmapFile, mixedGmapFile, theChromo);
				}

			}

			String mixedPhasingHapsFile = commonFilesInfo.getPhasingHapsFile(chr);
			String mixedPhasingSampleFile = commonFilesInfo.getPhasingSampleFile(chr);
			String mixedPhasingLogFile = commonFilesInfo.getPhasingLogFile(chr);
			String mixedExcludedSnpsFile = commonFilesInfo.getExcludedSnpsFile(chr);

			String mixedFilteredHaplotypesLogFile = commonFilesInfo.getFilteredHaplotypesLogFile(chr);
			String mixedFilteredHaplotypesVcfFile = commonFilesInfo.getFilteredHaplotypesVcfFile(chr);
			String mixedFilteredHaplotypesVcfFileBgzip = commonFilesInfo.getFilteredHaplotypesVcfFileBgzip(chr);
			String mixedFilteredHaplotypesVcfFileBgzipIndexed = commonFilesInfo
					.getFilteredHaplotypesVcfFileBgzipIndexed(chr);

			// Check if the input are in GEN o BED format
			inputFormat = parsingArgs.getInputFormat();
			if (inputFormat.equals("BED")) {

				bedFile = commonFilesInfo.getBedFile();
				bimFile = commonFilesInfo.getBimFile();
				famFile = commonFilesInfo.getFamFile();

				String mixedSplitChr23MalesLogFile = "";
				String mixedSplitChr23FemalesLogFile = "";
				String mixedChr23BedMalesFile = "";
				String mixedChr23BimMalesFile = "";
				String mixedChr23FamMalesFile = "";
				String mixedChr23BedFemalesFile = "";
				String mixedChr23BimFemalesFile = "";
				String mixedChr23FamFemalesFile = "";

				mixedBedFile = commonFilesInfo.getByChrBedFile(chr);
				mixedBimFile = commonFilesInfo.getByChrBimFile(chr);
				mixedFamFile = commonFilesInfo.getByChrFamFile(chr);
				mixedBedToBedLogFile = commonFilesInfo.getBedToBedLogFile(chr);

				doConvertFromBedToBed(parsingArgs, bedFile, bimFile, famFile, mixedBedFile, mixedBimFile, mixedFamFile,
						mixedBedToBedLogFile, theChromo);

				if (theChromo.equals("23")) {

					mixedSplitChr23MalesLogFile = commonFilesInfo.getSplitChr23MalesLogFile();
					mixedSplitChr23FemalesLogFile = commonFilesInfo.getSplitChr23FemalesLogFile();
					mixedChr23BedMalesFile = commonFilesInfo.getChr23BedMalesFile();
					mixedChr23BimMalesFile = commonFilesInfo.getChr23BimMalesFile();
					mixedChr23FamMalesFile = commonFilesInfo.getChr23FamMalesFile();
					mixedChr23BedFemalesFile = commonFilesInfo.getChr23BedFemalesFile();
					mixedChr23BimFemalesFile = commonFilesInfo.getChr23BimFemalesFile();
					mixedChr23FamFemalesFile = commonFilesInfo.getChr23FamFemalesFile();

					doSplitChr23(parsingArgs, mixedBedFile, mixedBimFile, mixedFamFile, mixedSplitChr23MalesLogFile,
							mixedChr23BedMalesFile, mixedChr23BimMalesFile, mixedChr23FamMalesFile, SEX1, theChromo);

					doSplitChr23(parsingArgs, mixedBedFile, mixedBimFile, mixedFamFile, mixedSplitChr23FemalesLogFile,
							mixedChr23BedFemalesFile, mixedChr23BimFemalesFile, mixedChr23FamFemalesFile, SEX2,
							theChromo);

				}

				// Create the RsId list of SNPs that are AT, TA, CG, or GC
				// In that case, because inputType is BED we pass the newBimFile

				doCreateRsIdList(parsingArgs, mixedBimFile, exclCgatFlag, mixedPairsFile, inputFormat);

				String responseVar = parsingArgs.getResponseVar(0);
				String covariables = parsingArgs.getCovariables(0);

				if (!theChromo.equals("23")) {

					doPhasingBed(parsingArgs, theChromo, mixedBedFile, mixedBimFile, mixedFamFile, gmapFile,
							mixedSampleFile, covariables, responseVar, mixedExcludedSnpsFile, mixedPhasingHapsFile,
							mixedPhasingSampleFile, mixedPhasingLogFile, mixedFilteredHaplotypesLogFile,
							mixedFilteredHaplotypesVcfFile, mixedFilteredHaplotypesVcfFileBgzip,
							mixedFilteredHaplotypesVcfFileBgzipIndexed, exclCgatFlag, NO_SEX);

				} else if (theChromo.equals("23")) {

					String mixedPhasingHapsMalesFile = commonFilesInfo.getPhasingHapsMalesFile();
					String mixedPhasingHapsFemalesFile = commonFilesInfo.getPhasingHapsFemalesFile();

					String mixedPhasingSampleMalesFile = commonFilesInfo.getPhasingSampleMalesFile();
					String mixedPhasingSampleFemalesFile = commonFilesInfo.getPhasingSampleFemalesFile();

					String mixedPhasingLogMalesFile = commonFilesInfo.getPhasingLogMalesFile();
					String mixedPhasingLogFemalesFile = commonFilesInfo.getPhasingLogFemalesFile();

					String mixedExcludedSnpsMalesFile = commonFilesInfo.getExcludedSnpsMalesFile();
					String mixedExcludedSnpsFemalesFile = commonFilesInfo.getExcludedSnpsFemalesFile();

					String mixedFilteredHaplotypesLogMalesFile = commonFilesInfo.getFilteredHaplotypesLogMalesFile();
					String mixedFilteredHaplotypesLogFemalesFile = commonFilesInfo
							.getFilteredHaplotypesLogFemalesFile();

					String mixedFilteredHaplotypesVcfMalesFile = commonFilesInfo.getFilteredHaplotypesVcfMalesFile();
					String mixedFilteredHaplotypesVcfFemalesFile = commonFilesInfo
							.getFilteredHaplotypesVcfFemalesFile();

					String mixedFilteredHaplotypesVcfMalesFileBgzip = commonFilesInfo
							.getFilteredHaplotypesVcfMalesFileBgzip();
					String mixedFilteredHaplotypesVcfFemalesFileBgzip = commonFilesInfo
							.getFilteredHaplotypesVcfFemalesFileBgzip();

					String mixedFilteredHaplotypesVcfMalesFileBgzipIndexed = commonFilesInfo
							.getFilteredHaplotypesVcfMalesFileBgzipIndexed();
					String mixedFilteredHaplotypesVcfFemalesFileBgzipIndexed = commonFilesInfo
							.getFilteredHaplotypesVcfFemalesFileBgzipIndexed();

					doPhasingBed(parsingArgs, theChromo, mixedChr23BedMalesFile, mixedChr23BimMalesFile,
							mixedChr23FamMalesFile, gmapFile, mixedSampleFile, covariables, responseVar,
							mixedExcludedSnpsMalesFile, mixedPhasingHapsMalesFile, mixedPhasingSampleMalesFile,
							mixedPhasingLogMalesFile, mixedFilteredHaplotypesLogMalesFile,
							mixedFilteredHaplotypesVcfMalesFile, mixedFilteredHaplotypesVcfMalesFileBgzip,
							mixedFilteredHaplotypesVcfMalesFileBgzipIndexed, exclCgatFlag, SEX1);

					doPhasingBed(parsingArgs, theChromo, mixedChr23BedFemalesFile, mixedChr23BimFemalesFile,
							mixedChr23FamFemalesFile, gmapFile, mixedSampleFile, covariables, responseVar,
							mixedExcludedSnpsFemalesFile, mixedPhasingHapsFemalesFile, mixedPhasingSampleFemalesFile,
							mixedPhasingLogFemalesFile, mixedFilteredHaplotypesLogFemalesFile,
							mixedFilteredHaplotypesVcfFemalesFile, mixedFilteredHaplotypesVcfFemalesFileBgzip,
							mixedFilteredHaplotypesVcfFemalesFileBgzipIndexed, exclCgatFlag, SEX2);

				}

			} else { // The inputFormat is GEN

				mixedGenFile = commonFilesInfo.getGenFile(chr);
				// Task
				doCreateRsIdList(parsingArgs, mixedGenFile, exclCgatFlag, mixedPairsFile, inputFormat);

				doPhasing(parsingArgs, theChromo, mixedGenFile, mixedSampleFile, gmapFile, mixedExcludedSnpsFile,
						mixedPhasingHapsFile, mixedPhasingSampleFile, mixedPhasingLogFile,
						mixedFilteredHaplotypesLogFile, mixedFilteredHaplotypesVcfFile,
						mixedFilteredHaplotypesVcfFileBgzip, mixedFilteredHaplotypesVcfFileBgzipIndexed, exclCgatFlag);

			} // End of inputFormat GEN

			flushCommands();
		}

		flushCommands();
		// COMPSs.barrier();

		// The number of Chromos to process is endChr - startChr + 1;
		for (int chr = startChr; chr <= endChr; chr++) {

			String gmapFile = parsingArgs.getGmapDir() + "/" + parsingArgs.getGmapFileName(chr);
			String mixedGmapFile = commonFilesInfo.getGmapFile(chr);

			// ***********************************************
			// * Now perform the imputation tasks *
			// * over the data for the different reference *
			// * panels that the user wants to include. The *
			// * list of the panels to be used are stored in *
			// * the List array rpanelTypes *
			// ***********************************************
			for (int panel = 0; panel < rpanelTypes.size(); panel++) {
				// String rpanelDir = parsingArgs.getRpanelDir(kk);

				int maxSize = ChromoInfo.getMaxSize(chr);
				int minSize = ChromoInfo.getMinSize(chr);
				int lim1 = minSize;
				int lim2 = lim1 + chunkSize - 1;
				for (int j = minSize; j < maxSize; j = j + chunkSize) {
					makeImputationPerChunk(parsingArgs, chr, lim1, lim2, panel, gmapFile, mixedGmapFile,
							imputationFilesInfo, commonFilesInfo);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
					flushCommands();
				}
			} // End for panel types

			flushCommands();

		} // End for chromosomes

		// *******************************************************************
		// * COMPSs API Call to wait for all tasks *
		// * comment out following line to include the synchronization point *
		// *******************************************************************
		// COMPSs.waitForAllTasks();
		// COMPSs.barrier();

		// Now we continue with the association
		int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
		for (int test = 0; test < numberOfTestTypes; test++) {
			for (int panel = 0; panel < rpanelTypes.size(); panel++) {
				for (int chr = startChr; chr <= endChr; chr++) {
					int minSize = ChromoInfo.getMinSize(chr);
					int maxSize = ChromoInfo.getMaxSize(chr);
					int lim1 = minSize;
					int lim2 = lim1 + chunkSize - 1;
					String panelName = rpanelTypes.get(panel);
					/*
					 * for(j=minSize; j<maxSize; j = j + chunkSize) {
					 * makeAssociationPerChunk(parsingArgs,tt, kk, i, lim1, lim2,
					 * imputationFilesInfo, commonFilesInfo, assocFilesInfo, listOfCommands);
					 * 
					 * lim1 = lim1 + chunkSize; lim2 = lim2 + chunkSize; }
					 * 
					 * // Now we perform the merge of chunks for each chromosome
					 * makeMergeOfChunks(parsingArgs, listOfCommands, tt, kk, i, minSize, maxSize,
					 * chunkSize, assocFilesInfo, mergeFilesInfo);
					 * 
					 * // Then, we have to filterByAll the last merged File of each chromomose.
					 * String theLastReducedFile = mergeFilesInfo.getTheLastReducedFile(tt, kk, i);
					 * String filteredByAllFile = mergeFilesInfo.getFilteredByAllFile(tt, kk, i);
					 * String condensedFile = mergeFilesInfo.getCondensedFile(tt, kk,i);
					 * 
					 * doFilterByAll(parsingArgs, listOfCommands, theLastReducedFile,
					 * filteredByAllFile, condensedFile);
					 */

					for (int j = minSize; j < maxSize; j = j + chunkSize) {
						if (DEBUG) {
							String testName = parsingArgs.getTestTypeName(test);
							LOGGER.debug("Making association for TEST " + testName + ", PANEL " + panelName
									+ ", CHROMO " + chr + ", CHUNK " + lim1 + " - " + lim2);
						}
						String testName = parsingArgs.getTestTypeName(test);
						System.out.println("Making association for TEST " + testName + ", PANEL " + panelName
								+ ", CHROMO " + chr + ", CHUNK " + lim1 + " - " + lim2);
						makeAssociationPerChunk(parsingArgs, test, panel, panelName, chr, lim1, lim2,
								imputationFilesInfo, commonFilesInfo, assocFilesInfo);

						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}

					System.out.println("Association per chunk finished successfully for chromosome " + chr
							+ " and panel " + panelName);

					flushCommands();

					if (chr == 23) {
						// Now we perform the merge of chunks for each chromosome
						if (parsingArgs.getStageStatus("jointFilteredByAllFiles") == 1) {
							makeMergeOfChunksSex(parsingArgs, test, SEX1, panel, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, FILTERED);

							makeMergeOfChunksSex(parsingArgs, test, SEX2, panel, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, FILTERED);
						}
						if (parsingArgs.getStageStatus("jointCondensedFiles") == 1) {
							makeMergeOfChunksSex(parsingArgs, test, SEX1, panel, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, CONDENSED);

							makeMergeOfChunksSex(parsingArgs, test, SEX2, panel, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, CONDENSED);
						}
					} else {
						// Now we perform the merge of chunks for each chromosome
						if (parsingArgs.getStageStatus("jointFilteredByAllFiles") == 1) {
							makeMergeOfChunks(parsingArgs, test, panel, chr, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, FILTERED);
						}
						if (parsingArgs.getStageStatus("jointCondensedFiles") == 1) {
							makeMergeOfChunks(parsingArgs, test, panel, chr, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, CONDENSED);
						}
					}
				} // End for Chromo

				// Now we have to joint the condensedFiles of each chromosome. There is not
				// problem if chr23 is being
				// processed, because
				// the format of condensedFiles is the same for all chromosome.
				makeJointCondensedFiles(parsingArgs, test, panel, startChr, endChr, mergeFilesInfo);

				// Now we have to joint the filteredByAllFiles of each chromosome. Here there is
				// an additional
				// complexity due to the chr23,
				// because the file format changes if chr23 is being processed. This situation
				// is taken into account
				// inside the next function.
				String rpanelName = rpanelTypes.get(panel);
				makeJointFilteredByAllFiles(parsingArgs, test, panel, rpanelName, startChr, endChr, mergeFilesInfo);

				// Generate TopHits
				String lastCondensedFile = mergeFilesInfo.getFinalCondensedFile(test, panel);
				String lastFilteredByAllFile = mergeFilesInfo.getFinalFilteredByAllFile(test, panel);
				String topHitsResults = resultsFilesInfo.getTopHitsFile(test, panel);
				String filteredByAllXMalesFile;
				String filteredByAllXFemalesFile;
				if (endChr == 23) {
					filteredByAllXMalesFile = mergeFilesInfo.getAdditionalFilteredByAllXMalesFile(test, panel, 0);
					filteredByAllXFemalesFile = mergeFilesInfo.getAdditionalFilteredByAllXFemalesFile(test, panel, 0);
				} else {
					filteredByAllXMalesFile = lastFilteredByAllFile;
					filteredByAllXFemalesFile = lastFilteredByAllFile;
				}

				// doGenerateTopHits(parsingArgs, lastFilteredByAllFile,
				// filteredByAllXMalesFile,
				// filteredByAllXFemalesFile, topHitsResults, PVA_THRESHOLD_STR);

				// Generate QQManhattan Plots
				String qqPlotPdfFile = resultsFilesInfo.getQqPlotPdfFile(test, panel);
				String qqPlotTiffFile = resultsFilesInfo.getQqPlotTiffFile(test, panel);
				String manhattanPlotPdfFile = resultsFilesInfo.getManhattanPdfFile(test, panel);
				String manhattanPlotTiffFile = resultsFilesInfo.getManhattanTiffFile(test, panel);
				// String correctedPvaluesFile = resultsFilesInfo.getCorrectedPvaluesFile(test,
				// panel);

				// doGenerateQQManhattanPlots(parsingArgs, lastCondensedFile, qqPlotPdfFile,
				// manhattanPlotPdfFile,
				// qqPlotTiffFile, manhattanPlotTiffFile);

				flushCommands();
				// COMPSs.barrier();

			} // End for refPanels

			// Now we continue with the combining of the results of the different reference
			// panels.
			// It is done if the refPanelCombine flag is true.
			makeCombinePanels(parsingArgs, assocFilesInfo, mergeFilesInfo, combinedPanelsFilesInfo, rpanelTypes, test);

		} // End for test types

		if (1 < numberOfTestTypes) {
			// makePhenotypeAnalysis(parsingArgs, mergeFilesInfo, resultsFilesInfo,
			// phenomeAnalysisFilesInfo, rpanelTypes);
		} else {
			LOGGER.info("\n[Guidance] No cross-phenotype analysis. Only one phenotype available");
		}

		// makePhenotypeAnalysis(parsingArgs, mergeFilesInfo, resultsFilesInfo,
		// phenomeAnalysisFilesInfo, rpanelTypes,
		// listOfCommands);

		// -------------------------->
		// makeMergingAndSummaryResults(parsingArgs, rpanelTypes, assocFilesInfo,
		// startChr, endChr, listOfCommands);

		// makeCombinePanels(parsingArgs, rpanelTypes, assocFilesInfo,
		// numOfChromosToProcess,
		// refPanelCombine, mafThreshold, infoThreshold, hweCohortThreshold,
		// hweCasesThreshold, hweControlsThreshold, listOfCommands);

		// Final part: Summary, qqplot and Manhatan plots generation.
		/*
		 * makeMergingAndDrawingResults(parsingArgs, rpanelTypes, assocFilesInfo,
		 * numOfChromosToProcess, refPanelCombine, mafThreshold, infoThreshold,
		 * hweCohortThreshold, hweCasesThreshold, hweControlsThreshold, startChr,
		 * endChr, listOfCommands);
		 */
		LOGGER.info("\n[Guidance] All tasks are in execution, please wait...");

		/**
		 * Now it is a good moment to start with the cleaning and compression of the
		 * temporal files. It should be done if variable removeTemporalFiles and
		 * compressFiles where enabled. We compress and clean in this order: -
		 * commonFilesInfo - imputationFilesInfo - assocFilesInfo Due to that currently
		 * COMPSs does not allow us to delete de temporal files, then we will compress
		 * them and then we create a text file that list the files that we want to
		 * clean. After the end of the execution, the user can delete them.
		 */
		/*
		 * try{ compressCommonFiles(parsingArgs, commonFilesInfo); } catch (Exception
		 * e){ LOGGER.error("[Guidance] Exception compressing commonFilesInfo."); }
		 * 
		 * try{ compressImputationFiles(parsingArgs, generalChromoInfo, rpanelTypes,
		 * imputationFilesInfo); } catch (Exception e){
		 * LOGGER.error("[Guidance] Exception compressing imputationFilesInfo."); }
		 * 
		 * try{ compressAssocFiles(parsingArgs, generalChromoInfo, rpanelTypes,
		 * assocFilesInfo); } catch (Exception e){
		 * LOGGER.error("[Guidance] Exception compressing imputationFilesInfo."); }
		 * 
		 * // Now we delete files try{ deleteCommonFiles(parsingArgs, commonFilesInfo);
		 * } catch (Exception e){
		 * LOGGER.error("[Guidance] Exception deleting commonFilesInfo."); }
		 * 
		 * try{ deleteImputationFiles(parsingArgs, generalChromoInfo, rpanelTypes,
		 * imputationFilesInfo); } catch (Exception e){
		 * LOGGER.error("[Guidance] Exception deleting imputationFilesInfo."); }
		 * 
		 * try{ deleteAssocFiles(parsingArgs, generalChromoInfo, rpanelTypes,
		 * assocFilesInfo); } catch (Exception e){
		 * LOGGER.error("[Guidance] Exception deleting imputationFilesInfo."); }
		 */
	}

	/**
	 * Method that generates all the tasks for imputation and association.
	 * 
	 * @param parsingArgs
	 * @param chrNumber
	 * @param lim1
	 * @param lim2
	 * @param panelIndex
	 * @param gmapFile
	 * @param imputationFilesInfo
	 * @param commonFilesInfo
	 * @param listOfCommands
	 */
	private static void makeImputationPerChunk(ParseCmdLine parsingArgs, int chrNumber, int lim1, int lim2,
			int panelIndex, String gmapFile, String mixedGmapFile, ImputationFiles imputationFilesInfo,
			CommonFiles commonFilesInfo) {

		int chunkSize = parsingArgs.getChunkSize();
		double infoThreshold = parsingArgs.getInfoThreshold();

		String rpanelDir = parsingArgs.getRpanelDir(panelIndex);

		String lim1S = Integer.toString(lim1);
		String lim2S = Integer.toString(lim2);
		String chrS = Integer.toString(chrNumber);
		String infoThresholdS = Double.toString(infoThreshold);
		String imputationTool = parsingArgs.getImputationTool();
		String phasingTool = parsingArgs.getPhasingTool();

		if (imputationTool.equals("impute")) {
			String knownHapFileName = parsingArgs.getRpanelHapFileName(panelIndex, chrNumber);
			String knownHapFile = rpanelDir + File.separator + knownHapFileName;

			String legendFileName = parsingArgs.getRpanelLegFileName(panelIndex, chrNumber);
			String legendFile = rpanelDir + File.separator + legendFileName;

			String mixedSampleFile = "";
			String mixedPhasingHapsFile = "";
			String mixedPhasingSampleFile = "";
			String mixedPairsFile = commonFilesInfo.getPairsFile(chrNumber);
			String mixedImputeFile = "";
			String mixedImputeFileInfo = "";
			String mixedImputeFileSummary = "";
			String mixedImputeFileWarnings = "";
			// String mixedImputeLogFile = "";
			String mixedFilteredFile = "";
			String mixedFilteredLogFile = "";
			String mixedFilteredRsIdFile = "";

			if (!chrS.equals("23")) {
				mixedSampleFile = commonFilesInfo.getSampleFile(chrNumber);
				mixedPhasingHapsFile = commonFilesInfo.getPhasingHapsFile(chrNumber);
				mixedPhasingSampleFile = commonFilesInfo.getPhasingSampleFile(chrNumber);
				mixedImputeFile = imputationFilesInfo.getImputedFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
				mixedImputeFileInfo = imputationFilesInfo.getImputedInfoFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				mixedImputeFileSummary = imputationFilesInfo.getImputedSummaryFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				mixedImputeFileWarnings = imputationFilesInfo.getImputedWarningsFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				/*
				 * mixedImputeLogFile = imputationFilesInfo.getImputedLogFile(panelIndex,
				 * chrNumber, lim1, lim2, chunkSize);
				 */
				mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
				mixedFilteredLogFile = imputationFilesInfo.getFilteredLogFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				mixedFilteredRsIdFile = imputationFilesInfo.getFilteredRsIdFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
			}

			// We create the list of rsId that are greater than or equal to the
			// infoThreshold value
			if (phasingTool.equals("eagle")) {

				if (!chrS.equals("23")) {
					doImputationWithImpute(parsingArgs, chrS, mixedGmapFile, knownHapFile, legendFile,
							mixedPhasingHapsFile, mixedSampleFile, lim1S, lim2S, mixedPairsFile, mixedImputeFile,
							mixedImputeFileInfo, mixedImputeFileSummary, mixedImputeFileWarnings, NO_SEX);
					doFilterByInfo(parsingArgs, mixedImputeFileInfo, mixedFilteredRsIdFile, infoThresholdS, NO_SEX);
					doQctoolS(parsingArgs, mixedImputeFile, mixedFilteredRsIdFile, mixedFilteredFile,
							mixedFilteredLogFile, NO_SEX);

				} else if (chrS.equals("23")) {

					String mixedPhasingHapsMalesFile = commonFilesInfo.getPhasingHapsMalesFile();
					String mixedPhasingHapsFemalesFile = commonFilesInfo.getPhasingHapsFemalesFile();

					String mixedPhasingSampleMalesFile = commonFilesInfo.getPhasingSampleMalesFile();
					String mixedPhasingSampleFemalesFile = commonFilesInfo.getPhasingSampleFemalesFile();

					String mixedImputeMalesFile = imputationFilesInfo.getImputedMalesFile(panelIndex, chrNumber, lim1,
							lim2, chunkSize);
					String mixedImputeFemalesFile = imputationFilesInfo.getImputedFemalesFile(panelIndex, chrNumber,
							lim1, lim2, chunkSize);

					String mixedImputeMalesFileInfo = imputationFilesInfo.getImputedInfoMalesFile(panelIndex, chrNumber,
							lim1, lim2, chunkSize);
					String mixedImputeFemalesFileInfo = imputationFilesInfo.getImputedInfoFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					String mixedImputeMalesFileSummary = imputationFilesInfo.getImputedSummaryMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedImputeFemalesFileSummary = imputationFilesInfo.getImputedSummaryFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					String mixedImputeMalesFileWarnings = imputationFilesInfo.getImputedWarningsMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedImputeFemalesFileWarnings = imputationFilesInfo
							.getImputedWarningsFemalesFile(panelIndex, chrNumber, lim1, lim2, chunkSize);

					/*
					 * String mixedImputeLogMalesFile =
					 * imputationFilesInfo.getImputedLogMalesFile(panelIndex, chrNumber, lim1, lim2,
					 * chunkSize); String mixedImputeLogFemalesFile =
					 * imputationFilesInfo.getImputedLogFemalesFile(panelIndex, chrNumber, lim1,
					 * lim2, chunkSize);
					 */

					String mixedFilteredMalesFile = imputationFilesInfo.getFilteredMalesFile(panelIndex, lim1, lim2,
							chunkSize);
					String mixedFilteredFemalesFile = imputationFilesInfo.getFilteredFemalesFile(panelIndex, lim1, lim2,
							chunkSize);

					String mixedFilteredLogMalesFile = imputationFilesInfo.getFilteredLogMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedFilteredLogFemalesFile = imputationFilesInfo.getFilteredLogFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					String mixedFilteredRsIdMalesFile = imputationFilesInfo.getFilteredRsIdMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedFilteredRsIdFemalesFile = imputationFilesInfo.getFilteredRsIdFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					doImputationWithImpute(parsingArgs, chrS, mixedGmapFile, knownHapFile, legendFile,
							mixedPhasingHapsMalesFile, mixedPhasingSampleMalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeMalesFile, mixedImputeMalesFileInfo, mixedImputeMalesFileSummary,
							mixedImputeMalesFileWarnings, SEX1);

					doFilterByInfo(parsingArgs, mixedImputeMalesFileInfo, mixedFilteredRsIdMalesFile, infoThresholdS,
							SEX1);

					doQctoolS(parsingArgs, mixedImputeMalesFile, mixedFilteredRsIdMalesFile, mixedFilteredMalesFile,
							mixedFilteredLogMalesFile, SEX1);

					doImputationWithImpute(parsingArgs, chrS, mixedGmapFile, knownHapFile, legendFile,
							mixedPhasingHapsFemalesFile, mixedPhasingSampleFemalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeFemalesFile, mixedImputeFemalesFileInfo, mixedImputeFemalesFileSummary,
							mixedImputeFemalesFileWarnings, SEX2);

					doFilterByInfo(parsingArgs, mixedImputeFemalesFileInfo, mixedFilteredRsIdFemalesFile,
							infoThresholdS, SEX2);

					doQctoolS(parsingArgs, mixedImputeFemalesFile, mixedFilteredRsIdFemalesFile,
							mixedFilteredFemalesFile, mixedFilteredLogFemalesFile, SEX2);
				}

			} else if (phasingTool.equals("shapeit")) {
				if (!chrS.equals("23")) {
					doImputationWithImpute(parsingArgs, chrS, gmapFile, knownHapFile, legendFile, mixedPhasingHapsFile,
							mixedPhasingSampleFile, lim1S, lim2S, mixedPairsFile, mixedImputeFile, mixedImputeFileInfo,
							mixedImputeFileSummary, mixedImputeFileWarnings, NO_SEX);
					doFilterByInfo(parsingArgs, mixedImputeFileInfo, mixedFilteredRsIdFile, infoThresholdS, NO_SEX);
					doQctoolS(parsingArgs, mixedImputeFile, mixedFilteredRsIdFile, mixedFilteredFile,
							mixedFilteredLogFile, NO_SEX);

				} else if (chrS.equals("23")) {
					String mixedPhasingHapsMalesFile = commonFilesInfo.getPhasingHapsMalesFile();
					String mixedPhasingHapsFemalesFile = commonFilesInfo.getPhasingHapsFemalesFile();

					String mixedPhasingSampleMalesFile = commonFilesInfo.getPhasingSampleMalesFile();
					String mixedPhasingSampleFemalesFile = commonFilesInfo.getPhasingSampleFemalesFile();

					String mixedImputeMalesFile = imputationFilesInfo.getImputedMalesFile(panelIndex, chrNumber, lim1,
							lim2, chunkSize);
					String mixedImputeFemalesFile = imputationFilesInfo.getImputedFemalesFile(panelIndex, chrNumber,
							lim1, lim2, chunkSize);

					String mixedImputeMalesFileInfo = imputationFilesInfo.getImputedInfoMalesFile(panelIndex, chrNumber,
							lim1, lim2, chunkSize);
					String mixedImputeFemalesFileInfo = imputationFilesInfo.getImputedInfoFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					String mixedImputeMalesFileSummary = imputationFilesInfo.getImputedSummaryMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedImputeFemalesFileSummary = imputationFilesInfo.getImputedSummaryFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					String mixedImputeMalesFileWarnings = imputationFilesInfo.getImputedWarningsMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedImputeFemalesFileWarnings = imputationFilesInfo
							.getImputedWarningsFemalesFile(panelIndex, chrNumber, lim1, lim2, chunkSize);

					/*
					 * String mixedImputeLogMalesFile =
					 * imputationFilesInfo.getImputedLogMalesFile(panelIndex, chrNumber, lim1, lim2,
					 * chunkSize); String mixedImputeLogFemalesFile =
					 * imputationFilesInfo.getImputedLogFemalesFile(panelIndex, chrNumber, lim1,
					 * lim2, chunkSize);
					 */

					String mixedFilteredMalesFile = imputationFilesInfo.getFilteredMalesFile(panelIndex, lim1, lim2,
							chunkSize);
					String mixedFilteredFemalesFile = imputationFilesInfo.getFilteredFemalesFile(panelIndex, lim1, lim2,
							chunkSize);

					String mixedFilteredLogMalesFile = imputationFilesInfo.getFilteredLogMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedFilteredLogFemalesFile = imputationFilesInfo.getFilteredLogFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					String mixedFilteredRsIdMalesFile = imputationFilesInfo.getFilteredRsIdMalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);
					String mixedFilteredRsIdFemalesFile = imputationFilesInfo.getFilteredRsIdFemalesFile(panelIndex,
							chrNumber, lim1, lim2, chunkSize);

					doImputationWithImpute(parsingArgs, chrS, gmapFile, knownHapFile, legendFile,
							mixedPhasingHapsMalesFile, mixedPhasingSampleMalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeMalesFile, mixedImputeMalesFileInfo, mixedImputeMalesFileSummary,
							mixedImputeMalesFileWarnings, SEX1);

					doFilterByInfo(parsingArgs, mixedImputeMalesFileInfo, mixedFilteredRsIdMalesFile, infoThresholdS,
							SEX1);

					doQctoolS(parsingArgs, mixedImputeMalesFile, mixedFilteredRsIdMalesFile, mixedFilteredMalesFile,
							mixedFilteredLogMalesFile, SEX1);

					doImputationWithImpute(parsingArgs, chrS, gmapFile, knownHapFile, legendFile,
							mixedPhasingHapsFemalesFile, mixedPhasingSampleFemalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeFemalesFile, mixedImputeFemalesFileInfo, mixedImputeFemalesFileSummary,
							mixedImputeFemalesFileWarnings, SEX2);

					doFilterByInfo(parsingArgs, mixedImputeFemalesFileInfo, mixedFilteredRsIdFemalesFile,
							infoThresholdS, SEX2);

					doQctoolS(parsingArgs, mixedImputeFemalesFile, mixedFilteredRsIdFemalesFile,
							mixedFilteredFemalesFile, mixedFilteredLogFemalesFile, SEX2);
				}
			}

		} else if (imputationTool.equals("minimac")) {
			String refVcfFileName = parsingArgs.getRpanelVCFFileName(panelIndex, chrNumber);
			String refVcfFile = rpanelDir + "/" + refVcfFileName;

			/*
			 * String mixedFilteredHaplotypesVcfFileBgzip =
			 * commonFilesInfo.getFilteredHaplotypesVcfFileBgzip(chrNumber); String
			 * mixedImputeMMDoseVCFFile =
			 * imputationFilesInfo.getImputedMMDoseVCFFile(panelIndex, chrNumber, lim1,
			 * lim2, chunkSize); String mixedImputeMMInfoFile =
			 * imputationFilesInfo.getImputedMMInfoFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedImputeMMErateFile =
			 * imputationFilesInfo.getImputedMMErateFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedImputeMMRecFile =
			 * imputationFilesInfo.getImputedMMRecFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedImputeMMM3VCFFile =
			 * imputationFilesInfo.getImputedMMM3VCFFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedImputeMMLogFile =
			 * imputationFilesInfo.getImputedMMLogFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedFilteredFile =
			 * imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedFilteredLogFile =
			 * imputationFilesInfo.getFilteredLogFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedFilteredRsIdFile =
			 * imputationFilesInfo.getFilteredRsIdFile(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedImputeFileBgzip =
			 * imputationFilesInfo.getImputedFileBgzip(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); String mixedImputeMMM3VCFFileBgzip =
			 * imputationFilesInfo.getImputedMMM3VCFFileBgzip(panelIndex, chrNumber, lim1,
			 * lim2, chunkSize); String mixedImputeFileTbi =
			 * imputationFilesInfo.getImputedFileTbi(panelIndex, chrNumber, lim1, lim2,
			 * chunkSize); // String mixedImputeMMM3VCFFileTbi =
			 * imputationFilesInfo.getImputedMMM3VCFFileTbi(panelIndex, chrNumber, lim1,
			 * lim2, chunkSize);
			 */
			if (!chrS.equals("23")) {

				String mixedFilteredHaplotypesVcfFileBgzip = commonFilesInfo
						.getFilteredHaplotypesVcfFileBgzip(chrNumber);
				String mixedImputeMMDoseVCFFile = imputationFilesInfo.getImputedMMDoseVCFFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedImputeMMInfoFile = imputationFilesInfo.getImputedMMInfoFile(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String mixedImputeMMErateFile = imputationFilesInfo.getImputedMMErateFile(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String mixedImputeMMRecFile = imputationFilesInfo.getImputedMMRecFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				String mixedImputeMMM3VCFFile = imputationFilesInfo.getImputedMMM3VCFFile(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String mixedImputeMMLogFile = imputationFilesInfo.getImputedMMLogFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				String mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				String mixedFilteredLogFile = imputationFilesInfo.getFilteredLogFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				String mixedFilteredRsIdFile = imputationFilesInfo.getFilteredRsIdFile(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String mixedImputeFileBgzip = imputationFilesInfo.getImputedFileBgzip(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				String mixedImputeMMM3VCFFileBgzip = imputationFilesInfo.getImputedMMM3VCFFileBgzip(panelIndex,
						chrNumber, lim1, lim2, chunkSize);
				String mixedImputeFileTbi = imputationFilesInfo.getImputedFileTbi(panelIndex, chrNumber, lim1, lim2,
						chunkSize);

				doImputationWithMinimac(parsingArgs, refVcfFile, mixedFilteredHaplotypesVcfFileBgzip, chrS, lim1S,
						lim2S, mixedImputeMMDoseVCFFile, mixedImputeMMInfoFile, mixedImputeMMErateFile,
						mixedImputeMMRecFile, mixedImputeMMM3VCFFile, mixedImputeMMLogFile, mixedImputeFileBgzip,
						mixedImputeMMM3VCFFileBgzip, mixedImputeFileTbi, NO_SEX);

				doFilterByInfo(parsingArgs, mixedImputeMMInfoFile, mixedFilteredRsIdFile, infoThresholdS, NO_SEX);

				doQctoolS(parsingArgs, mixedImputeFileBgzip, mixedFilteredRsIdFile, mixedFilteredFile,
						mixedFilteredLogFile, NO_SEX);

			}

			else if (chrS.equals("23")) {
				String mixedFilteredHaplotypesVcfMalesFileBgzip = commonFilesInfo
						.getFilteredHaplotypesVcfMalesFileBgzip();
				String mixedFilteredHaplotypesVcfFemalesFileBgzip = commonFilesInfo
						.getFilteredHaplotypesVcfFemalesFileBgzip();

				String mixedImputeMMDoseVCFMalesFile = imputationFilesInfo.getImputedMMDoseVCFMalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMDoseVCFFemalesFile = imputationFilesInfo.getImputedMMDoseVCFFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMInfoMalesFile = imputationFilesInfo.getImputedMMInfoMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedImputeMMInfoFemalesFile = imputationFilesInfo.getImputedMMInfoFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMErateMalesFile = imputationFilesInfo.getImputedMMErateMalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);
				String mixedImputeMMErateFemalesFile = imputationFilesInfo.getImputedMMErateFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMRecMalesFile = imputationFilesInfo.getImputedMMRecMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedImputeMMRecFemalesFile = imputationFilesInfo.getImputedMMRecFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMM3VCFMalesFile = imputationFilesInfo.getImputedMMM3VCFMalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);
				String mixedImputeMMM3VCFFemalesFile = imputationFilesInfo.getImputedMMM3VCFFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMLogMalesFile = imputationFilesInfo.getImputedMMLogMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedImputeMMLogFemalesFile = imputationFilesInfo.getImputedMMLogFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedFilteredMalesFile = imputationFilesInfo.getFilteredMalesFile(panelIndex, lim1, lim2,
						chunkSize);
				String mixedFilteredFemalesFile = imputationFilesInfo.getFilteredFemalesFile(panelIndex, lim1, lim2,
						chunkSize);

				String mixedFilteredLogMalesFile = imputationFilesInfo.getFilteredLogMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedFilteredLogFemalesFile = imputationFilesInfo.getFilteredLogFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedFilteredRsIdMalesFile = imputationFilesInfo.getFilteredRsIdMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedFilteredRsIdFemalesFile = imputationFilesInfo.getFilteredRsIdFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMalesFileBgzip = imputationFilesInfo.getImputedMalesFileBgzip(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedImputeFemalesFileBgzip = imputationFilesInfo.getImputedFemalesFileBgzip(panelIndex,
						chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMMM3VCFMalesFileBgzip = imputationFilesInfo
						.getImputedMMM3VCFMalesFileBgzip(panelIndex, chrNumber, lim1, lim2, chunkSize);
				String mixedImputeMMM3VCFFemalesFileBgzip = imputationFilesInfo
						.getImputedMMM3VCFFemalesFileBgzip(panelIndex, chrNumber, lim1, lim2, chunkSize);

				String mixedImputeMalesFileTbi = imputationFilesInfo.getImputedMalesFileTbi(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String mixedImputeFemalesFileTbi = imputationFilesInfo.getImputedFemalesFileTbi(panelIndex, chrNumber,
						lim1, lim2, chunkSize);

				doImputationWithMinimac(parsingArgs, refVcfFile, mixedFilteredHaplotypesVcfMalesFileBgzip, chrS, lim1S,
						lim2S, mixedImputeMMDoseVCFMalesFile, mixedImputeMMInfoMalesFile, mixedImputeMMErateMalesFile,
						mixedImputeMMRecMalesFile, mixedImputeMMM3VCFMalesFile, mixedImputeMMLogMalesFile,
						mixedImputeMalesFileBgzip, mixedImputeMMM3VCFMalesFileBgzip, mixedImputeMalesFileTbi, SEX1);

				doFilterByInfo(parsingArgs, mixedImputeMMInfoMalesFile, mixedFilteredRsIdMalesFile, infoThresholdS,
						SEX1);

				doQctoolS(parsingArgs, mixedImputeMalesFileBgzip, mixedFilteredRsIdMalesFile, mixedFilteredMalesFile,
						mixedFilteredLogMalesFile, SEX1);

				doImputationWithMinimac(parsingArgs, refVcfFile, mixedFilteredHaplotypesVcfFemalesFileBgzip, chrS,
						lim1S, lim2S, mixedImputeMMDoseVCFFemalesFile, mixedImputeMMInfoFemalesFile,
						mixedImputeMMErateFemalesFile, mixedImputeMMRecFemalesFile, mixedImputeMMM3VCFFemalesFile,
						mixedImputeMMLogFemalesFile, mixedImputeFemalesFileBgzip, mixedImputeMMM3VCFFemalesFileBgzip,
						mixedImputeFemalesFileTbi, SEX2);

				doFilterByInfo(parsingArgs, mixedImputeMMInfoFemalesFile, mixedFilteredRsIdFemalesFile, infoThresholdS,
						SEX2);

				doQctoolS(parsingArgs, mixedImputeFemalesFileBgzip, mixedFilteredRsIdFemalesFile,
						mixedFilteredFemalesFile, mixedFilteredLogFemalesFile, SEX2);

			}

		} else {
			System.err.println(
					"\t[makeImputationPerChunk]: Error, the imputation tool " + imputationTool + " is not allowed...");
			System.exit(1);
		}

	}

	/**
	 * Method that generates all the tasks for association.
	 * 
	 * @param parsingArgs
	 * @param testTypeIndex
	 * @param panelIndex
	 * @param chrNumber
	 * @param lim1
	 * @param lim2
	 * @param imputationFilesInfo
	 * @param commonFilesInfo
	 * @param assocFilesInfo
	 * @param listOfCommands
	 */
	private static void makeAssociationPerChunk(ParseCmdLine parsingArgs, int testTypeIndex, int panelIndex,
			String rpanelName, int chrNumber, int lim1, int lim2, ImputationFiles imputationFilesInfo,
			CommonFiles commonFilesInfo, AssocFiles assocFilesInfo) {

		String imputationTool = parsingArgs.getImputationTool();

		int chunkSize = parsingArgs.getChunkSize();
		double mafThreshold = parsingArgs.getMafThreshold();
		double infoThreshold = parsingArgs.getInfoThreshold();
		double hweCohortThreshold = parsingArgs.getHweCohortThreshold();
		double hweCasesThreshold = parsingArgs.getHweCasesThreshold();
		double hweControlsThreshold = parsingArgs.getHweControlsThreshold();

		String responseVar = parsingArgs.getResponseVar(testTypeIndex);
		String covariables = parsingArgs.getCovariables(testTypeIndex);

		String chrS = Integer.toString(chrNumber);

		String mafThresholdS = Double.toString(mafThreshold);
		String infoThresholdS = Double.toString(infoThreshold);
		String hweCohortThresholdS = Double.toString(hweCohortThreshold);
		String hweCasesThresholdS = Double.toString(hweCasesThreshold);
		String hweControlsThresholdS = Double.toString(hweControlsThreshold);

		// String inputFormat = parsingArgs.getInputFormat();

		// String mixedBimFile = null;
		// String mixedGenFile = null;
		// String exclCgatFlag = parsingArgs.getExclCgatSnp();

		// String inputFormat = parsingArgs.getInputFormat();
		// String exclCgatFlag = parsingArgs.getExclCgatSnp();

		if (chrNumber == 23) {
			String snptestOutMalesFile = assocFilesInfo.getSnptestOutMalesFile(testTypeIndex, panelIndex, lim1, lim2,
					chunkSize);
			String snptestOutFemalesFile = assocFilesInfo.getSnptestOutFemalesFile(testTypeIndex, panelIndex, lim1,
					lim2, chunkSize);
			String snptestLogMalesFile = assocFilesInfo.getSnptestLogMalesFile(testTypeIndex, panelIndex, lim1, lim2,
					chunkSize);
			String snptestLogFemalesFile = assocFilesInfo.getSnptestLogFemalesFile(testTypeIndex, panelIndex, lim1,
					lim2, chunkSize);
			String mixedFilteredMalesFile = imputationFilesInfo.getFilteredMalesFile(panelIndex, lim1, lim2, chunkSize);
			String mixedFilteredFemalesFile = imputationFilesInfo.getFilteredFemalesFile(panelIndex, lim1, lim2,
					chunkSize);
			String mixedPhasingSampleMalesFile = commonFilesInfo.getPhasingSampleMalesFile();
			String mixedPhasingSampleFemalesFile = commonFilesInfo.getPhasingSampleFemalesFile();

			doSnptest(parsingArgs, chrS, mixedFilteredMalesFile, mixedPhasingSampleMalesFile, snptestOutMalesFile,
					snptestLogMalesFile, responseVar, covariables);

			doSnptest(parsingArgs, chrS, mixedFilteredFemalesFile, mixedPhasingSampleFemalesFile, snptestOutFemalesFile,
					snptestLogFemalesFile, responseVar, covariables);

			if (imputationTool.equals("impute")) {
				String mixedImputeMalesFileInfo = imputationFilesInfo.getImputedInfoMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String mixedImputeFemalesFileInfo = imputationFilesInfo.getImputedInfoFemalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String summaryMalesFile = assocFilesInfo.getSummaryMalesFile(testTypeIndex, panelIndex, lim1, lim2,
						chunkSize);
				String summaryFemalesFile = assocFilesInfo.getSummaryFemalesFile(testTypeIndex, panelIndex, lim1, lim2,
						chunkSize);

				doCollectSummary(parsingArgs, chrS, mixedImputeMalesFileInfo, snptestOutMalesFile, summaryMalesFile,
						mafThresholdS, infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);

				doCollectSummary(parsingArgs, chrS, mixedImputeFemalesFileInfo, snptestOutFemalesFile,
						summaryFemalesFile, mafThresholdS, infoThresholdS, hweCohortThresholdS, hweCasesThresholdS,
						hweControlsThresholdS);

				String assocMalesFilteredByAll = assocFilesInfo.getSummaryFilteredMalesFile(testTypeIndex, panelIndex,
						lim1, lim2, chunkSize);
				String assocMalesCondensed = assocFilesInfo.getSummaryCondensedMalesFile(testTypeIndex, panelIndex,
						lim1, lim2, chunkSize);

				doFilterByAll(parsingArgs, summaryFemalesFile, assocMalesFilteredByAll, assocMalesCondensed, SEX1,
						rpanelName);

				String assocFemalesFilteredByAll = assocFilesInfo.getSummaryFilteredFemalesFile(testTypeIndex,
						panelIndex, lim1, lim2, chunkSize);
				String assocFemalesCondensed = assocFilesInfo.getSummaryCondensedFemalesFile(testTypeIndex, panelIndex,
						lim1, lim2, chunkSize);

				doFilterByAll(parsingArgs, summaryFemalesFile, assocFemalesFilteredByAll, assocFemalesCondensed, SEX2,
						rpanelName);

			} else if (imputationTool.equals("minimac")) {
				String mixedImputedMMInfoMalesFile = imputationFilesInfo.getImputedMMInfoMalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);
				String mixedImputedMMInfoFemalesFile = imputationFilesInfo.getImputedMMInfoFemalesFile(panelIndex,
						chrNumber, lim1, lim2, chunkSize);
				String summaryMalesFile = assocFilesInfo.getSummaryMalesFile(testTypeIndex, panelIndex, lim1, lim2,
						chunkSize);
				String summaryFemalesFile = assocFilesInfo.getSummaryFemalesFile(testTypeIndex, panelIndex, lim1, lim2,
						chunkSize);

				doCollectSummary(parsingArgs, chrS, mixedImputedMMInfoMalesFile, snptestOutMalesFile, summaryMalesFile,
						mafThresholdS, infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);

				doCollectSummary(parsingArgs, chrS, mixedImputedMMInfoFemalesFile, snptestOutFemalesFile,
						summaryFemalesFile, mafThresholdS, infoThresholdS, hweCohortThresholdS, hweCasesThresholdS,
						hweControlsThresholdS);

				String assocMalesFilteredByAll = assocFilesInfo.getSummaryFilteredMalesFile(testTypeIndex, panelIndex,
						lim1, lim2, chunkSize);
				String assocMalesCondensed = assocFilesInfo.getSummaryCondensedMalesFile(testTypeIndex, panelIndex,
						lim1, lim2, chunkSize);

				doFilterByAll(parsingArgs, summaryMalesFile, assocMalesFilteredByAll, assocMalesCondensed, SEX1,
						rpanelName);

				String assocFemalesFilteredByAll = assocFilesInfo.getSummaryFilteredFemalesFile(testTypeIndex,
						panelIndex, lim1, lim2, chunkSize);
				String assocFemalesCondensed = assocFilesInfo.getSummaryCondensedFemalesFile(testTypeIndex, panelIndex,
						lim1, lim2, chunkSize);

				doFilterByAll(parsingArgs, summaryFemalesFile, assocFemalesFilteredByAll, assocFemalesCondensed, SEX2,
						rpanelName);
			}

		} else {
			String snptestOutFile = assocFilesInfo.getSnptestOutFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String snptestLogFile = assocFilesInfo.getSnptestLogFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedPhasingSampleFile = commonFilesInfo.getPhasingSampleFile(chrNumber);

			doSnptest(parsingArgs, chrS, mixedFilteredFile, mixedPhasingSampleFile, snptestOutFile, snptestLogFile,
					responseVar, covariables);

			if (imputationTool.equals("impute")) {
				String mixedImputeFileInfo = imputationFilesInfo.getImputedInfoFile(panelIndex, chrNumber, lim1, lim2,
						chunkSize);
				String summaryFile = assocFilesInfo.getSummaryFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
						chunkSize);

				doCollectSummary(parsingArgs, chrS, mixedImputeFileInfo, snptestOutFile, summaryFile, mafThresholdS,
						infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);

				String assocFilteredByAll = assocFilesInfo.getSummaryFilteredFile(testTypeIndex, panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String assocCondensed = assocFilesInfo.getSummaryCondensedFile(testTypeIndex, panelIndex, chrNumber,
						lim1, lim2, chunkSize);

				doFilterByAll(parsingArgs, summaryFile, assocFilteredByAll, assocCondensed, NO_SEX, rpanelName);

			} else if (imputationTool.equals("minimac")) {
				String mixedImputedMMInfoFile = imputationFilesInfo.getImputedMMInfoFile(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String summaryFile = assocFilesInfo.getSummaryFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
						chunkSize);

				doCollectSummary(parsingArgs, chrS, mixedImputedMMInfoFile, snptestOutFile, summaryFile, mafThresholdS,
						infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);

				String assocFilteredByAll = assocFilesInfo.getSummaryFilteredFile(testTypeIndex, panelIndex, chrNumber,
						lim1, lim2, chunkSize);
				String assocCondensed = assocFilesInfo.getSummaryCondensedFile(testTypeIndex, panelIndex, chrNumber,
						lim1, lim2, chunkSize);

				doFilterByAll(parsingArgs, summaryFile, assocFilteredByAll, assocCondensed, NO_SEX, rpanelName);
			}

		}

	}

	/**
	 * Method to perform the merging of chunks for chromosome 23.
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param ttIndex
	 * @param rpanelIndex
	 * @param chr
	 * @param minSize
	 * @param maxSize
	 * @param chunkSize
	 * @param assocFilesInfo
	 * @param mergeFilesInfo
	 * @param type
	 */
	private static void makeMergeOfChunksSex(ParseCmdLine parsingArgs, int ttIndex, String sex, int rpanelIndex,
			int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo,
			String type) {
		int lim1 = minSize;
		int lim2 = lim1 + chunkSize - 1;

		int numberOfChunks = maxSize / chunkSize;
		int module = maxSize % chunkSize;
		if (module != 0)
			numberOfChunks++;

		int indexA = 0;
		int indexC = 0;
		String reducedA = null;
		String reducedB = null;
		String reducedC = null;
		if (sex.equals(SEX1)) {
			if (type.equals(FILTERED)) {
				String filteredByAllFile = mergeFilesInfo.getFilteredByAllMalesFile(ttIndex, rpanelIndex);
				for (int processedChunks = 0; processedChunks < 2 * numberOfChunks
						- 2; processedChunks = processedChunks + 2) {
					if (processedChunks < numberOfChunks) {
						reducedA = assocFilesInfo.getSummaryFilteredMalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedA = mergeFilesInfo.getReducedFilteredMalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}
					if (processedChunks < numberOfChunks - 1) {
						reducedB = assocFilesInfo.getSummaryFilteredMalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedB = mergeFilesInfo.getReducedFilteredMalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}

					if (processedChunks == 2 * numberOfChunks - 4) {
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, filteredByAllFile);
						indexC++;
					} else {
						reducedC = mergeFilesInfo.getReducedFilteredMalesFile(ttIndex, rpanelIndex, indexC);
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
						indexC++;
					}

					// Clean intermediate files
					// File fA = new File(reducedA);
					// fA.delete();
					// File fB = new File(reducedB);
					// fB.delete();

				} // End for Chunks

			} else if (type.equals(CONDENSED)) {
				String condensedFile = mergeFilesInfo.getCondensedMalesFile(ttIndex, rpanelIndex);
				for (int processedChunks = 0; processedChunks < 2 * numberOfChunks
						- 2; processedChunks = processedChunks + 2) {
					if (processedChunks < numberOfChunks) {
						reducedA = assocFilesInfo.getSummaryCondensedMalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedA = mergeFilesInfo.getReducedCondensedMalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}
					if (processedChunks < numberOfChunks - 1) {
						reducedB = assocFilesInfo.getSummaryCondensedMalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedB = mergeFilesInfo.getReducedCondensedMalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}

					if (processedChunks == 2 * numberOfChunks - 4) {
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, condensedFile);
						indexC++;
					} else {
						reducedC = mergeFilesInfo.getReducedCondensedMalesFile(ttIndex, rpanelIndex, indexC);
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
						indexC++;
					}

					// Clean intermediate files
					// File fA = new File(reducedA);
					// fA.delete();
					// File fB = new File(reducedB);
					// fB.delete();

				} // End of for Chunks
			}
		} else if (sex.equals(SEX2)) {
			if (type.equals(FILTERED)) {
				// LOGGER.info("Number of chunks for testType " + ttIndex + " | rpanel " +
				// rpanelIndex + " |chr " +
				// chr + " " + numberOfChunks);
				String filteredByAllFile = mergeFilesInfo.getFilteredByAllFemalesFile(ttIndex, rpanelIndex);
				for (int processedChunks = 0; processedChunks < 2 * numberOfChunks
						- 2; processedChunks = processedChunks + 2) {
					if (processedChunks < numberOfChunks) {
						reducedA = assocFilesInfo.getSummaryFilteredFemalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedA = mergeFilesInfo.getReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}
					if (processedChunks < numberOfChunks - 1) {
						reducedB = assocFilesInfo.getSummaryFilteredFemalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedB = mergeFilesInfo.getReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}

					if (processedChunks == 2 * numberOfChunks - 4) {
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, filteredByAllFile);
						indexC++;
					} else {
						reducedC = mergeFilesInfo.getReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexC);
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
						indexC++;
					}

					// Clean intermediate files
					// File fA = new File(reducedA);
					// fA.delete();
					// File fB = new File(reducedB);
					// fB.delete();

				} // End for Chunks

			} else if (type.equals(CONDENSED)) {
				String condensedFile = mergeFilesInfo.getCondensedFemalesFile(ttIndex, rpanelIndex);
				for (int processedChunks = 0; processedChunks < 2 * numberOfChunks
						- 2; processedChunks = processedChunks + 2) {
					if (processedChunks < numberOfChunks) {
						reducedA = assocFilesInfo.getSummaryCondensedFemalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedA = mergeFilesInfo.getReducedCondensedFemalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}
					if (processedChunks < numberOfChunks - 1) {
						reducedB = assocFilesInfo.getSummaryCondensedFemalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedB = mergeFilesInfo.getReducedCondensedFemalesFile(ttIndex, rpanelIndex, indexA);
						indexA++;
					}

					if (processedChunks == 2 * numberOfChunks - 4) {
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, condensedFile);
						indexC++;
					} else {
						reducedC = mergeFilesInfo.getReducedCondensedFemalesFile(ttIndex, rpanelIndex, indexC);
						doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
						indexC++;
					}

					// Clean intermediate files
					// File fA = new File(reducedA);
					// fA.delete();
					// File fB = new File(reducedB);
					// fB.delete();

				} // End of for Chunks
			}
		}
	}

	/**
	 * Method to perform the merging of chunks for each chromosome.
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param ttIndex
	 * @param rpanelIndex
	 * @param chr
	 * @param minSize
	 * @param maxSize
	 * @param chunkSize
	 * @param assocFilesInfo
	 * @param mergeFilesInfo
	 * @param type
	 */
	private static void makeMergeOfChunks(ParseCmdLine parsingArgs, int ttIndex, int rpanelIndex, int chr, int minSize,
			int maxSize, int chunkSize, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo, String type) {

		int lim1 = minSize;
		int lim2 = lim1 + chunkSize - 1;

		int numberOfChunks = maxSize / chunkSize;
		int module = maxSize % chunkSize;
		if (module != 0)
			numberOfChunks++;

		int indexA = 0;
		int indexC = 0;
		String reducedA = null;
		String reducedB = null;
		String reducedC = null;

		if (type.equals(FILTERED)) {
			// LOGGER.info("Number of chunks for testType " + ttIndex + " | rpanel " +
			// rpanelIndex + " |chr " +
			// chr + " " + numberOfChunks);
			String filteredByAllFile = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, chr);
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getSummaryFilteredFile(ttIndex, rpanelIndex, chr, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getReducedFilteredFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getSummaryFilteredFile(ttIndex, rpanelIndex, chr, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getReducedFilteredFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, filteredByAllFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getReducedFilteredFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End for Chunks

		} else if (type.equals(CONDENSED)) {
			String condensedFile = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, chr);
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getSummaryCondensedFile(ttIndex, rpanelIndex, chr, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getReducedCondensedFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getSummaryCondensedFile(ttIndex, rpanelIndex, chr, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getReducedCondensedFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, condensedFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getReducedCondensedFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End of for Chunks
		}
	}

	/**
	 * Method to perform the merging of chunks for each chromosome.
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param ttIndex
	 * @param rpanelIndex
	 * @param chr
	 * @param minSize
	 * @param maxSize
	 * @param chunkSize
	 * @param assocFilesInfo
	 * @param mergeFilesInfo
	 * @param type
	 */
	private static void makeMergeOfChunksCombined(ParseCmdLine parsingArgs, int ttIndex, int rpanelIndex, int chr,
			int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo,
			String type) {

		int lim1 = minSize;
		int lim2 = lim1 + chunkSize - 1;

		int numberOfChunks = maxSize / chunkSize;
		int module = maxSize % chunkSize;
		if (module != 0)
			numberOfChunks++;

		int indexA = 0;
		int indexC = 0;
		String reducedA = null;
		String reducedB = null;
		String reducedC = null;

		String filteredByAllFile;
		filteredByAllFile = mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, rpanelIndex, chr);
		String condensedFile = mergeFilesInfo.getCombinedCondensedFile(ttIndex, rpanelIndex, chr);

		if (type.equals(FILTERED)) {
			// LOGGER.info("Number of chunks for testType " + ttIndex + " | rpanel " +
			// rpanelIndex + " |chr " +
			// chr + " " + numberOfChunks);
			System.out.println("Number of chunks for testType " + ttIndex + " | rpanel " + rpanelIndex + " |chr " + chr
					+ " numberOfChunks: " + numberOfChunks);
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				System.out.println("Number of chunks for testType " + ttIndex + " | rpanel " + rpanelIndex + " |chr "
						+ chr + " processedChunks: " + processedChunks);
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getCombinedFilteredFile(ttIndex, rpanelIndex, chr, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getCombinedReducedFilteredFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getCombinedFilteredFile(ttIndex, rpanelIndex, chr, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getCombinedReducedFilteredFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, filteredByAllFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedFilteredFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End for Chunks

		} else if (type.equals(CONDENSED)) {
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getCombinedCondensedFile(ttIndex, rpanelIndex, chr, lim1, lim2,
							chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getCombinedReducedCondensedFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getCombinedCondensedFile(ttIndex, rpanelIndex, chr, lim1, lim2,
							chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getCombinedReducedCondensedFile(ttIndex, rpanelIndex, chr, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, condensedFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedCondensedFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End of for Chunks
		}
	}

	/**
	 * Method to perform the merging of chunks for each chromosome.
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param ttIndex
	 * @param rpanelIndex
	 * @param chr
	 * @param minSize
	 * @param maxSize
	 * @param chunkSize
	 * @param assocFilesInfo
	 * @param mergeFilesInfo
	 * @param type
	 */
	private static void makeMergeOfChunksCombinedSex(ParseCmdLine parsingArgs, int ttIndex, int rpanelIndex,
			int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo,
			String type) {

		int lim1 = minSize;
		int lim2 = lim1 + chunkSize - 1;

		int numberOfChunks = maxSize / chunkSize;
		int module = maxSize % chunkSize;
		if (module != 0)
			numberOfChunks++;

		int indexA = 0;
		int indexC = 0;
		String reducedA = null;
		String reducedB = null;
		String reducedC = null;

		String filteredByAllFileMales = mergeFilesInfo.getCombinedFilteredByAllMalesFile(ttIndex, rpanelIndex);
		String filteredByAllFileFemales = mergeFilesInfo.getCombinedFilteredByAllFemalesFile(ttIndex, rpanelIndex);
		String condensedFileMales = mergeFilesInfo.getCombinedCondensedMalesFile(ttIndex, rpanelIndex);

		if (type.equals(FILTERED)) {
			// LOGGER.info("Number of chunks for testType " + ttIndex + " | rpanel " +
			// rpanelIndex + " |chr " +
			// chr + " " + numberOfChunks);
			System.out.println("Number of chunks for testType " + ttIndex + " | rpanel " + rpanelIndex + " |chr "
					+ " numberOfChunks: " + numberOfChunks);
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				System.out.println("Number of chunks for testType " + ttIndex + " | rpanel " + rpanelIndex + " |chr "
						+ " processedChunks: " + processedChunks);
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getCombinedFilteredMalesFile(ttIndex, rpanelIndex, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getCombinedReducedFilteredMalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getCombinedFilteredMalesFile(ttIndex, rpanelIndex, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getCombinedReducedFilteredMalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, filteredByAllFileMales);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedFilteredMalesFile(ttIndex, rpanelIndex, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End for Chunks
			
			indexA = 0;
			indexC = 0;
			
			lim1 = minSize;
			lim2 = lim1 + chunkSize - 1;
			
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				System.out.println("Number of chunks for testType " + ttIndex + " | rpanel " + rpanelIndex + " |chr "
						+ " processedChunks: " + processedChunks);
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getCombinedFilteredFemalesFile(ttIndex, rpanelIndex, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getCombinedReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getCombinedFilteredFemalesFile(ttIndex, rpanelIndex, lim1, lim2, chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getCombinedReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, filteredByAllFileFemales);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End for Chunks

		} else if (type.equals(CONDENSED)) {
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getCombinedCondensedMalesFile(ttIndex, rpanelIndex, lim1, lim2,
							chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getCombinedReducedCondensedMalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getCombinedCondensedMalesFile(ttIndex, rpanelIndex, lim1, lim2,
							chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedB = mergeFilesInfo.getCombinedReducedCondensedMalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}

				if (processedChunks == 2 * numberOfChunks - 4) {
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, condensedFileMales);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedCondensedMalesFile(ttIndex, rpanelIndex, indexC);
					doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
					indexC++;
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End of for Chunks
		}
	}

	/**
	 * Method to perform the joint of condensed files of each rpanel
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param ttIndex
	 * @param rpanelIndex
	 * @param startChr
	 * @param endChr
	 * @param mergeFilesInfo
	 * @throws GuidanceTaskException
	 */
	private static void makeJointCondensedFiles(ParseCmdLine parsingArgs, int ttIndex, int rpanelIndex, int startChr,
			int endChr, MergeFiles mergeFilesInfo) throws GuidanceTaskException {

		int indexA = 0;
		int indexC = 0;
		int i = 0;
		String condensedA = null;
		String condensedB = null;
		String condensedC = null;
		int numberOfChrs = endChr - startChr + 1;

		if (numberOfChrs == 1) { // There is only one chr to process.
			// DO something.
			condensedA = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, startChr);

			condensedC = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexA);

			try {
				GuidanceImpl.copyFile(condensedA, condensedC);
			} catch (Exception e) {
				throw new GuidanceTaskException("Error copying the file " + e);
			}

			// doJointCondenseFiles(parsingArgs, condensedA, condensedB, condensedC);
		} else {
			if (endChr == 23)
				++numberOfChrs;
			System.out.println("From " + startChr + " to " + endChr + " with " + numberOfChrs + " chromosomes");
			for (int processedCondensed = 0; processedCondensed < 2 * numberOfChrs
					- 2; processedCondensed = processedCondensed + 2) {
				System.out.println("I'm in processedCondensed with processed " + processedCondensed);
				if (processedCondensed < numberOfChrs) {
					if (numberOfChrs == processedCondensed + 1 && endChr == 23) {
						condensedA = mergeFilesInfo.getCondensedFemalesFile(ttIndex, rpanelIndex);
					} else {
						i = startChr + processedCondensed;
						condensedA = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, i);
					}
				} else {
					condensedA = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}

				System.out.println(
						"I'm in processedCondensed 2 with processed " + processedCondensed + " and " + numberOfChrs);
				if (processedCondensed < numberOfChrs - 1) {
					if (numberOfChrs == processedCondensed + 2 && endChr == 23) {
						condensedB = mergeFilesInfo.getCondensedFemalesFile(ttIndex, rpanelIndex);
					} else {
						i = startChr + processedCondensed + 1;
						condensedB = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, i);
					}
				} else {
					System.out.println("Asking for additional condensed file and obtaining " + condensedB);
					condensedB = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}
				condensedC = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexC);

				doJointCondenseFiles(parsingArgs, condensedA, condensedB, condensedC);
				indexC++;
			} // End for(int processedCondensed=0; processedCondensed<= 2*numberOfChrs -2;
				// processedCondensed =
				// processedCondensed +2)
		}

		// Clean intermediate files
		// File fA = new File(condensedA);
		// fA.delete();
		// File fB = new File(condensedB);
		// fB.delete();
	}

	/**
	 * Method to perform the joint of filteredByAll files of each rpanel.
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param ttIndex
	 * @param rpanelIndex
	 * @param rpanelName
	 * @param startChr
	 * @param endChr
	 * @param mergeFilesInfo
	 */
	private static void makeJointFilteredByAllFiles(ParseCmdLine parsingArgs, int ttIndex, int rpanelIndex,
			String rpanelName, int startChr, int endChr, MergeFiles mergeFilesInfo) {

		int endChrNormal = endChr;
		if (startChr < endChr) {
			if (endChr != 23) {
				endChrNormal = endChr;
			} else {
				endChrNormal = endChr - 1;
			}
		}

		int indexA = 0;
		int indexC = 0;
		int processedFiltered = 0;
		String rpanelFlag = "NO";
		String filteredA = null;
		String filteredB = null;
		String filteredC = null;
		int numberOfChrs = endChrNormal - startChr + 1;

		if (numberOfChrs == 1) { // There is only one chr to process.
			rpanelFlag = "YES";
			// DO something.
			filteredA = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, startChr);
			filteredB = filteredA;

			filteredC = mergeFilesInfo.getAdditionalFilteredByAllFile(ttIndex, rpanelIndex, indexC);
			if (startChr != 23) {
				doJointFilteredByAllFiles(parsingArgs, filteredA, filteredB, filteredC, rpanelName, rpanelFlag);
			}
		} else {
			for (processedFiltered = 0; processedFiltered < 2 * numberOfChrs - 2; processedFiltered = processedFiltered
					+ 2) {
				if (processedFiltered < numberOfChrs) {
					int i = startChr + processedFiltered;
					filteredA = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, i);
				} else {
					filteredA = mergeFilesInfo.getAdditionalFilteredByAllFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}

				if (processedFiltered < numberOfChrs - 1) {
					int i = startChr + processedFiltered + 1;
					filteredB = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, i);
				} else {
					filteredB = mergeFilesInfo.getAdditionalFilteredByAllFile(ttIndex, rpanelIndex, indexA);
					indexA++;
					rpanelFlag = "YES";
				}
				filteredC = mergeFilesInfo.getAdditionalFilteredByAllFile(ttIndex, rpanelIndex, indexC);

				doJointFilteredByAllFiles(parsingArgs, filteredA, filteredB, filteredC, rpanelName, rpanelFlag);
				indexC++;
			} // End for fitlered chromosomes
				// processedFiltered +2)
		}

		// Now we process the chr 23 if this is defined in this execution
		if (endChr == 23) {
			rpanelFlag = "YES";
			filteredA = mergeFilesInfo.getFilteredByAllMalesFile(ttIndex, rpanelIndex);
			filteredC = mergeFilesInfo.getAdditionalFilteredByAllMalesFile(ttIndex, rpanelIndex);

			doJointFilteredByAllFiles(parsingArgs, filteredA, filteredA, filteredC, rpanelName, rpanelFlag);

			filteredA = mergeFilesInfo.getFilteredByAllFemalesFile(ttIndex, rpanelIndex);
			filteredC = mergeFilesInfo.getAdditionalFilteredByAllFemalesFile(ttIndex, rpanelIndex);

			doJointFilteredByAllFiles(parsingArgs, filteredA, filteredA, filteredC, rpanelName, rpanelFlag);
		}

		// Clean intermediate files
		// File fA = new File(filteredA);
		// fA.delete();
		// File fB = new File(filteredB);
		// fB.delete();
	}

	/**
	 * Method that performs the last part of the work flow corresponding to the
	 * merging, combining and drawing of results
	 * 
	 * @param parsingArgs
	 * @param assocFilesInfo
	 * @param mergeFilesInfo
	 * @param combinedPanelsFilesInfo
	 * @param rpanelTypes
	 * @param ttIndex
	 * @param listOfCommands
	 * @param generalChromoInfo
	 * @throws IOException
	 * @throws GuidanceTaskException
	 */

	public static void makeCombinePanels(ParseCmdLine parsingArgs, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo,
			CombinedPanelsFiles combinedPanelsFilesInfo, List<String> rpanelTypes, int ttIndex)
			throws IOException, GuidanceTaskException {

		final boolean refPanelCombine = parsingArgs.getRefPanelCombine();
		if (!refPanelCombine || parsingArgs.getStageStatus("combinePanelsComplex") != 1) {
			// We are not asked to combine the panels. End
			return;
		}

		// We combine the panels per chromosome
		final int startChr = parsingArgs.getStart();
		final int endChr = parsingArgs.getEnd();
		final int chunkSize = parsingArgs.getChunkSize();

		// PLACES TO STORE ALL COMBINED INFORMATION
		final String filteredCombineAll = combinedPanelsFilesInfo.getCombinedFilteredByAllFile(ttIndex);

		// We use same file when we do not need to compute the chromo 23 (doTopHits task
		// uses this information)
		String filteredCombineAllXMales = filteredCombineAll;
		String filteredCombineAllXFemales = filteredCombineAll;
		if (endChr == 23) {
			filteredCombineAllXMales = combinedPanelsFilesInfo.getCombinedFilteredByAllXMalesFile(ttIndex);
			filteredCombineAllXFemales = combinedPanelsFilesInfo.getCombinedFilteredByAllXFemalesFile(ttIndex);
		}

		final String condensedCombineAll = combinedPanelsFilesInfo.getCombinedCondensedFile(ttIndex);

		// INITIALIZE THE FILES TO STORE ALL COMBINED INFORMATION (ADD HEADER AND
		// COMPRESS)
		String filteredHeader = null;
		if (startChr < 23) {
			filteredHeader = Headers.constructHeader();
			final String plainfilteredCombineAll = filteredCombineAll.substring(0, filteredCombineAll.length() - 3);
			try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainfilteredCombineAll))) {
				writer.write(filteredHeader);
				writer.newLine();
				writer.flush();
			} catch (IOException ioe) {
				LOGGER.error("[Guidance] Exception when initializing makeCombinePanel filtered ALL file", ioe);
			}
			FileUtils.gzipFile(plainfilteredCombineAll, filteredCombineAll);
			new File(plainfilteredCombineAll).delete();
		}

		String filteredXHeader = null;
		if (endChr == 23) {
			filteredXHeader = Headers.constructHeaderX();

			final String plainfilteredCombineAllXMales = filteredCombineAllXMales.substring(0,
					filteredCombineAllXMales.length() - 3);
			try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainfilteredCombineAllXMales))) {
				writer.write(filteredXHeader);
				writer.newLine();
				writer.flush();
			} catch (IOException ioe) {
				LOGGER.error("[Guidance] Exception when initializing makeCombinePanel filteredX ALL file", ioe);
			}
			FileUtils.gzipFile(plainfilteredCombineAllXMales, filteredCombineAllXMales);
			new File(plainfilteredCombineAllXMales).delete();

			final String plainfilteredCombineAllXFemales = filteredCombineAllXFemales.substring(0,
					filteredCombineAllXFemales.length() - 3);
			try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainfilteredCombineAllXFemales))) {
				writer.write(filteredXHeader);
				writer.newLine();
				writer.flush();
			} catch (IOException ioe) {
				LOGGER.error("[Guidance] Exception when initializing makeCombinePanel filteredX ALL file", ioe);
			}
			FileUtils.gzipFile(plainfilteredCombineAllXFemales, filteredCombineAllXFemales);
			new File(plainfilteredCombineAllXFemales).delete();
		}

		final String condensedHeader = Headers.constructCondensedHeader();
		final String plainCondensedCombineAll = condensedCombineAll.substring(0, condensedCombineAll.length() - 3);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainCondensedCombineAll))) {
			writer.write(condensedHeader);
			writer.newLine();
			writer.flush();
		} catch (IOException ioe) {
			LOGGER.error("[Guidance] Exception when initializing makeCombinePanel condensed ALL file", ioe);
		}
		FileUtils.gzipFile(plainCondensedCombineAll, condensedCombineAll);
		new File(plainCondensedCombineAll).delete();

		// TODO Until here we initialize the filtered and condensed files

		LinkedList<String> filteredCombined = new LinkedList<>();
		LinkedList<String> condensedCombined = new LinkedList<>();

		// CHR LOOP
		for (int chr = startChr; chr <= endChr; chr++) {
			int minSize = ChromoInfo.getMinSize(chr);
			int maxSize = ChromoInfo.getMaxSize(chr);

			int lim1 = minSize;
			int lim2 = lim1 + chunkSize - 1;
			// int indexFC = 0;
			// int indexCC = 0;
			// int indexXFC = 0;

			// CHUNK LOOP
			for (int j = minSize; j < maxSize; j = j + chunkSize) {
				// -- FILTERED PART --
				// Construct a queue with all the filtered panels to combine
				LinkedList<String> filteredPanelsToCombine = new LinkedList<>();
				LinkedList<String> filteredPanelsToCombineMales = new LinkedList<>();
				LinkedList<String> filteredPanelsToCombineFemales = new LinkedList<>();

				for (int k = 0; k < rpanelTypes.size(); ++k) {
					if (chr == 23) {
						String filteredPanelMales = assocFilesInfo.getSummaryFilteredMalesFile(ttIndex, k, lim1, lim2,
								chunkSize);
						if (DEBUG) {
							System.out.println("[Guidance] Adding filtered file to combine " + filteredPanelMales);
						}
						filteredPanelsToCombineMales.add(filteredPanelMales);

						String filteredPanelFemales = assocFilesInfo.getSummaryFilteredFemalesFile(ttIndex, k, lim1, lim2,
								chunkSize);
						System.out.println("[Guidance] Adding filtered file to combine " + filteredPanelFemales);
						if (DEBUG) {
							System.out.println("[Guidance] Adding filtered file to combine " + filteredPanelFemales);
						}
						filteredPanelsToCombineFemales.add(filteredPanelFemales);
					} else {
						String filteredPanel = assocFilesInfo.getSummaryFilteredFile(ttIndex, k, chr, lim1, lim2,
								chunkSize);
						if (DEBUG) {
							LOGGER.debug("[Guidance] Adding filtered file to combine " + filteredPanel);
						}
						filteredPanelsToCombine.add(filteredPanel);
					}
				}

				flushCommands();

				// Combine all the filtered panels 2 by 2 until there are no remaining panels
				if (chr == 23) {
					while (!filteredPanelsToCombineMales.isEmpty()) {
						String filteredPanelA = filteredPanelsToCombineMales.poll();
						if (!filteredPanelsToCombineMales.isEmpty()) {
							String filteredPanelB = filteredPanelsToCombineMales.poll();
							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelA, lim1,
									lim2);
							// Adds A to the queue again
							filteredPanelsToCombineMales.add(filteredPanelA);

							// Deletes B since it is no longer needed
							// new File(filteredPanelB).delete();
						} else {
							String destFilteredPanelMales = assocFilesInfo.getCombinedFilteredMalesFile(ttIndex, 0,
									lim1, lim2, chunkSize);
							doCopyFile(parsingArgs, filteredPanelA, destFilteredPanelMales);
						}
					}

					while (!filteredPanelsToCombineFemales.isEmpty()) {
						String filteredPanelA = filteredPanelsToCombineFemales.poll();
						if (!filteredPanelsToCombineFemales.isEmpty()) {
							String filteredPanelB = filteredPanelsToCombineFemales.poll();
							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelA, lim1,
									lim2);
							// Adds A to the queue again
							filteredPanelsToCombineFemales.add(filteredPanelA);

							// Deletes B since it is no longer needed
							// new File(filteredPanelB).delete();
						} else {
							String destFilteredPanelFemales = assocFilesInfo.getCombinedFilteredFemalesFile(ttIndex, 0,
									lim1, lim2, chunkSize);
							doCopyFile(parsingArgs, filteredPanelA, destFilteredPanelFemales);
						}
					}
				} else {
					while (!filteredPanelsToCombine.isEmpty()) {
						String filteredPanelA = filteredPanelsToCombine.poll();
						if (!filteredPanelsToCombine.isEmpty()) {
							String filteredPanelB = filteredPanelsToCombine.poll();
							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelA, lim1,
									lim2);
							// Adds A to the queue again
							filteredPanelsToCombine.add(filteredPanelA);

							// Deletes B since it is no longer needed
							// new File(filteredPanelB).delete();
						} else {
							String destFilteredPanel = assocFilesInfo.getCombinedFilteredFile(ttIndex, 0, chr, lim1,
									lim2, chunkSize);
							doCopyFile(parsingArgs, filteredPanelA, destFilteredPanel);
						}
					}
				}

				flushCommands();
				// TODO : think how to do the condensed part
				/*
				 * 
				 * // -- CONDENSED PART -- // Construct a queue with all the condensed panels to
				 * combine LinkedList<String> condensedPanelsToCombine = new LinkedList<>();
				 * LinkedList<String> condensedPanelsToCombineMales = new LinkedList<>();
				 * LinkedList<String> condensedPanelsToCombineFemales = new LinkedList<>(); for
				 * (int k = 0; k < rpanelTypes.size(); ++k) { if (chr == 23) { String
				 * condensedPanelMales = assocFilesInfo.getSummaryCondensedMalesFile(ttIndex, k,
				 * lim1, lim2, chunkSize); if (DEBUG) {
				 * LOGGER.debug("[Guidance] Adding condensed file to combine " +
				 * condensedPanelMales); }
				 * condensedPanelsToCombineMales.add(condensedPanelMales);
				 * 
				 * String condensedPanelFemales =
				 * assocFilesInfo.getSummaryCondensedFemalesFile(ttIndex, k, lim1, lim2,
				 * chunkSize); if (DEBUG) {
				 * LOGGER.debug("[Guidance] Adding condensed file to combine " +
				 * condensedPanelFemales); }
				 * condensedPanelsToCombineFemales.add(condensedPanelFemales); } else { String
				 * condensedPanel = assocFilesInfo.getSummaryCondensedFile(ttIndex, k, chr,
				 * lim1, lim2, chunkSize); if (DEBUG) {
				 * LOGGER.debug("[Guidance] Adding condensed file to combine " +
				 * condensedPanel); } condensedPanelsToCombine.add(condensedPanel); } }
				 * 
				 * if (chr == 23) { while (!condensedPanelsToCombineMales.isEmpty()) { String
				 * condensedPanelA = condensedPanelsToCombineMales.poll(); if
				 * (!condensedPanelsToCombineMales.isEmpty()) { String condensedPanelB =
				 * condensedPanelsToCombineMales.poll(); // Filtered part: combines A and B into
				 * A if (DEBUG) { LOGGER.debug("[Guidance] Combining " + condensedPanelA +
				 * " and " + condensedPanelB); } doCombinePanelsComplex(parsingArgs,
				 * condensedPanelA, condensedPanelB, lim1, lim2); // Adds A to the queue again
				 * condensedPanelsToCombine.add(condensedPanelA);
				 * 
				 * // Deletes B since it is no longer needed // new
				 * File(condensedPanelB).delete(); } else { String destCondensedPanelMales =
				 * assocFilesInfo.getCombinedCondensedMalesFile(ttIndex, 0, lim1, lim2,
				 * chunkSize); doCopyFile(parsingArgs, condensedPanelA,
				 * destCondensedPanelMales); } }
				 * 
				 * while (!condensedPanelsToCombineFemales.isEmpty()) { String condensedPanelA =
				 * condensedPanelsToCombineFemales.poll(); if
				 * (!condensedPanelsToCombineFemales.isEmpty()) { String condensedPanelB =
				 * condensedPanelsToCombineFemales.poll(); // Filtered part: combines A and B
				 * into A if (DEBUG) { LOGGER.debug("[Guidance] Combining " + condensedPanelA +
				 * " and " + condensedPanelB); } doCombinePanelsComplex(parsingArgs,
				 * condensedPanelA, condensedPanelB, lim1, lim2); // Adds A to the queue again
				 * condensedPanelsToCombine.add(condensedPanelA);
				 * 
				 * // Deletes B since it is no longer needed // new
				 * File(condensedPanelB).delete(); } else { String destCondensedPanelFemales =
				 * assocFilesInfo.getCombinedCondensedFemalesFile(ttIndex, 0, lim1, lim2,
				 * chunkSize); doCopyFile(parsingArgs, condensedPanelA,
				 * destCondensedPanelFemales); } }
				 * 
				 * } else { // Combine all the filtered panels 2 by 2 until there are no
				 * remaining panels while (!condensedPanelsToCombine.isEmpty()) { String
				 * condensedPanelA = condensedPanelsToCombine.poll(); if
				 * (!condensedPanelsToCombine.isEmpty()) { String condensedPanelB =
				 * condensedPanelsToCombine.poll(); // Filtered part: combines A and B into A if
				 * (DEBUG) { LOGGER.debug("[Guidance] Combining " + condensedPanelA + " and " +
				 * condensedPanelB); } doCombinePanelsComplex(parsingArgs, condensedPanelA,
				 * condensedPanelB, lim1, lim2); // Adds A to the queue again
				 * condensedPanelsToCombine.add(condensedPanelA);
				 * 
				 * // Deletes B since it is no longer needed // new
				 * File(condensedPanelB).delete(); } else { String destCondensedPanel =
				 * assocFilesInfo.getCombinedCondensedFile(ttIndex, 0, chr, lim1, lim2,
				 * chunkSize); doCopyFile(parsingArgs, condensedPanelA, destCondensedPanel); } }
				 * }
				 * 
				 */

				flushCommands();
				// TODO at this point we have generated the combined files (both the condensed
				// and the filtered) for each chunk

				// -- COMBINE FILTERED TO FINAL CHUNK AND ALL FILES
				// Files are merged to the first or the last file depending on the panel Size

				/*
				 * int mergeIndex = (rpanelTypes.size() % 2 == 0) ? 0 : rpanelTypes.size() - 1;
				 * String chunkResultsFiltered = assocFilesInfo.getSummaryFilteredFile(ttIndex,
				 * mergeIndex, chr, lim1, lim2, chunkSize); if (chr != 23) { // We merge in the
				 * chunk files String filteredCombinePerChunk =
				 * combinedPanelsFilesInfo.getCombinedFilteredByAllChromoFile(ttIndex, indexFC);
				 * if (!new File(filteredCombinePerChunk).exists()) { // Initialize filtered per
				 * chunk file to default header and compress it String
				 * plainfilteredCombinePerChunk = filteredCombinePerChunk.substring(0,
				 * filteredCombinePerChunk.length() - 3); try (BufferedWriter writer = new
				 * BufferedWriter(new FileWriter(plainfilteredCombinePerChunk))) {
				 * writer.write(filteredHeader); writer.newLine(); writer.flush(); } catch
				 * (IOException ioe) { LOGGER.error(
				 * "[Guidance] Exception when initializing makeCombinePanel filteredCombinePerChunk file"
				 * , ioe); } FileUtils.gzipFile(plainfilteredCombinePerChunk,
				 * filteredCombinePerChunk); new File(plainfilteredCombinePerChunk).delete(); }
				 * // We merge in the all files if (DEBUG) {
				 * LOGGER.debug("[Guidance] Merging filtered into ALL: " + filteredCombineAll +
				 * " and " + chunkResultsFiltered); } ++indexFC; } else { // We merge in the
				 * chunk files String filteredCombineXPerChunk = combinedPanelsFilesInfo
				 * .getCombinedFilteredByAllXMalesChromoFile(ttIndex, indexXFC); if (!new
				 * File(filteredCombineXPerChunk).exists()) { // Initialize filtered per chunk
				 * file to default header and compress it String plainfilteredCombineXPerChunk =
				 * filteredCombineXPerChunk.substring(0, filteredCombineXPerChunk.length() - 3);
				 * try (BufferedWriter writer = new BufferedWriter( new
				 * FileWriter(plainfilteredCombineXPerChunk))) { writer.write(filteredXHeader);
				 * writer.newLine(); writer.flush(); } catch (IOException ioe) { LOGGER.error(
				 * "[Guidance] Exception when initializing makeCombinePanel filteredCombineXPerChunk file"
				 * , ioe); } FileUtils.gzipFile(plainfilteredCombineXPerChunk,
				 * filteredCombineXPerChunk); new File(plainfilteredCombineXPerChunk).delete();
				 * }
				 * 
				 * // We merge in the all files ++indexXFC; }
				 * 
				 * // Clean partial results // new File(chunkResultsFiltered).delete();
				 * 
				 * // -- COMBINE CONDENSED TO FINAL CHUNK AND ALL FILES String
				 * chunkResultsCondensed = assocFilesInfo.getSummaryCondensedFile(ttIndex,
				 * mergeIndex, chr, lim1, lim2, chunkSize); // We merge in the chunk files
				 * String condensedCombinePerChunk =
				 * combinedPanelsFilesInfo.getCombinedCondensedChromoFile(ttIndex, indexCC); if
				 * (!new File(condensedCombinePerChunk).exists()) { // Initialize filtered per
				 * chunk file to default header and compress it String
				 * plainCondensedCombinePerChunk = condensedCombinePerChunk.substring(0,
				 * condensedCombinePerChunk.length() - 3); try (BufferedWriter writer = new
				 * BufferedWriter(new FileWriter(plainCondensedCombinePerChunk))) {
				 * writer.write(condensedHeader); writer.newLine(); writer.flush(); } catch
				 * (IOException ioe) { LOGGER.error(
				 * "[Guidance] Exception when initializing makeCombinePanel plainCondensedCombinePerChunk file"
				 * , ioe); } FileUtils.gzipFile(plainCondensedCombinePerChunk,
				 * condensedCombinePerChunk); new File(plainCondensedCombinePerChunk).delete();
				 * }
				 * 
				 * // We merge in the all files if (DEBUG) {
				 * LOGGER.debug("[Guidance] Merging combined into ALL: " + condensedCombineAll +
				 * " and " + chunkResultsCondensed); }
				 * 
				 * ++indexCC;
				 */

				// Clean partial results
				// new File(chunkResultsCondensed).delete();

				// Increase loop variables
				lim1 = lim1 + chunkSize;
				lim2 = lim2 + chunkSize;
			} // End for chunk

			if (chr == 23) {
				makeMergeOfChunksCombinedSex(parsingArgs, ttIndex, 0, minSize, maxSize, chunkSize, assocFilesInfo,
						mergeFilesInfo, FILTERED);

				// makeMergeOfChunksCombinedSex(parsingArgs, ttIndex, 0, minSize, maxSize,
				// chunkSize, assocFilesInfo,
				// mergeFilesInfo, CONDENSED);
			} else {
				makeMergeOfChunksCombined(parsingArgs, ttIndex, 0, chr, minSize, maxSize, chunkSize, assocFilesInfo,
						mergeFilesInfo, FILTERED);

				// makeMergeOfChunksCombined(parsingArgs, ttIndex, 0, chr, minSize, maxSize,
				// chunkSize, assocFilesInfo,
				// mergeFilesInfo, CONDENSED);
			}
			// TODO look at the condensed part
			if (chr != 23) {
				String filteredByAllCurrentFile = mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, 0, chr);
				// String condensedCurrentFile =
				// mergeFilesInfo.getCombinedCondensedFile(ttIndex, 0, chr);

				filteredCombined.add(filteredByAllCurrentFile);
				// condensedCombined.add(condensedCurrentFile);
			} else {
				String filteredByAllCurrentMalesFile = mergeFilesInfo.getCombinedFilteredByAllMalesFile(ttIndex, 0);
				doCopyFile(parsingArgs, filteredByAllCurrentMalesFile, filteredCombineAllXMales);
				String filteredByAllCurrentFemalesFile = mergeFilesInfo.getCombinedFilteredByAllFemalesFile(ttIndex, 0);
				doCopyFile(parsingArgs, filteredByAllCurrentFemalesFile, filteredCombineAllXFemales);

				// String condensedCurrentMalesFile =
				// mergeFilesInfo.getCombinedCondensedMalesFile(ttIndex, 0);
				// condensedCombined.add(condensedCurrentMalesFile);

				// String condensedCurrentFemalesFile =
				// mergeFilesInfo.getCombinedCondensedFemalesFile(ttIndex, 0);
				// condensedCombined.add(condensedCurrentFemalesFile);
			}

		} // End for chromosomes

		if (filteredCombined.size() == 1) {
			String singleFilteredFile = filteredCombined.peek();
			String destinationFilteredFile = filteredCombineAll;
			doCopyFile(parsingArgs, singleFilteredFile, destinationFilteredFile);
		}

		int reduceCounter = 0;
		while (filteredCombined.size() > 1) {
			String originFilteredFileA = filteredCombined.poll();
			String originFilteredFileB = filteredCombined.poll();

			String destinationFilteredFile;

			if (filteredCombined.isEmpty()) {
				destinationFilteredFile = filteredCombineAll;
			} else {
				destinationFilteredFile = filteredCombineAll.substring(0, filteredCombineAll.length() - 7) + "_reduce_"
						+ Integer.toString(reduceCounter) + ".txt.gz";
			}

			doMergeTwoChunksUnconditional(parsingArgs, originFilteredFileA, originFilteredFileB,
					destinationFilteredFile);
			filteredCombined.add(destinationFilteredFile);

			reduceCounter += 1;
		}
		// TODO look at the condensed part
		/*
		 * if (condensedCombined.size() == 1) { String singleCondensedFile =
		 * condensedCombined.poll(); String destinationCondensedFile =
		 * condensedCombineAll; doCopyFile(parsingArgs, singleCondensedFile,
		 * destinationCondensedFile); }
		 * 
		 * reduceCounter = 0; while (condensedCombined.size() > 1) { String
		 * originCondensedFileA = condensedCombined.poll(); String originCondensedFileB
		 * = condensedCombined.poll();
		 * 
		 * String destinationCondensedFile;
		 * 
		 * if (condensedCombined.isEmpty()) { destinationCondensedFile =
		 * condensedCombineAll; } else { destinationCondensedFile =
		 * condensedCombineAll.substring(0, condensedCombineAll.length() - 7) +
		 * "_reduce_" + Integer.toString(reduceCounter) + ".txt.gz"; }
		 * 
		 * doMergeTwoChunksUnconditional(parsingArgs, originCondensedFileA,
		 * originCondensedFileB, destinationCondensedFile);
		 * condensedCombined.add(destinationCondensedFile);
		 * 
		 * reduceCounter += 1; }
		 */

		// Finally, we create topHits from filteredCombined, and QQ and Manhattan plots
		// from condensedCombined
		// TODO wait for the correct binary
		String topHitsCombinedResults = combinedPanelsFilesInfo.getTopHitsFile(ttIndex);
		// doGenerateTopHits(parsingArgs, filteredCombineAll, filteredCombineAllXMales,
		// filteredCombineAllXFemales,
		// topHitsCombinedResults, PVA_THRESHOLD_STR);

		String combinedQqPlotPdfFile = combinedPanelsFilesInfo.getQqPlotPdfFile(ttIndex);
		String combinedQqPlotTiffFile = combinedPanelsFilesInfo.getQqPlotTiffFile(ttIndex);
		String combinedManhattanPlotPdfFile = combinedPanelsFilesInfo.getManhattanPdfFile(ttIndex);
		String combinedManhattanPlotTiffFile = combinedPanelsFilesInfo.getManhattanTiffFile(ttIndex);
		// String combinedCorrectedPvaluesFile =
		// combinedPanelsFilesInfo.getCorrectedPvaluesFile(ttIndex);

		// doGenerateQQManhattanPlots(parsingArgs, condensedCombineAll,
		// combinedQqPlotPdfFile,
		// combinedManhattanPlotPdfFile, combinedQqPlotTiffFile,
		// combinedManhattanPlotTiffFile);

	}

	/**
	 * Method that performs the last part of the work flow corresponding to the
	 * phenome Analysis, combining the results of each phenoType analysis
	 * 
	 * @param parsingArgs
	 * @param mergeFilesInfo
	 * @param resultsFilesInfo
	 * @param phenomeAnalysisFilesInfo
	 * @param rpanelTypes
	 * @param listOfCommands
	 */
	private static void makePhenotypeAnalysis(ParseCmdLine parsingArgs, MergeFiles mergeFilesInfo,
			ResultsFiles resultsFilesInfo, PhenomeAnalysisFiles phenomeAnalysisFilesInfo, List<String> rpanelTypes) {

		int endChr = parsingArgs.getEnd();
		String endChrS = Integer.toString(endChr);

		String topHitsFile = null;

		String filteredByAllFile = null;
		String filteredByAllXMalesFile = null;
		String filteredByAllXFemalesFile = null;

		String phenomeFileA = null;
		String phenomeFileB = null;
		String phenomeFileC = null;

		int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
		int numberOfRpanelsTypes = rpanelTypes.size();

		// int maxPhenoIndex = numberOfTestTypes * numberOfRpanelsTypes - 1;
		int phenoIndex = 0;

		int ttIndex = 0;
		int rpIndex = 0;
		String ttName = parsingArgs.getTestTypeName(ttIndex);
		String rpName = rpanelTypes.get(rpIndex);

		topHitsFile = resultsFilesInfo.getTopHitsFile(ttIndex, rpIndex);
		// filteredByAllFile = mergeFilesInfo.getFinalFilteredByAllFile(ttIndex,
		// rpIndex);
		// phenomeFileA = phenomeAnalysisFilesInfo.getPhenotypeFile(phenoIndex);
		phenomeFileA = phenomeAnalysisFilesInfo.getPhenotypeIntermediateFile(phenoIndex);
		phenoIndex++;

		doInitPhenoMatrix(parsingArgs, topHitsFile, ttName, rpName, phenomeFileA);

		for (ttIndex = 0; ttIndex < numberOfTestTypes; ttIndex++) {
			int startRp = 0;
			if (ttIndex == 0) {
				startRp = 1;
			}
			ttName = parsingArgs.getTestTypeName(ttIndex);
			for (rpIndex = startRp; rpIndex < numberOfRpanelsTypes; rpIndex++) {
				rpName = rpanelTypes.get(rpIndex);

				topHitsFile = resultsFilesInfo.getTopHitsFile(ttIndex, rpIndex);
				// filteredByAllFile = mergeFilesInfo.getFinalFilteredByAllFile(ttIndex,
				// rpIndex);
				// condensedFile = mergeFilesInfo.getFinalCondensedFile(ttIndex, rpIndex);

				// phenomeFileB = phenomeAnalysisFilesInfo.getPhenotypeFile(phenoIndex);
				phenomeFileB = phenomeAnalysisFilesInfo.getPhenotypeIntermediateFile(phenoIndex);

				doAddToPhenoMatrix(parsingArgs, phenomeFileA, topHitsFile, ttName, rpName, phenomeFileB);

				phenomeFileA = phenomeFileB;
				phenoIndex++;
			}
		}

		phenoIndex = 0;
		// Lets do the fillinout of the phenomeAnalysis.
		for (ttIndex = 0; ttIndex < numberOfTestTypes; ttIndex++) {
			ttName = parsingArgs.getTestTypeName(ttIndex);

			for (rpIndex = 0; rpIndex < numberOfRpanelsTypes; rpIndex++) {
				rpName = rpanelTypes.get(rpIndex);

				// condensedFile = mergeFilesInfo.getFinalCondensedFile(ttIndex, rpIndex);
				filteredByAllFile = mergeFilesInfo.getFinalFilteredByAllFile(ttIndex, rpIndex);
				if (endChrS.equals("23")) {
					filteredByAllXMalesFile = mergeFilesInfo.getAdditionalFilteredByAllXMalesFile(ttIndex, rpIndex, 0);
					filteredByAllXFemalesFile = mergeFilesInfo.getAdditionalFilteredByAllXFemalesFile(ttIndex, rpIndex,
							0);
				} else {
					filteredByAllXMalesFile = mergeFilesInfo.getFinalFilteredByAllFile(ttIndex, rpIndex);
					filteredByAllXFemalesFile = mergeFilesInfo.getFinalFilteredByAllFile(ttIndex, rpIndex);
				}

				phenomeFileB = phenomeAnalysisFilesInfo.getPhenotypeFile(phenoIndex);

				doFilloutPhenoMatrix(parsingArgs, phenomeFileA, filteredByAllFile, filteredByAllXMalesFile,
						filteredByAllXFemalesFile, endChrS, ttName, rpName, phenomeFileB);

				// phenomeFileA = phenomeFileB;
				phenoIndex++;
			}
		}

		// Last round to generate final results
		phenoIndex = 0;
		phenomeFileA = phenomeAnalysisFilesInfo.getPhenotypeFile(phenoIndex);

		phenoIndex++;

		// Lets do the finalization of the phenomeAnalysis.
		for (ttIndex = 0; ttIndex < numberOfTestTypes; ttIndex++) {
			int startRp = 0;
			if (ttIndex == 0) {
				startRp = 1;
			}

			ttName = parsingArgs.getTestTypeName(ttIndex);
			for (rpIndex = startRp; rpIndex < numberOfRpanelsTypes; rpIndex++) {
				rpName = rpanelTypes.get(rpIndex);
				// LOGGER.info("i\n\n\t[makeImputationPerChunk]: ttName " + ttName + " | rpName
				// " + rpName + " |
				// phenoIndex " + phenoIndex);
				// LOGGER.info("\t[makeImputationPerChunk]: phenomeFileA " + phenomeFileA);
				phenomeFileB = phenomeAnalysisFilesInfo.getPhenotypeFile(phenoIndex);
				// LOGGER.info("\t[makeImputationPerChunk]: phenomeFileB " + phenomeFileB);

				phenomeFileC = phenomeAnalysisFilesInfo.getPhenotypeFinalFile(phenoIndex);
				// LOGGER.info("\t[makeImputationPerChunk]: phenomeFileC " + phenomeFileC);

				doFinalizePhenoMatrix(parsingArgs, phenomeFileA, phenomeFileB, ttName, rpName, phenomeFileC);

				phenomeFileA = phenomeFileC;
				phenoIndex++;
			}
		}

	}

	/**
	 * Method the wraps the execution of createSplitedFiles tasks and store the
	 * command in the listOfCommands
	 */

	private static void doCreateSplitedFiles(ParseCmdLine parsingArgs, String gmapFile, String gmapFileChr,
			String theChromo) {

		String cmdToStore = null;

		// if(parsingArgs.getStageStatus("imputeWithImpute")==1) {

		cmdToStore = JAVA_HOME + "/java createSplitedFiles.jar " + gmapFile + " " + gmapFileChr + " " + theChromo;
		listOfCommands.add(cmdToStore);

		try {
			GuidanceImpl.createSplitedFiles(gmapFile, gmapFileChr, theChromo, cmdToStore);
			System.out.println("Executing command " + cmdToStore);
		} catch (Exception e) {
			System.err.println("[Guidance] Exception trying the execution of createSplitedFiles task");
			System.err.println(e.getMessage());
		}
		// }
	}

	/**
	 * Method the wraps the execution of convertFromBedToBed tasks and store the
	 * command in the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param bedFile
	 * @param bimFile
	 * @param famFile
	 * @param mixedBedFile
	 * @param mixedBimFile
	 * @param mixedFamFile
	 * @param mixedBedToBedLogFile
	 * @param theChromo
	 */
	private static void doConvertFromBedToBed(ParseCmdLine parsingArgs, String bedFile, String bimFile, String famFile,
			String mixedBedFile, String mixedBimFile, String mixedFamFile, String mixedBedToBedLogFile,
			String theChromo) {

		if (parsingArgs.getStageStatus("convertFromBedToBed") == 1) {
			// Task
			String cmdToStore = JAVA_HOME + "/java convertFromBedToBed.jar " + bedFile + " " + bimFile + " " + famFile
					+ " " + mixedBedFile + " " + mixedBimFile + " " + mixedFamFile + " " + mixedBedToBedLogFile + " "
					+ theChromo;
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.convertFromBedToBed(bedFile, bimFile, famFile, mixedBedFile, mixedBimFile, mixedFamFile,
						mixedBedToBedLogFile, theChromo, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of convertFromBedToBed task", gte);
			}
		}
	}

	/**
	 * Method the wraps the execution of doSplitChr23 tasks and store the command in
	 * the listOfCommands
	 */

	private static void doSplitChr23(ParseCmdLine parsingArgs, String bedFile, String bimFile, String famFile,
			String logFile, String bedChr23File, String bimChr23File, String famChr23File, String sex,
			String theChromo) {

		String cmdToStore = null;
		String myPrefix = null;

		if (parsingArgs.getStageStatus("convertFromBedToBed") == 1) {

			if (sex.equals(SEX1)) {
				myPrefix = bedChr23File.split("\\.")[0];

				cmdToStore = PLINKBINARY + " --noweb --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
						+ " --filter-males --out " + myPrefix + " --make-bed";
			} else if (sex.equals(SEX2)) {
				myPrefix = bedChr23File.split("\\.")[0];

				cmdToStore = PLINKBINARY + " --noweb --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
						+ " --filter-females --out " + myPrefix + " --make-bed";
			}

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.splitChr23(myPrefix, bedFile, bimFile, famFile, bedChr23File, bimChr23File, famChr23File,
						logFile, sex, theChromo, cmdToStore);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of splitChr23 task");
				System.err.println(e.getMessage());
			}
		}

	}

	/**
	 * Method the wraps the execution of createRsIdList task and store the command
	 * in the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param mixedBimOrGenFile
	 * @param exclCgatFlag
	 * @param mixedPairsFile
	 * @param inputFormat
	 */
	private static void doCreateRsIdList(ParseCmdLine parsingArgs, String mixedBimOrGenFile, String exclCgatFlag,
			String mixedPairsFile, String inputFormat) {

		if (parsingArgs.getStageStatus("createRsIdList") == 1) {
			// Task
			String cmdToStore = JAVA_HOME + "/java createRsIdList " + mixedBimOrGenFile + " " + exclCgatFlag + " "
					+ mixedPairsFile + " " + inputFormat;
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.createRsIdList(mixedBimOrGenFile, exclCgatFlag, mixedPairsFile, inputFormat, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of createRsIdList task", gte);
			}
		}
	}

	/**
	 * Method that wraps the execution of phasing task with bed input formats and
	 * store the command in the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param theChromo
	 * @param bedFile
	 * @param bimFile
	 * @param famFile
	 * @param gmapFile
	 * @param excludedSnpsFile
	 * @param shapeitHapsFile
	 * @param shapeitSampleFile
	 * @param shapeitLogFile
	 * @param filteredHaplotypesFile
	 * @param filteredHaplotypesSampleFile
	 * @param filteredHaplotypesLogFile
	 * @param filteredHaplotypesVcfFile
	 * @param listOfSnpsFile
	 * @param exclCgatFlag
	 * @param exclSVFlag
	 */
	private static void doPhasingBed(ParseCmdLine parsingArgs, String theChromo, String bedFile, String bimFile,
			String famFile, String gmapFile, String sampleFile, String covariables, String responseVar,
			String excludedSnpsFile, String phasingHapsFile, String phasingSampleFile, String phasingLogFile,
			String filteredHaplotypesLogFile, String filteredHaplotypesVcfFile, String filteredHaplotypesVcfFileBgzip,
			String filteredHaplotypesVcfFileBgzipIndexed, String exclCgatFlag, String sex) {

		String phasingTool = parsingArgs.getPhasingTool();
		String myPrefix = phasingHapsFile.split("\\.")[0];
		String cmd;

		if (parsingArgs.getStageStatus("phasingBed") == 1) {

			if (phasingTool.equals("shapeit")) {

				// If we process chromoso X (23) then we change the cmdToStore
				if (theChromo.equals("23")) {
					cmd = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map "
							+ gmapFile + " --chrX --output-max " + phasingHapsFile + " " + phasingSampleFile
							+ " --thread 16 --effective-size 20000 --output-log " + phasingLogFile;

					if (sex.equals(SEX1)) {
						cmd = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile
								+ " --input-map " + gmapFile + " --chrX --output-max " + phasingHapsFile + " "
								+ phasingSampleFile + " --thread 16 --effective-size 20000 --output-log "
								+ phasingLogFile;
					} else if (sex.equals(SEX2)) {
						cmd = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile
								+ " --input-map " + gmapFile + " --chrX --output-max " + phasingHapsFile + " "
								+ phasingSampleFile + " --thread 16 --effective-size 20000 --output-log "
								+ phasingLogFile;
					}

				} else {
					cmd = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map "
							+ gmapFile + " --output-max " + phasingHapsFile + " " + phasingSampleFile
							+ " --thread 16 --effective-size 20000 --output-log " + phasingLogFile;
				}

				listOfCommands.add(new String(cmd));
				try {
					GuidanceImpl.phasingBed(theChromo, bedFile, bimFile, famFile, gmapFile, phasingHapsFile,
							phasingSampleFile, phasingLogFile, myPrefix, sex, phasingTool, cmd);
				} catch (Exception e) {
					System.err.println("[Guidance] Exception trying the execution of phasing task");
					System.err.println(e.getMessage());
				}

			} else if (phasingTool.equals("eagle")) {

				if (parsingArgs.getStageStatus("phasingBed") == 1) {
					// If we process chromoso X (23) then we change the cmdToStore
					if (theChromo.equals("23")) {
						cmd = EAGLEBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
								+ " --chrom " + theChromo + " --geneticMapFile " + gmapFile
								+ " --numThreads 16 --outPrefix " + myPrefix;
						if (sex.equals(SEX1)) {
							cmd = EAGLEBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
									+ " --chrom " + theChromo + " --geneticMapFile " + gmapFile
									+ " --numThreads 16 --outPrefix " + myPrefix;
						} else if (sex.equals(SEX2)) {
							cmd = EAGLEBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
									+ " --chrom " + theChromo + " --geneticMapFile " + gmapFile
									+ " --numThreads 16 --outPrefix " + myPrefix;
						}
					} else {
						cmd = EAGLEBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
								+ " --chrom " + theChromo + " --geneticMapFile " + gmapFile
								+ " --numThreads 16 --outPrefix " + myPrefix;
					}

					listOfCommands.add(cmd);
					try {
						GuidanceImpl.phasingBed(theChromo, bedFile, bimFile, famFile, gmapFile, phasingHapsFile,
								phasingSampleFile, phasingLogFile, myPrefix, sex, phasingTool, cmd);
					} catch (Exception e) {
						System.err.println("[Guidance] Exception trying the execution of phasing task");
						System.err.println(e.getMessage());
					}
				}
			} else {
				System.err.println("[Guidance] Only Eagle and Shapeit are available for phasing haplotypes");

			}

			cmd = JAVA_HOME + " newSample.jar " + sampleFile + " " + phasingSampleFile + " " + covariables + " "
					+ responseVar;

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.newSample(sampleFile, phasingSampleFile, responseVar, covariables, cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of phasing task");
				System.err.println(e.getMessage());
			}

		}

		if (parsingArgs.getStageStatus("createListOfExcludedSnps") == 1) {
			cmd = JAVA_HOME + "/java createListOfExcludedSnps.jar " + phasingHapsFile + " " + excludedSnpsFile + " "
					+ exclCgatFlag;

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.createListOfExcludedSnps(phasingHapsFile, excludedSnpsFile, exclCgatFlag, cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of createListOfExcludedSnps task");
				System.err.println(e.getMessage());
			}
		}

		if (parsingArgs.getStageStatus("filterHaplotypes") == 1) {

			cmd = SHAPEIT_BINARY + " -convert --input-haps " + phasingHapsFile + " " + phasingSampleFile
					+ " --exclude-snp " + excludedSnpsFile + " " + " --output-log " + filteredHaplotypesLogFile
					+ " --output-vcf " + filteredHaplotypesVcfFile;

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.filterHaplotypes(phasingHapsFile, phasingSampleFile, excludedSnpsFile,
						filteredHaplotypesLogFile, filteredHaplotypesVcfFile, cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of filterHaplotypes task");
				System.err.println(e.getMessage());
			}

			cmd = SAMTOOLSBINARY + "/bgzip " + filteredHaplotypesVcfFile + " " + filteredHaplotypesVcfFile + ".gz";

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.samtoolsBgzip(filteredHaplotypesVcfFile, filteredHaplotypesVcfFile + ".gz", cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of samtoolsBgzip task");
				System.err.println(e.getMessage());
			}
		}

	}

	/**
	 * Method that wraps the execution of phasing task with gen input formats and
	 * store the command in the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param theChromo
	 * @param genFile
	 * @param sampleFile
	 * @param gmapFile
	 * @param excludedSnpsFile
	 * @param shapeitHapsFile
	 * @param shapeitSampleFile
	 * @param shapeitLogFile
	 * @param filteredHaplotypesFile
	 * @param filteredHaplotypesSampleFile
	 * @param filteredHaplotypesLogFile
	 * @param filteredHaplotypesVcfFile
	 * @param listOfSnpsFile
	 * @param exclCgatFlag
	 * @param exclSVFlag
	 */
	private static void doPhasing(ParseCmdLine parsingArgs, String theChromo, String genFile, String sampleFile,
			String gmapFile, String excludedSnpsFile, String phasingHapsFile, String phasingSampleFile,
			String phasingLogFile, String filteredHaplotypesLogFile, String filteredHaplotypesVcfFile,
			String filteredHaplotypesVcfFileBgzip, String filteredHaplotypesVcfFileBgzipIndexed, String exclCgatFlag) {

		String phasingTool = parsingArgs.getPhasingTool();
		// String myPrefix = phasingHapsFile.split("\\.")[0];
		String cmd;

		if (parsingArgs.getStageStatus("pashing") == 1) {
			if (phasingTool.equals("shapeit")) {
				String cmdToStore = null;
				if (theChromo.equals("23")) {
					// If we process chromoso X (23) then we change the cmdToStore
					cmdToStore = SHAPEIT_BINARY + " --input-gen " + genFile + " " + sampleFile + " --input-map "
							+ gmapFile + " --chrX --output-max " + phasingHapsFile + " " + phasingSampleFile
							+ " --thread 16 --effective-size 20000 --output-log " + phasingLogFile;
				} else {
					cmdToStore = SHAPEIT_BINARY + " --input-gen " + genFile + " " + sampleFile + " --input-map "
							+ gmapFile + " --output-max " + phasingHapsFile + " " + phasingSampleFile
							+ " --thread 16 --effective-size 20000 --output-log " + phasingLogFile;
				}
				listOfCommands.add(cmdToStore);
				try {
					GuidanceImpl.phasing(theChromo, genFile, sampleFile, gmapFile, phasingHapsFile, phasingSampleFile,
							phasingLogFile, phasingTool, cmdToStore);

				} catch (Exception e) {
					System.err.println("[Guidance] Exception trying the execution of phasing task");
					System.err.println(e.getMessage());
				}
			} else if (phasingTool.equals("eagle")) {
				System.err.println("[Guidance] Eagle does not take input in GEN format");
			} else {
				System.err.println("[Guidance] Only Eagle and Shapeit are available for phasing haplotypes");
			}
		}

		if (parsingArgs.getStageStatus("createListOfExcludedSnps") == 1) {
			cmd = JAVA_HOME + "/java createListOfExcludedSnps " + phasingHapsFile + " " + excludedSnpsFile + " "
					+ exclCgatFlag;

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.createListOfExcludedSnps(phasingHapsFile, excludedSnpsFile, exclCgatFlag, cmd);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of createListOfExcludedSnps task", gte);
			}
		}

		if (parsingArgs.getStageStatus("filterHaplotypes") == 1) {

			cmd = SHAPEIT_BINARY + " -convert --input-haps " + phasingHapsFile + " " + phasingSampleFile
					+ " --exclude-snp " + excludedSnpsFile + " " + " --output-log " + filteredHaplotypesLogFile
					+ " --output-vcf " + filteredHaplotypesVcfFile;

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.filterHaplotypes(phasingHapsFile, phasingSampleFile, excludedSnpsFile,
						filteredHaplotypesLogFile, filteredHaplotypesVcfFile, cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of filterHaplotypes task");
				System.err.println(e.getMessage());
			}

			cmd = SAMTOOLSBINARY + "/bgzip " + filteredHaplotypesVcfFile + " " + filteredHaplotypesVcfFile + ".gz";

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.samtoolsBgzip(filteredHaplotypesVcfFile, filteredHaplotypesVcfFile + ".gz", cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of samtoolsBgzip task");
				System.err.println(e.getMessage());
			}
		}

	}

	/**
	 * Method that wraps the execution of impute task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param chrS
	 * @param gmapFile
	 * @param knownHapFile
	 * @param legendFile
	 * @param shapeitHapsFile
	 * @param shapeitSampleFile
	 * @param lim1S
	 * @param lim2S
	 * @param pairsFile
	 * @param imputeFile
	 * @param imputeFileInfo
	 * @param imputeFileSummary
	 * @param imputeFileWarnings
	 */
	private static void doImputationWithImpute(ParseCmdLine parsingArgs, String chrS, String gmapFile,
			String knownHapFile, String legendFile, String phasingHapsFile, String phasingSampleFile, String lim1S,
			String lim2S, String pairsFile, String imputeFile, String imputeFileInfo, String imputeFileSummary,
			String imputeFileWarnings, String sex) {

		String cmdToStore = null;
		if (parsingArgs.getStageStatus("imputeWithImpute") == 1) {
			// Submitting the impute task per chunk

			if (chrS.equals("23")) {
				cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
						+ legendFile + " -known_haps_g " + phasingHapsFile + " -sample_g " + phasingSampleFile
						+ " -int " + lim1S + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile
						+ " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r "
						+ imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz ";

				if (sex.equals("males")) {
					cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
							+ legendFile + " -known_haps_g " + phasingHapsFile + " -sample_g " + phasingSampleFile
							+ " -int " + lim1S + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile
							+ " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r "
							+ imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz ";
				} else if (sex.equals("females")) {
					cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
							+ legendFile + " -known_haps_g " + phasingHapsFile + " -sample_g " + phasingSampleFile
							+ " -int " + lim1S + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile
							+ " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r "
							+ imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz ";
				}

			} else {
				cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
						+ legendFile + " -known_haps_g " + phasingHapsFile + " -int " + lim1S + " " + lim2S
						+ " -exclude_snps_g " + pairsFile + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i "
						+ imputeFileInfo + " -r " + imputeFileSummary + " -w " + imputeFileWarnings
						+ " -no_sample_qc_info -o_gz";
			}
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.imputeWithImpute(gmapFile, knownHapFile, legendFile, phasingHapsFile, phasingSampleFile,
						lim1S, lim2S, pairsFile, imputeFile, imputeFileInfo, imputeFileSummary, imputeFileWarnings,
						chrS, sex, cmdToStore);

			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of impute task");
				System.err.println(e.getMessage());
			}
		}

	}

	/**
	 * Method that wraps the execution of Imputation task with Minimac2 and store
	 * the command in the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param knownHapFile
	 * @param filteredHapsFile
	 * @param filteredSampleFile
	 * @param filteredListOfSnpsFile
	 * @param imputedMMFileName
	 * @param imputedMMInfoFile
	 * @param imputedMMErateFile
	 * @param imputedMMRecFile
	 * @param imputedMMDoseFile
	 * @param imputedMMLogFile
	 * @param chrS
	 * @param lim1S
	 * @param lim2S
	 */
	private static void doImputationWithMinimac(ParseCmdLine parsingArgs, String refVcfFile,
			String filteredHaplotypesVcfFileBgzip, String chrS, String lim1S, String lim2S, String imputeFile,
			String imputeFileInfo, String imputeFileErate, String imputeFileRec, String imputeFileM3vcf,
			String imputeFileLog, String imputeFileBgzip, String imputeMMM3VCFFileBgzip, String imputeFileTbi,
			String sex) {

		if (parsingArgs.getStageStatus("imputeWithMinimac") == 1) {
			// Submitting the impute task per chunk
			// We don't distinguish chrS 23 since the cmdToStore is the same
			String cmdToStore = null;
			String myPrefix = imputeFile.split("\\.")[0];

			if (chrS.equals("23")) {
				if (sex.equals(SEX1)) {
					cmdToStore = MINIMAC_BINARY + " --refHaps " + refVcfFile + " --haps "
							+ filteredHaplotypesVcfFileBgzip + " --start " + lim1S + " --end " + lim2S
							+ " --chr X  --window 500000 --prefix " + myPrefix
							+ " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
				} else if (sex.equals(SEX2)) {
					cmdToStore = MINIMAC_BINARY + " --refHaps " + refVcfFile + " --haps "
							+ filteredHaplotypesVcfFileBgzip + " --start " + lim1S + " --end " + lim2S
							+ " --chr X --window 500000 --prefix " + myPrefix
							+ " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
				}
			} else {

				cmdToStore = MINIMAC_BINARY + " --refHaps " + refVcfFile + " --haps " + filteredHaplotypesVcfFileBgzip
						+ " --start " + lim1S + " --end " + lim2S + " --chr " + chrS + " --window 500000 --prefix "
						+ myPrefix + " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
			}

			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.imputeWithMinimac(refVcfFile, filteredHaplotypesVcfFileBgzip, imputeFile, imputeFileInfo,
						imputeFileErate, imputeFileRec, imputeFileM3vcf, imputeFileLog, chrS, lim1S, lim2S, myPrefix,
						sex, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of imputationWithMinimac task", gte);
			}

			cmdToStore = SAMTOOLSBINARY + "/bgzip " + imputeFile + " " + imputeFile;

			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.samtoolsBgzip(imputeFile, imputeFileBgzip, cmdToStore);
			} catch (Exception e) {
				LOGGER.error("[Guidance] Exception trying the execution of samtoolsBgzip task");
			}

			System.out.println("Handling samtoolsTabix for " + chrS + " from " + lim1S + " to " + lim2S);

			cmdToStore = SAMTOOLSBINARY + "/tabix -p vcf -f " + imputeFileBgzip;

			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.samtoolsTabix(imputeFileBgzip, imputeFileTbi, cmdToStore);
			} catch (Exception e) {
				LOGGER.error("[Guidance] Exception trying the execution of samtoolsTabix task");
			}
			System.out.println("Ended minimac imputation for crhomo " + chrS + " from " + lim1S + " to " + lim2S);
		}
	}

	/**
	 * Method that wraps the execution of filterByInfo task and store the command in
	 * the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param imputeFileInfo
	 * @param filteredRsIdFile
	 * @param infoThresholdS
	 */
	private static void doFilterByInfo(ParseCmdLine parsingArgs, String imputeFileInfo, String filteredRsIdFile,
			String infoThresholdS, String sex) {

		if (parsingArgs.getStageStatus("filterByInfo") == 1) {
			String imputationTool = parsingArgs.getImputationTool();
			String cmdToStore = null;
			// We create the list of rsId that are greater than or equal to the
			// infoThreshold value
			cmdToStore = JAVA_HOME + "/java filterByInfo " + imputeFileInfo + " " + filteredRsIdFile + " "
					+ infoThresholdS;
			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.filterByInfo(imputationTool, imputeFileInfo, filteredRsIdFile, infoThresholdS, cmdToStore);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of filterByInfo task");
				System.err.println(e.getMessage());
			}

		}
	}

	/**
	 * Method that wraps the execution of qctoolS task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param imputeFile
	 * @param filteredRsIdFile
	 * @param filteredFile
	 * @param filteredLogFile
	 */
	private static void doQctoolS(ParseCmdLine parsingArgs, String imputeFile, String filteredRsIdFile,
			String filteredFile, String filteredLogFile, String sex) {

		double mafThreshold = parsingArgs.getMafThreshold();
		String mafThresholdS = Double.toString(mafThreshold);
		String imputationTool = parsingArgs.getImputationTool();

		if (parsingArgs.getStageStatus("qctoolS") == 1) {
			String cmdToStore = null;

			if (imputationTool.equals("impute")) {

				cmdToStore = QCTOOL_BINARY + " -g " + imputeFile + " -og " + filteredFile + " -incl-rsids "
						+ filteredRsIdFile + " -omit-chromosome -force -log " + filteredLogFile + " -maf "
						+ mafThresholdS + " 1";

			} else if (imputationTool.equals("minimac")) {

				cmdToStore = QCTOOL_BINARY + " -g " + imputeFile + " -og " + filteredFile + " -incl-rsids "
						+ filteredRsIdFile + " -omit-chromosome -force -log " + filteredLogFile + " -maf "
						+ mafThresholdS + " 1 " + "-vcf-genotype-field GP";
			}

			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.qctoolS(imputationTool, imputeFile, filteredRsIdFile, mafThresholdS, filteredFile,
						filteredLogFile, cmdToStore);
			} catch (Exception gte) {
				LOGGER.error("[Guidance] Exception trying the execution of qctoolS tasks for controls", gte);
			}

		}
	}

	/**
	 * Method that wraps the execution of snptest task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param chrS
	 * @param mergedGenFile
	 * @param mergedSampleFile
	 * @param snptestOutFile
	 * @param snptestLogFile
	 * @param responseVar
	 * @param covariables
	 */
	private static void doSnptest(ParseCmdLine parsingArgs, String chrS, String mergedGenFile, String mergedSampleFile,
			String snptestOutFile, String snptestLogFile, String responseVar, String covariables) {

		String newStr = covariables.replace(',', ' ');
		if (parsingArgs.getStageStatus("snptest") == 1) {
			String cmdToStore;
			if (covariables.equals("none")) {
				cmdToStore = SNPTEST_BINARY + " -data " + mergedGenFile + " " + mergedSampleFile + " -o "
						+ snptestOutFile + " -pheno " + responseVar + " -hwe -log " + snptestLogFile;
			} else {
				cmdToStore = SNPTEST_BINARY + " -data " + mergedGenFile + " " + mergedSampleFile + " -o "
						+ snptestOutFile + " -pheno " + responseVar + " -cov_names " + newStr + " -hwe -log "
						+ snptestLogFile;
			}

			// Different parameters for chromo 23 (X) and the rest.
			if (chrS.equals("23")) {
				cmdToStore = cmdToStore + " -method newml -assume_chromosome X -stratify_on sex -frequentist 1";
			} else {
				cmdToStore = cmdToStore + " -method em -frequentist 1 2 3 4 5";
			}

			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}

			// Submitting the snptest task per this chunk
			try {
				GuidanceImpl.snptest(mergedGenFile, mergedSampleFile, snptestOutFile, snptestLogFile, responseVar,
						covariables, chrS, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of snptest task", gte);
			}
		}
	}

	/**
	 * Method that wraps the execution of collectSummary task and store the command
	 * in the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param chrS
	 * @param imputeFileInfo
	 * @param snptestOutFile
	 * @param summaryFile
	 * @param mafThresholdS
	 * @param infoThresholdS
	 * @param hweCohortThresholdS
	 * @param hweCasesThresholdS
	 * @param hweControlsThresholdS
	 */
	private static void doCollectSummary(ParseCmdLine parsingArgs, String chrS, String imputeFileInfo,
			String snptestOutFile, String summaryFile, String mafThresholdS, String infoThresholdS,
			String hweCohortThresholdS, String hweCasesThresholdS, String hweControlsThresholdS) {

		if (parsingArgs.getStageStatus("collectSummary") == 1) {
			// Submitting the collect_summary task per this chunk
			String cmdToStore = JAVA_HOME + "/java collectSummary " + chrS + " " + imputeFileInfo + " " + snptestOutFile
					+ " " + summaryFile + " " + mafThresholdS + " " + infoThresholdS + " " + hweCohortThresholdS + " "
					+ hweCasesThresholdS + " " + hweControlsThresholdS;

			listOfCommands.add(cmdToStore);

			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}

			try {
				GuidanceImpl.collectSummary(chrS, imputeFileInfo, snptestOutFile, summaryFile, mafThresholdS,
						infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of collectSummary task", gte);
			}
		}
	}

	/**
	 * Method that wraps the jointCondensedFiles task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param condensedA
	 * @param condensedB
	 * @param condensedC
	 */
	private static void doJointCondenseFiles(ParseCmdLine parsingArgs, String condensedA, String condensedB,
			String condensedC) {

		if (parsingArgs.getStageStatus("jointCondensedFiles") == 1) {
			String cmdToStore = JAVA_HOME + "/java jointCondensedFiles " + condensedA + " " + condensedB + " "
					+ condensedC;
			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.jointCondensedFiles(condensedA, condensedB, condensedC, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of jointCondensedFiles task", gte);
			}
		}
	}

	/**
	 * Method that wraps the doJointFilteredByAllFiles task and store the command in
	 * the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param filteredByAllA
	 * @param filteredByAllB
	 * @param filteredByAllC
	 * @param rpanelName
	 * @param rpanelFlag
	 */
	private static void doJointFilteredByAllFiles(ParseCmdLine parsingArgs, String filteredByAllA,
			String filteredByAllB, String filteredByAllC, String rpanelName, String rpanelFlag) {

		if (parsingArgs.getStageStatus("jointFilteredByAllFiles") == 1) {
			String cmdToStore = JAVA_HOME + "/java jointFilteredByAllFiles " + filteredByAllA + " " + filteredByAllB
					+ " " + filteredByAllC + " " + rpanelName + " " + rpanelFlag;
			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.jointFilteredByAllFiles(filteredByAllA, filteredByAllB, filteredByAllC, rpanelName,
						rpanelFlag, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of getBestSnps task", gte);
			}
		}
	}

	/**
	 * Method that wraps the generateTopHits task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param filteredFile
	 * @param filteredXFile
	 * @param topHitsResults
	 * @param pvaThrS
	 */
	private static void doGenerateTopHits(ParseCmdLine parsingArgs, String filteredFile, String filteredXMalesFile,
			String filteredXFemalesFile, String topHitsResults, String pvaThrS) {

		if (parsingArgs.getStageStatus("generateTopHits") == 1) {
			String cmdToStore = JAVA_HOME + "/java generateTopHits " + filteredFile + " " + filteredXMalesFile + " "
					+ filteredXFemalesFile + " " + topHitsResults + " " + pvaThrS;
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.generateTopHitsAll(filteredFile, filteredXMalesFile, filteredXFemalesFile, topHitsResults,
						pvaThrS, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of generateTopHits task", gte);
			}
		}
	}

	/**
	 * Method that wraps the generateQQManhattanPlots task and store the command in
	 * the listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param condensedFile
	 * @param qqPlotFile
	 * @param manhattanPlotFile
	 * @param qqPlotTiffFile
	 * @param manhattanPlotTiffFile
	 * @param correctedPvaluesFile
	 */
	private static void doGenerateQQManhattanPlots(ParseCmdLine parsingArgs, String condensedFile, String qqPlotFile,
			String manhattanPlotFile, String qqPlotTiffFile, String manhattanPlotTiffFile) {

		if (parsingArgs.getStageStatus("generateQQManhattanPlots") == 1) {
			String cmdToStore = R_SCRIPT_BIN_DIR + "/Rscript " + R_SCRIPT_DIR + "/qqplot_manhattan.R " + condensedFile
					+ " " + qqPlotFile + " " + manhattanPlotFile + " " + qqPlotTiffFile + " " + manhattanPlotTiffFile;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.generateQQManhattanPlots(condensedFile, qqPlotFile, manhattanPlotFile, qqPlotTiffFile,
						manhattanPlotTiffFile, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of generateQQManhattanPlots task", gte);
			}
		}
	}

	/**
	 * Method that wraps the combinePanelsComplex task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param resultsPanelA
	 * @param resultsPanelB
	 * @param lastResultFile
	 * @param lim1
	 * @param lim2
	 */
	private static void doCombinePanelsComplex(ParseCmdLine parsingArgs, String resultsPanelA, String resultsPanelB,
			String resultsPanelC, int lim1, int lim2) {

		String cmdToStore = JAVA_HOME + "/java combinePanelsComplex " + resultsPanelA + " " + resultsPanelB + " "
				+ resultsPanelC + " " + lim1 + " " + lim2;
		listOfCommands.add(cmdToStore);

		try {
			GuidanceImpl.combinePanelsComplex(resultsPanelA, resultsPanelB, resultsPanelC, lim1, lim2, cmdToStore);
		} catch (GuidanceTaskException gte) {
			LOGGER.error("[Guidance] Exception trying the execution of combinePanelsComplex task", gte);
		}

	}

	/**
	 * Method that wraps the mergeTwoChunks task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param reduceA
	 * @param reduceB
	 * @param reduceC
	 * @param theChromo
	 * @param type
	 */
	private static void doMergeTwoChunks(ParseCmdLine parsingArgs, String reduceA, String reduceB, String reduceC) {

		if (parsingArgs.getStageStatus("mergeTwoChunks") == 1) {
			// Task
			String cmdToStore = JAVA_HOME + "/java mergeTwoChunks " + reduceA + " " + reduceB + " " + reduceC;
			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.mergeTwoChunks(reduceA, reduceB, reduceC, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of mergeTwoChunks task", gte);
			}

		}
	}

	/**
	 * Method that wraps the mergeTwoChunks task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param reduceA
	 * @param reduceB
	 * @param reduceC
	 * @param theChromo
	 * @param type
	 */
	private static void doMergeTwoChunksUnconditional(ParseCmdLine parsingArgs, String reduceA, String reduceB,
			String reduceC) {

		String cmdToStore = JAVA_HOME + "/java mergeTwoChunksUnconditional " + reduceA + " " + reduceB + " " + reduceC;
		listOfCommands.add(cmdToStore);
		try {
			flushCommands();
		} catch (IOException e) {
			LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
		}
		try {
			GuidanceImpl.mergeTwoChunks(reduceA, reduceB, reduceC, cmdToStore);
		} catch (GuidanceTaskException gte) {
			LOGGER.error("[Guidance] Exception trying the execution of mergeTwoChunks task", gte);
		}

	}

	/**
	 * Method that wraps the mergeTwoChunks task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param reduceA
	 * @param reduceB
	 * @param reduceC
	 * @param theChromo
	 * @param type
	 */
	/*
	 * private static void doMergeTwoChunksInTheFirst(ParseCmdLine parsingArgs,
	 * String reduceA, String reduceB, String theChromo, String type) {
	 * 
	 * if (parsingArgs.getStageStatus("mergeTwoChunks") == 1) { // Task String
	 * cmdToStore = JAVA_HOME + "/java mergeTwoChunksInTheFirst " + reduceA + " " +
	 * reduceB + " " + theChromo; listOfCommands.add(cmdToStore);
	 * 
	 * try { GuidanceImpl.mergeTwoChunksInTheFirst(reduceA, reduceB, cmdToStore); }
	 * catch (GuidanceTaskException gte) { LOGGER.
	 * error("[Guidance] Exception trying the execution of mergeTwoChunks task",
	 * gte); }
	 * 
	 * } }
	 */

	/**
	 * Method that wraps the filterByAll task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param inputFile
	 * @param outputFile
	 * @param outputCondensedFile
	 */
	private static void doFilterByAll(ParseCmdLine parsingArgs, String inputFile, String outputFile,
			String outputCondensedFile, String sex, String rpanelName) {

		String imputationTool = parsingArgs.getImputationTool();
		double mafThreshold = parsingArgs.getMafThreshold();
		double infoThreshold = parsingArgs.getInfoThreshold();
		double hweCohortThreshold = parsingArgs.getHweCohortThreshold();
		double hweCasesThreshold = parsingArgs.getHweCasesThreshold();
		double hweControlsThreshold = parsingArgs.getHweControlsThreshold();

		String mafThresholdS = Double.toString(mafThreshold);
		String infoThresholdS = Double.toString(infoThreshold);
		String hweCohortThresholdS = Double.toString(hweCohortThreshold);
		String hweCasesThresholdS = Double.toString(hweCasesThreshold);
		String hweControlsThresholdS = Double.toString(hweControlsThreshold);

		// Task
		if (parsingArgs.getStageStatus("filterByAll") == 1) {
			String cmdToStore = JAVA_HOME + "/java filterByAll " + imputationTool + " " + inputFile + " " + outputFile + " "
					+ outputCondensedFile + " " + mafThresholdS + " " + infoThresholdS + " " + hweCohortThresholdS + " "
					+ hweCasesThresholdS + " " + hweControlsThresholdS + " " + sex + " " + rpanelName;

			listOfCommands.add(cmdToStore);
			try {
				flushCommands();
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}

			try {
				GuidanceImpl.filterByAll(imputationTool, inputFile, outputFile, outputCondensedFile, mafThresholdS,
						infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, sex, rpanelName,
						cmdToStore);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of filterByAll task");
				System.err.println(e.getMessage());
			}
		}

	}

	/**
	 * Method that wraps the initPhenoMatrix task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param topHitsFile
	 * @param ttName
	 * @param rpName
	 * @param phenomeFile
	 */
	private static void doInitPhenoMatrix(ParseCmdLine parsingArgs, String topHitsFile, String ttName, String rpName,
			String phenomeFile) {

		if (parsingArgs.getStageStatus("initPhenoMatrix") == 1) {
			String cmdToStore = JAVA_HOME + "/java initPhenoMatrix " + topHitsFile + " " + ttName + " " + rpName + " "
					+ phenomeFile;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.initPhenoMatrix(topHitsFile, ttName, rpName, phenomeFile, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of initPhenoMatrix task", gte);
			}
		}
	}

	/**
	 * Method that wraps the addToPhenoMatrix task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param phenomeFileA
	 * @param topHitsFile
	 * @param ttName
	 * @param rpName
	 * @param phenomeFileB
	 */
	private static void doAddToPhenoMatrix(ParseCmdLine parsingArgs, String phenomeFileA, String topHitsFile,
			String ttName, String rpName, String phenomeFileB) {

		if (parsingArgs.getStageStatus("addToPhenoMatrix") == 1) {
			String cmdToStore = JAVA_HOME + "/java addToPhenoMatrix " + phenomeFileA + " " + topHitsFile + " " + ttName
					+ " " + rpName + " " + phenomeFileB;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.addToPhenoMatrix(phenomeFileA, topHitsFile, ttName, rpName, phenomeFileB, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of addToPhenoMatrix task", gte);
			}
		}
	}

	/**
	 * Method that wraps the filloutPhenoMatrix task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param phenomeFileA
	 * @param filteredByAllFile
	 * @param filteredByAllXFile
	 * @param endChrS
	 * @param ttName
	 * @param rpName
	 * @param phenomeFileB
	 */
	private static void doFilloutPhenoMatrix(ParseCmdLine parsingArgs, String phenomeFileA, String filteredByAllFile,
			String filteredByAllXMalesFile, String filteredByAllXFemalesFile, String endChrS, String ttName,
			String rpName, String phenomeFileB) {

		if (parsingArgs.getStageStatus("filloutPhenoMatrix") == 1) {
			String cmdToStore = JAVA_HOME + "/java filloutPhenoMatrix " + phenomeFileA + " " + filteredByAllFile + " "
					+ filteredByAllXMalesFile + " " + filteredByAllXFemalesFile + " " + endChrS + " " + ttName + " "
					+ rpName + " " + phenomeFileB;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.filloutPhenoMatrix(phenomeFileA, filteredByAllFile, filteredByAllXMalesFile,
						filteredByAllXFemalesFile, endChrS, ttName, rpName, phenomeFileB, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of filloutPhenoMatrix task", gte);
			}
		}
	}

	/**
	 * Method that wraps the finalizePhenoMatrix task and store the command in the
	 * listOfCommands
	 * 
	 * @param parsingArgs
	 * @param listOfCommands
	 * @param phenomeFileA
	 * @param phenomeFileB
	 * @param ttName
	 * @param rpName
	 * @param phenomeFileC
	 */
	private static void doFinalizePhenoMatrix(ParseCmdLine parsingArgs, String phenomeFileA, String phenomeFileB,
			String ttName, String rpName, String phenomeFileC) {

		if (parsingArgs.getStageStatus("finalizePhenoMatrix") == 1) {
			String cmdToStore = JAVA_HOME + "/java finalizePhenoMatrix " + phenomeFileA + " " + phenomeFileB + " "
					+ ttName + " " + rpName + " " + phenomeFileC;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.finalizePhenoMatrix(phenomeFileA, phenomeFileB, ttName, rpName, phenomeFileC, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of finalizePhenoMatrix task", gte);
			}
		}
	}

	private static void doCopyFile(ParseCmdLine parsingArgs, String originPath, String destinationPath)
			throws IOException, GuidanceTaskException {
		if (parsingArgs.getStageStatus("combinePanelsComplex") == 1) {
			String cmdToStore = JAVA_HOME + "/java copyFile " + originPath + " " + destinationPath;

			listOfCommands.add(cmdToStore);

			GuidanceImpl.copyFile(originPath, destinationPath);
		}
	}

	/**
	 * Method to set the final status of each file that has been generated in the
	 * execution So far, the three main status are: uncompressed, compressed, and
	 * deleted. The default status is: uncompressed. That is the initial status
	 * assigned to each file
	 * 
	 * @param parsingArgs
	 * @param commonFilesInfo
	 */
	private static void setFinalStatusForCommonFiles(ParseCmdLine parsingArgs, CommonFiles commonFilesInfo) {
		String compressFiles = parsingArgs.getCompressFiles();

		// This status finalStatus is used only for files that can be compressed,
		// uncompressed or deleted.
		// There are sevaral files that will be kept as they have been generated
		// (uncompressed)
		String finalStatus = null;
		// if( removeTempFiles.equals("YES") ) {
		// finalStatus = "deleted";
		// } else
		if (compressFiles.equals("NO")) {
			finalStatus = "uncompressed";
		} else {
			finalStatus = "compressed";
		}

		int startChr = parsingArgs.getStart();
		int endChr = parsingArgs.getEnd();
		for (int i = startChr; i <= endChr; i++) {
			commonFilesInfo.setPairsFileFinalStatus(i, finalStatus);
			commonFilesInfo.setPhasingHapsFileFinalStatus(i, finalStatus);
			commonFilesInfo.setPhasingSampleFileFinalStatus(i, finalStatus);
		}
	}

	/**
	 * Method to set the final status of each file that has been generated in the
	 * execution So far, the three main status are: uncompressed, compressed, and
	 * deleted. The default status is: uncompressed. That is the initial status
	 * assigned to each file
	 * 
	 * @param parsingArgs
	 * @param imputationFilesInfo
	 * @param generalChromoInfo
	 * @param refPanels
	 */
	private static void setFinalStatusForImputationFiles(ParseCmdLine parsingArgs, ImputationFiles imputationFilesInfo,
			List<String> refPanels) {

		int startChr = parsingArgs.getStart();
		int endChr = parsingArgs.getEnd();
		int chunkSize = parsingArgs.getChunkSize();
		String imputationTool = parsingArgs.getImputationTool();
		String compressFiles = parsingArgs.getCompressFiles();

		// This status finalStatus is used only for files that can be compressed,
		// uncompressed or deleted.
		// There are sevaral files that will be kept as they have been generated
		// (uncompressed)
		String finalStatus = null;
		// if( removeTempFiles.equals("YES") ) {
		// finalStatus = "deleted";
		// } else
		if (compressFiles.equals("NO")) {
			finalStatus = "uncompressed";
		} else {
			finalStatus = "compressed";
		}

		if (imputationTool.equals("impute")) {
			for (int j = 0; j < refPanels.size(); j++) {
				// String rPanel = refPanels.get(j);
				for (int i = startChr; i <= endChr; i++) {
					int chromo = i;
					int lim1 = 1;
					int lim2 = lim1 + chunkSize - 1;
					int numberOfChunks = ChromoInfo.getMaxSize(chromo) / chunkSize;
					int module = ChromoInfo.getMaxSize(chromo) % chunkSize;
					if (module != 0)
						numberOfChunks++;

					for (int k = 0; k < numberOfChunks; k++) {
						if (chromo == 23) {
							imputationFilesInfo.setImputedMalesFileFinalStatus(j, lim1, lim2, chunkSize, finalStatus);
							imputationFilesInfo.setFilteredMalesFileFinalStatus(j, lim1, lim2, chunkSize, finalStatus);
							imputationFilesInfo.setImputedMalesInfoFileFinalStatus(j, lim1, lim2, chunkSize,
									finalStatus);
							imputationFilesInfo.setImputedFemalesFileFinalStatus(j, lim1, lim2, chunkSize, finalStatus);
							imputationFilesInfo.setFilteredFemalesFileFinalStatus(j, lim1, lim2, chunkSize,
									finalStatus);
							imputationFilesInfo.setImputedFemalesInfoFileFinalStatus(j, lim1, lim2, chunkSize,
									finalStatus);
						} else {
							imputationFilesInfo.setImputedFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
									finalStatus);
							imputationFilesInfo.setFilteredFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
									finalStatus);
							imputationFilesInfo.setImputedInfoFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
									finalStatus);
						}
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}
				}
			}
		} else if (imputationTool.equals("minimac")) {
			for (int j = 0; j < refPanels.size(); j++) {
				// String rPanel = refPanels.get(j);
				for (int i = startChr; i <= endChr; i++) {
					int chromo = i;
					int lim1 = 1;
					int lim2 = lim1 + chunkSize - 1;
					int numberOfChunks = ChromoInfo.getMaxSize(chromo) / chunkSize;
					int module = ChromoInfo.getMaxSize(chromo) % chunkSize;
					if (module != 0)
						numberOfChunks++;

					for (int k = 0; k < numberOfChunks; k++) {
						if (chromo == 23) {
							imputationFilesInfo.setImputedMMInfoMalesFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
									finalStatus);
							imputationFilesInfo.setImputedMMInfoFemalesFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
									finalStatus);
						} else {
							imputationFilesInfo.setImputedMMInfoFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
									finalStatus);
						}
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}
				}
			}
		}

	}

	/**
	 * Method to set the final status of each file that has been generated in the
	 * execution So far, the three main status are: uncompressed, compressed, and
	 * deleted. The default status is: uncompressed. That is the initial status
	 * assigned to each file
	 * 
	 * @param parsingArgs
	 * @param assocFilesInfo
	 * @param generalChromoInfo
	 * @param refPanels
	 */
	private static void setFinalStatusForAssocFiles(ParseCmdLine parsingArgs, AssocFiles assocFilesInfo,
			List<String> refPanels) {

		int startChr = parsingArgs.getStart();
		int endChr = parsingArgs.getEnd();
		int chunkSize = parsingArgs.getChunkSize();

		String compressFiles = parsingArgs.getCompressFiles();

		// This status finalStatus is used only for files that can be compressed,
		// uncompressed or deleted.
		// There are sevaral files that will be kept as they have been generated
		// (uncompressed)
		String finalStatus = null;
		// if( removeTempFiles.equals("YES") ) {
		// finalStatus = "deleted";
		// } else
		if (compressFiles.equals("NO")) {
			finalStatus = "uncompressed";
		} else {
			finalStatus = "compressed";
		}

		// Now we continue with the association
		int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
		for (int tt = 0; tt < numberOfTestTypes; tt++) {
			for (int j = 0; j < refPanels.size(); j++) {
				// String rPanel = refPanels.get(j);
				for (int i = startChr; i <= endChr; i++) {
					int chromo = i;

					int maxSize = ChromoInfo.getMaxSize(chromo);
					int total_chunks = maxSize / chunkSize;
					int lim1 = 1;
					int lim2 = lim1 + chunkSize - 1;

					for (int k = 0; k < total_chunks; k++) {
						if (chromo == 23) {
							assocFilesInfo.setSnptestOutMalesFileFinalStatus(tt, j, lim1, lim2, chunkSize, finalStatus);
							assocFilesInfo.setSummaryMalesFileFinalStatus(tt, j, lim1, lim2, chunkSize, finalStatus);
							assocFilesInfo.setSnptestOutFemalesFileFinalStatus(tt, j, lim1, lim2, chunkSize,
									finalStatus);
							assocFilesInfo.setSummaryFemalesFileFinalStatus(tt, j, lim1, lim2, chunkSize, finalStatus);

						} else {
							assocFilesInfo.setSnptestOutFileFinalStatus(tt, j, chromo, lim1, lim2, chunkSize,
									finalStatus);
							assocFilesInfo.setSummaryFileFinalStatus(tt, j, chromo, lim1, lim2, chunkSize, finalStatus);

						}

						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}
					/*
					 * int index1=0; int index2=1; int
					 * last_index=assocFilesInfo.getTheLastReducedFileIndex(j,chromo);
					 * 
					 * for(int k=1; k<=last_index; k++) { //String theStatus =null; //String
					 * theStatusPrev= null; assocFilesInfo.setReducedFileFinalStatus(j, chromo, k,
					 * finalStatus); }
					 */
				}
			}
		}
	}

	/**
	 * Method to print the current status of each stage (0: unactive, 1:active) of
	 * Guidance workflow.
	 * 
	 * @param parsingArgs
	 */
	private static void printStagesStatus(ParseCmdLine parsingArgs) {
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
		LOGGER.info("[Guidance] combineCondensedFiles    " + parsingArgs.getStageStatus("combineCondensedFiles"));
		LOGGER.info("[Guidance] initPhenoMatrix          " + parsingArgs.getStageStatus("initPhenoMatrix"));
		LOGGER.info("[Guidance] addToPhenoMatrix         " + parsingArgs.getStageStatus("addToPhenoMatrix"));
		LOGGER.info("[Guidance] filloutPhenoMatrix       " + parsingArgs.getStageStatus("filloutPhenoMatrix"));
		LOGGER.info("[Guidance] finalizePhenoMatrix      " + parsingArgs.getStageStatus("finalizePhenoMatrix"));
		LOGGER.info("[Guidance] taskx                    " + parsingArgs.getStageStatus("taskx"));
		LOGGER.info("[Guidance] tasky                    " + parsingArgs.getStageStatus("tasky"));
		LOGGER.info("[Guidance] taskz                    " + parsingArgs.getStageStatus("taskz"));
	}

	/**
	 * Method to print the general information of Guidance
	 * 
	 */
	private static void printGuidancePackageVersion() {
		LOGGER.warn("[Guidance] *****************************************************************");
		LOGGER.warn("[Guidance] ** This is the Guidance framework to performing imputation,    **");
		LOGGER.warn("[Guidance] ** GWAS and Phenotype analysis of large scale GWAS datasets.   **");
		LOGGER.warn("[Guidance] ** Version: Guidance " + GUIDANCE_VERSION + " build " + GUIDANCE_BUILDNUMBER
				+ "           **");
		LOGGER.warn("[Guidance] ** Date release: 20-Jul-2016                                   **");
		LOGGER.warn("[Guidance] ** Contact: http://cg.bsc.es/guidance                          **");
		LOGGER.warn("[Guidance] *****************************************************************\n");
	}

	/**
	 * Method to verify that all the environment variables have been well defined
	 * 
	 */
	private static void verifyEnvVar() throws GuidanceEnvironmentException {
		verify(PLINK_BINARY);
		verify(GTOOL_BINARY);

		verifyEnvVarDefined(R_SCRIPT_BIN_DIR);
		verifyEnvVarDefined(R_SCRIPT_DIR);

		verify(QCTOOL_BINARY);
		verify(SHAPEIT_BINARY);
		verify(IMPUTE2_BINARY);
		verify(SNPTEST_BINARY);

		verifyEnvVarDefined(JAVA_HOME);
	}

	/**
	 * Verifies that the given environment variable is correctly defined
	 * 
	 * @param envVar
	 * @throws GuidanceEnvironmentException
	 */
	private static void verifyEnvVarDefined(String envVar) throws GuidanceEnvironmentException {
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
	private static void verify(String envVar) throws GuidanceEnvironmentException {
		verifyEnvVarDefined(envVar);

		File f = new File(envVar);
		if (!f.exists() || f.isDirectory()) {
			throw new GuidanceEnvironmentException("[Guidance] Error, " + envVar
					+ " does not exist or it is not a binary file. Please check your .bashrc");
		}
	}

}
