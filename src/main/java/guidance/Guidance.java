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
import guidance.utils.ParseCmdLine;
import guidance.utils.GeneralUtils;

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

	// Enable barriers
	private static final boolean BARRIERS = false;

	// Enable flushing to ease the debug
	private static final boolean FLUSH = true;

	// Threshold
	// private static final double PVA_THRESHOLD = 5e-8;
	// private static final String PVA_THRESHOLD_STR =
	// Double.toString(PVA_THRESHOLD);

	// Execution mode variables
	private static final String FILTERED = "filtered";
	private static final String CONDENSED = "condensed";

	// Environment variables that contain the information of the applications that
	// will be used in the Guidance
	// execution
	private static final String PLINK_BINARY = System.getenv("PLINKBINARY");
	private static final String R_SCRIPT_BIN_DIR = System.getenv("RSCRIPTBINDIR");
	private static final String R_SCRIPT_DIR = System.getenv("RSCRIPTDIR");
	private static final String BASHSCRIPTDIR = System.getenv("BASHSCRIPTDIR");
	private static final String QCTOOL_BINARY = System.getenv("QCTOOLBINARY");
	private static final String EAGLEBINARY = System.getenv("EAGLEBINARY");
	private static final String PLINKBINARY = System.getenv("PLINKBINARY");
	private static final String BCFTOOLSBINARY = System.getenv("BCFTOOLSBINARY");
	private static final String SHAPEIT_BINARY = System.getenv("SHAPEITBINARY");
	private static final String IMPUTE2_BINARY = System.getenv("IMPUTE2BINARY");
	private final static String MINIMAC3_BINARY = System.getenv("MINIMAC3BINARY");
	private final static String MINIMAC4_BINARY = System.getenv("MINIMAC4BINARY");
	private static final String SNPTEST_BINARY = System.getenv("SNPTESTBINARY");
	private static final String JAVA_HOME = System.getenv("JAVA_HOME");
	private static final String SAMTOOLSBINARY = System.getenv("SAMTOOLSBINARY");

	private static final String SEX1 = GuidanceImpl.getSex1();
	private static final String SEX2 = GuidanceImpl.getSex2();
	private static final String NO_SEX = GuidanceImpl.getNoSex();

	private static final boolean ERASE_FILES = false;

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
		GeneralUtils.printGuidancePackageVersion(LOGGER, GUIDANCE_VERSION, GUIDANCE_BUILDNUMBER);
		if (DEBUG) {
			GeneralUtils.printEnVariables(LOGGER);
		}
		// Get the input arguments
		ParseCmdLine parsingArgs = new ParseCmdLine(args);

		// Verify and print the status of each stage
		GeneralUtils.printStagesStatus(LOGGER, parsingArgs);
		LOGGER.info("\n[Guidance] Verifyed stages status.");

		// Get the file name where the list of commands is going to be saved
		// (listOfStagesFile)
		String listOfStagesFileName = parsingArgs.getListOfStagesFile();

		// Verify whether the file exists or not.
		listOfStages = GeneralUtils.createListOfStages(LOGGER, listOfStagesFileName);

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
		GeneralUtils.flushCommands(listOfStages, listOfCommands, true);
		LOGGER.info("[Guidance] Everything is working with Guidance, just wait...");
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

		PhenomeAnalysisFiles phenomeAnalysisFilesInfo = new PhenomeAnalysisFiles(parsingArgs, outDir);

		LOGGER.info("[Guidance] Directory structures and file names have been created...");

		if (BARRIERS) {
			COMPSs.barrier();
		}

		String phasingTool = parsingArgs.getPhasingTool();

		// The number of Chromos to process is endChr - startChr + 1;
		for (int chr = startChr; chr <= endChr; chr++) {
			// We get the output pairs file name for mixed
			String mixedPairsFile = commonFilesInfo.getPairsFile(chr);
			String mixedSampleFile = commonFilesInfo.getSampleFile(chr);

			String theChromo = Integer.toString(chr);
			String bedFile = null;
			String bimFile = null;
			String famFile = null;
			String mixedBedFile = null;
			String mixedBimFile = null;
			String mixedFamFile = null;
			String mixedBedToBedLogFile = null;

			String gmapFile = "";
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
			String mixedPhasingNewSampleFile = commonFilesInfo.getPhasingNewSampleFile(chr);
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

				String responseVar = parsingArgs.getAllResponseVars();
				String covariables = parsingArgs.getAllCovariables();

				if (parsingArgs.getStageStatus("phasingBed") == 1) {

					doPhasingBed(parsingArgs, theChromo, mixedBedFile, mixedBimFile, mixedFamFile, gmapFile,
							mixedSampleFile, covariables, responseVar, mixedExcludedSnpsFile, mixedPhasingHapsFile,
							mixedPhasingSampleFile, mixedPhasingLogFile, mixedFilteredHaplotypesLogFile,
							mixedFilteredHaplotypesVcfFile, mixedFilteredHaplotypesVcfFileBgzip,
							mixedFilteredHaplotypesVcfFileBgzipIndexed, exclCgatFlag, NO_SEX);

					String cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " + mixedPhasingSampleFile + " "
							+ mixedPhasingNewSampleFile + " " + covariables + " " + responseVar;

					listOfCommands.add(new String(cmd));
					GuidanceImpl.newSample(mixedSampleFile, mixedPhasingSampleFile, mixedPhasingNewSampleFile,
							responseVar, covariables, cmd);

					GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
				}

				if (theChromo.equals("23")) {

					if (parsingArgs.getStageStatus("phasingBed") == 1) {
						String mixedPhasingHapsMalesFile = commonFilesInfo.getPhasingHapsMalesFile();
						String mixedPhasingHapsFemalesFile = commonFilesInfo.getPhasingHapsFemalesFile();

						String mixedPhasingSampleMalesFile = commonFilesInfo.getPhasingSampleMalesFile();
						String mixedPhasingSampleFemalesFile = commonFilesInfo.getPhasingSampleFemalesFile();

						String mixedPhasingNewSampleMalesFile = commonFilesInfo.getPhasingNewSampleMalesFile();
						String mixedPhasingNewSampleFemalesFile = commonFilesInfo.getPhasingNewSampleFemalesFile();

						String mixedPhasingLogMalesFile = commonFilesInfo.getPhasingLogMalesFile();
						String mixedPhasingLogFemalesFile = commonFilesInfo.getPhasingLogFemalesFile();

						String mixedExcludedSnpsMalesFile = commonFilesInfo.getExcludedSnpsMalesFile();
						String mixedExcludedSnpsFemalesFile = commonFilesInfo.getExcludedSnpsFemalesFile();

						String mixedFilteredHaplotypesLogMalesFile = commonFilesInfo
								.getFilteredHaplotypesLogMalesFile();
						String mixedFilteredHaplotypesLogFemalesFile = commonFilesInfo
								.getFilteredHaplotypesLogFemalesFile();

						String mixedFilteredHaplotypesVcfMalesFile = commonFilesInfo
								.getFilteredHaplotypesVcfMalesFile();
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
								mixedExcludedSnpsFemalesFile, mixedPhasingHapsFemalesFile,
								mixedPhasingSampleFemalesFile, mixedPhasingLogFemalesFile,
								mixedFilteredHaplotypesLogFemalesFile, mixedFilteredHaplotypesVcfFemalesFile,
								mixedFilteredHaplotypesVcfFemalesFileBgzip,
								mixedFilteredHaplotypesVcfFemalesFileBgzipIndexed, exclCgatFlag, SEX2);

						String cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " + mixedPhasingSampleMalesFile
								+ " " + mixedPhasingNewSampleMalesFile + " " + covariables + " " + responseVar;

						listOfCommands.add(new String(cmd));
						GuidanceImpl.newSample(mixedSampleFile, mixedPhasingSampleMalesFile,
								mixedPhasingNewSampleMalesFile, responseVar, covariables, cmd);

						cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " + mixedPhasingSampleFemalesFile
								+ " " + mixedPhasingNewSampleFemalesFile + " " + covariables + " " + responseVar;

						listOfCommands.add(new String(cmd));

						GuidanceImpl.newSample(mixedSampleFile, mixedPhasingSampleFemalesFile,
								mixedPhasingNewSampleFemalesFile, responseVar, covariables, cmd);

						GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
					}

				}

			}

			GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
		}

		GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);

		if (BARRIERS) {
			COMPSs.barrier();
		}

		// The number of Chromos to process is endChr - startChr + 1;
		for (int panel = 0; panel < rpanelTypes.size(); panel++) {
			for (int chr = startChr; chr <= endChr; chr++) {

				String gmapFile = null;

				if (phasingTool.equals("shapeit")) {
					gmapFile = parsingArgs.getGmapDir() + "/" + parsingArgs.getGmapFileName(chr);
				} else if (phasingTool.equals("eagle")) {
					gmapFile = parsingArgs.getGmapDir() + "/" + parsingArgs.getGmapFileNameEagle();
				}

				String mixedGmapFile = commonFilesInfo.getGmapFile(chr);

				int maxSize = ChromoInfo.getMaxSize(chr);
				int minSize = ChromoInfo.getMinSize(chr);
				int lim1 = minSize;
				int lim2 = lim1 + chunkSize - 1;
				for (int j = minSize; j < maxSize; j = j + chunkSize) {
					makeImputationPerChunk(parsingArgs, chr, lim1, lim2, panel, gmapFile, mixedGmapFile,
							imputationFilesInfo, commonFilesInfo);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
					GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
				}
			} // End for panel types

			if (BARRIERS) {
				COMPSs.barrier();
			}

		} // End for chromosomes

		// *******************************************************************
		// * COMPSs API Call to wait for all tasks *
		// * comment out following line to include the synchronization point *
		// *******************************************************************
		// COMPSs.barrier();

		if (parsingArgs.getStageStatus("snptest") == 1) {

			String responseVar = parsingArgs.getAllResponseVars();
			String covariables = parsingArgs.getAllCovariables();

			for (int chr = startChr; chr <= endChr; chr++) {
				if (chr == 23) {

					String mixedPhasingNewSampleMalesFile = commonFilesInfo.getPhasingNewSampleMalesFile();
					String mixedPhasingNewSampleFemalesFile = commonFilesInfo.getPhasingNewSampleFemalesFile();

					String mixedSampleFile = commonFilesInfo.getSampleFile(chr);

					String mixedPhasingSampleMalesFile = commonFilesInfo.getPhasingSampleMalesFile();
					String mixedPhasingSampleFemalesFile = commonFilesInfo.getPhasingSampleFemalesFile();

					String cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " + mixedPhasingSampleMalesFile
							+ " " + mixedPhasingNewSampleMalesFile + " " + covariables + " " + responseVar;

					listOfCommands.add(new String(cmd));
					try {
						GuidanceImpl.newSample(mixedSampleFile, mixedPhasingSampleMalesFile,
								mixedPhasingNewSampleMalesFile, responseVar, covariables, cmd);
					} catch (Exception e) {
						System.err.println("[Guidance] Exception trying the execution of snptest task");
						System.err.println(e.getMessage());
					}

					cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " + mixedPhasingSampleFemalesFile + " "
							+ mixedPhasingNewSampleFemalesFile + " " + covariables + " " + responseVar;

					listOfCommands.add(new String(cmd));
					try {
						GuidanceImpl.newSample(mixedSampleFile, mixedPhasingSampleFemalesFile,
								mixedPhasingNewSampleFemalesFile, responseVar, covariables, cmd);
					} catch (Exception e) {
						System.err.println("[Guidance] Exception trying the execution of snptest task");
						System.err.println(e.getMessage());
					}
				} // else {
				String mixedPhasingNewSampleFile = commonFilesInfo.getPhasingNewSampleFile(chr);

				String mixedSampleFile = commonFilesInfo.getSampleFile(chr);

				String mixedPhasingSampleFile = commonFilesInfo.getPhasingSampleFile(chr);

				String cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " + mixedPhasingSampleFile + " "
						+ mixedPhasingNewSampleFile + " " + covariables + " " + responseVar;

				listOfCommands.add(new String(cmd));
				try {
					GuidanceImpl.newSample(mixedSampleFile, mixedPhasingSampleFile, mixedPhasingNewSampleFile,
							responseVar, covariables, cmd);
				} catch (Exception e) {
					System.err.println("[Guidance] Exception trying the execution of snptest task");
					System.err.println(e.getMessage());
				}
				// }
			}
		}

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

					for (int j = minSize; j < maxSize; j = j + chunkSize) {
						if (DEBUG) {
							String testName = parsingArgs.getTestTypeName(test);
							LOGGER.debug("Making association for TEST " + testName + ", PANEL " + panelName
									+ ", CHROMO " + chr + ", CHUNK " + lim1 + " - " + lim2);
						}
						makeAssociationPerChunk(parsingArgs, test, panel, panelName, chr, lim1, lim2,
								imputationFilesInfo, commonFilesInfo, assocFilesInfo);

						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}

					GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);

					if (chr == 23) {
						// Now we perform the merge of chunks for each chromosome
						if (parsingArgs.getStageStatus("jointFilteredByAllFiles") == 1) {
							makeMergeOfChunksSex(parsingArgs, test, SEX1, panel, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, FILTERED);

							makeMergeOfChunksSex(parsingArgs, test, SEX2, panel, minSize, maxSize, chunkSize,
									assocFilesInfo, mergeFilesInfo, FILTERED);
						}
					}
					// } else {
					// Now we perform the merge of chunks for each chromosome
					if (parsingArgs.getStageStatus("jointFilteredByAllFiles") == 1) {
						makeMergeOfChunks(parsingArgs, test, panel, chr, minSize, maxSize, chunkSize, assocFilesInfo,
								mergeFilesInfo, FILTERED);
					}
					// }
				} // End for Chromo

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
				String crossRangesResults = resultsFilesInfo.getCrossRangesFile(test, panel);
				String filteredByAllXMalesFile;
				String filteredByAllXFemalesFile;
				if (endChr == 23) {
					filteredByAllXMalesFile = mergeFilesInfo.getAdditionalFilteredByAllXMalesFile(test, panel, 0);
					filteredByAllXFemalesFile = mergeFilesInfo.getAdditionalFilteredByAllXFemalesFile(test, panel, 0);
					if (startChr == 23) {
						lastFilteredByAllFile = filteredByAllXMalesFile;
					}
				} else {
					filteredByAllXMalesFile = lastFilteredByAllFile;
					filteredByAllXFemalesFile = lastFilteredByAllFile;
				}

				if (parsingArgs.getStageStatus("jointCondensedFiles") == 1) {
					doGenerateCondensedAndTopHitsFile(parsingArgs, lastFilteredByAllFile, filteredByAllXMalesFile,
							filteredByAllXFemalesFile, lastCondensedFile, topHitsResults, crossRangesResults);
				}

				// Generate QQManhattan Plots
				String qqPlotPdfFile = resultsFilesInfo.getQqPlotPdfFile(test, panel);
				String qqPlotTiffFile = resultsFilesInfo.getQqPlotTiffFile(test, panel);
				String manhattanPlotPdfFile = resultsFilesInfo.getManhattanPdfFile(test, panel);
				String manhattanPlotTiffFile = resultsFilesInfo.getManhattanTiffFile(test, panel);

				if (parsingArgs.getStageStatus("generateQQManhattanPlots") == 1) {
					doGenerateQQManhattanPlots(parsingArgs, lastCondensedFile, qqPlotPdfFile, manhattanPlotPdfFile,
							qqPlotTiffFile, manhattanPlotTiffFile);
				}

				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);

				if (BARRIERS) {
					COMPSs.barrier();
				}

			} // End for refPanels

			// Now we continue with the combining of the results of the different reference
			// panels.
			// It is done if the refPanelCombine flag is true.
			makeCombinePanels(parsingArgs, assocFilesInfo, mergeFilesInfo, combinedPanelsFilesInfo, rpanelTypes, test);

		} // End for test types

		// GeneralUtils.flushCommands(listOfStages, listOfCommands, true);

		if (1 < numberOfTestTypes) {
			makePhenotypeAnalysis(parsingArgs, combinedPanelsFilesInfo, resultsFilesInfo, phenomeAnalysisFilesInfo);
			GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
		} else {
			LOGGER.info("\n[Guidance] No cross-phenotype analysis. Only one phenotype available");
		}

		LOGGER.info("\n[Guidance] All tasks are in execution, please wait...");

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

		String rpanelDir = parsingArgs.getRpanelDir(panelIndex);

		String lim1S = Integer.toString(lim1);
		String lim2S = Integer.toString(lim2);
		String chrS = Integer.toString(chrNumber);
		// String infoThresholdS = Double.toString(infoThreshold);
		String imputationTool = parsingArgs.getImputationTool();
		String phasingTool = parsingArgs.getPhasingTool();

		// if (imputationTool.equals("impute") || (imputationTool.equals("minimac") &&
		// chrNumber == 23)) {
		// boolean specialCase = (imputationTool.equals("minimac") && chrNumber == 23);
		if (imputationTool.equals("impute")) {

			String knownHapFileName = null;
			String knownHapFile = null;

			String legendFileName = null;
			String legendFile = null;
			// if (specialCase) {
			// knownHapFile = parsingArgs.getRpanelHap23FileName(panelIndex);
			// legendFile = parsingArgs.getRpanelLeg23FileName(panelIndex);
			// } else {
			knownHapFileName = parsingArgs.getRpanelHapFileName(panelIndex, chrNumber);
			knownHapFile = rpanelDir + File.separator + knownHapFileName;

			legendFileName = parsingArgs.getRpanelLegFileName(panelIndex, chrNumber);
			legendFile = rpanelDir + File.separator + legendFileName;
			// }

			String mixedSampleFile = "";
			String mixedPhasingHapsFile = "";
			String mixedPhasingNewSampleFile = "";
			String mixedPairsFile = commonFilesInfo.getPairsFile(chrNumber);
			String mixedImputeFile = "";
			String mixedImputeFileInfo = "";
			String mixedImputeFileSummary = "";
			String mixedImputeFileWarnings = "";
			// String mixedImputeLogFile = "";
			String mixedFilteredFile = "";
			String mixedFilteredLogFile = "";
			String mixedFilteredRsIdFile = "";

			// if (!chrS.equals("23")) {
			mixedSampleFile = commonFilesInfo.getSampleFile(chrNumber);
			mixedPhasingHapsFile = commonFilesInfo.getPhasingHapsFile(chrNumber);
			mixedPhasingNewSampleFile = commonFilesInfo.getPhasingNewSampleFile(chrNumber);
			mixedImputeFile = imputationFilesInfo.getImputedFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
			mixedImputeFileInfo = imputationFilesInfo.getImputedInfoFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
			mixedImputeFileSummary = imputationFilesInfo.getImputedSummaryFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			mixedImputeFileWarnings = imputationFilesInfo.getImputedWarningsFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			/*
			 * mixedImputeLogFile = imputationFilesInfo.getImputedLogFile(panelIndex,
			 * chrNumber, lim1, lim2, chunkSize);
			 */
			mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
			mixedFilteredLogFile = imputationFilesInfo.getFilteredLogFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
			mixedFilteredRsIdFile = imputationFilesInfo.getFilteredRsIdFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			// }

			// We create the list of rsId that are greater than or equal to the
			// infoThreshold value
			if (phasingTool.equals("eagle")) {

				// if (!chrS.equals("23")) {
				doImputationWithImpute(parsingArgs, chrS, mixedGmapFile, knownHapFile, legendFile, mixedPhasingHapsFile,
						mixedSampleFile, lim1S, lim2S, mixedPairsFile, mixedImputeFile, mixedImputeFileInfo,
						mixedImputeFileSummary, mixedImputeFileWarnings, NO_SEX, panelIndex);
				doFilterByInfo(parsingArgs, mixedImputeFileInfo, mixedFilteredRsIdFile, chrS);
				doQctoolS(parsingArgs, mixedImputeFile, mixedFilteredRsIdFile, mixedFilteredFile, mixedFilteredLogFile,
						chrS);

				// } else if (chrS.equals("23")) {
				if (chrS.equals("23")) {

					String mixedPhasingHapsMalesFile = commonFilesInfo.getPhasingHapsMalesFile();
					String mixedPhasingHapsFemalesFile = commonFilesInfo.getPhasingHapsFemalesFile();

					String mixedPhasingNewSampleMalesFile = commonFilesInfo.getPhasingNewSampleMalesFile();
					String mixedPhasingNewSampleFemalesFile = commonFilesInfo.getPhasingNewSampleFemalesFile();

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
							mixedPhasingHapsMalesFile, mixedPhasingNewSampleMalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeMalesFile, mixedImputeMalesFileInfo, mixedImputeMalesFileSummary,
							mixedImputeMalesFileWarnings, SEX1, panelIndex);

					doFilterByInfo(parsingArgs, mixedImputeMalesFileInfo, mixedFilteredRsIdMalesFile, chrS);

					doQctoolS(parsingArgs, mixedImputeMalesFile, mixedFilteredRsIdMalesFile, mixedFilteredMalesFile,
							mixedFilteredLogMalesFile, chrS);

					doImputationWithImpute(parsingArgs, chrS, mixedGmapFile, knownHapFile, legendFile,
							mixedPhasingHapsFemalesFile, mixedPhasingNewSampleFemalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeFemalesFile, mixedImputeFemalesFileInfo, mixedImputeFemalesFileSummary,
							mixedImputeFemalesFileWarnings, SEX2, panelIndex);

					doFilterByInfo(parsingArgs, mixedImputeFemalesFileInfo, mixedFilteredRsIdFemalesFile, chrS);

					doQctoolS(parsingArgs, mixedImputeFemalesFile, mixedFilteredRsIdFemalesFile,
							mixedFilteredFemalesFile, mixedFilteredLogFemalesFile, chrS);
				}

			} else if (phasingTool.equals("shapeit")) {
				// if (!chrS.equals("23")) {
				doImputationWithImpute(parsingArgs, chrS, gmapFile, knownHapFile, legendFile, mixedPhasingHapsFile,
						mixedPhasingNewSampleFile, lim1S, lim2S, mixedPairsFile, mixedImputeFile, mixedImputeFileInfo,
						mixedImputeFileSummary, mixedImputeFileWarnings, NO_SEX, panelIndex);
				doFilterByInfo(parsingArgs, mixedImputeFileInfo, mixedFilteredRsIdFile, chrS);
				doQctoolS(parsingArgs, mixedImputeFile, mixedFilteredRsIdFile, mixedFilteredFile, mixedFilteredLogFile,
						chrS);
				// } else if (chrS.equals("23")) {
				if (chrS.equals("23")) {
					String mixedPhasingHapsMalesFile = commonFilesInfo.getPhasingHapsMalesFile();
					String mixedPhasingHapsFemalesFile = commonFilesInfo.getPhasingHapsFemalesFile();

					String mixedPhasingNewSampleMalesFile = commonFilesInfo.getPhasingNewSampleMalesFile();
					String mixedPhasingNewSampleFemalesFile = commonFilesInfo.getPhasingNewSampleFemalesFile();

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
							mixedPhasingHapsMalesFile, mixedPhasingNewSampleMalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeMalesFile, mixedImputeMalesFileInfo, mixedImputeMalesFileSummary,
							mixedImputeMalesFileWarnings, SEX1, panelIndex);

					doFilterByInfo(parsingArgs, mixedImputeMalesFileInfo, mixedFilteredRsIdMalesFile, chrS);

					doQctoolS(parsingArgs, mixedImputeMalesFile, mixedFilteredRsIdMalesFile, mixedFilteredMalesFile,
							mixedFilteredLogMalesFile, chrS);

					doImputationWithImpute(parsingArgs, chrS, gmapFile, knownHapFile, legendFile,
							mixedPhasingHapsFemalesFile, mixedPhasingNewSampleFemalesFile, lim1S, lim2S, mixedPairsFile,
							mixedImputeFemalesFile, mixedImputeFemalesFileInfo, mixedImputeFemalesFileSummary,
							mixedImputeFemalesFileWarnings, SEX2, panelIndex);

					doFilterByInfo(parsingArgs, mixedImputeFemalesFileInfo, mixedFilteredRsIdFemalesFile, chrS);

					doQctoolS(parsingArgs, mixedImputeFemalesFile, mixedFilteredRsIdFemalesFile,
							mixedFilteredFemalesFile, mixedFilteredLogFemalesFile, chrS);
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
			// if (!chrS.equals("23")) {

			String mixedFilteredHaplotypesVcfFileBgzip = commonFilesInfo.getFilteredHaplotypesVcfFileBgzip(chrNumber);
			// String mixedImputeMMDoseVCFFile =
			// imputationFilesInfo.getImputedMMDoseVCFFile(panelIndex, chrNumber,
			// lim1, lim2, chunkSize);
			String mixedImputeMMInfoFile = imputationFilesInfo.getImputedMMInfoFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeMMErateFile = imputationFilesInfo.getImputedMMErateFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeMMRecFile = imputationFilesInfo.getImputedMMRecFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeMMM3VCFFile = imputationFilesInfo.getImputedMMM3VCFFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeMMLogFile = imputationFilesInfo.getImputedMMLogFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedFilteredLogFile = imputationFilesInfo.getFilteredLogFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedFilteredRsIdFile = imputationFilesInfo.getFilteredRsIdFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeFileBgzip = imputationFilesInfo.getImputedFileBgzip(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			// String mixedImputeMMM3VCFFileBgzip =
			// imputationFilesInfo.getImputedMMM3VCFFileBgzip(panelIndex,
			// chrNumber, lim1, lim2, chunkSize);
			String mixedImputeFileTbi = imputationFilesInfo.getImputedFileTbi(panelIndex, chrNumber, lim1, lim2,
					chunkSize);

			doImputationWithMinimac(parsingArgs, refVcfFile, mixedFilteredHaplotypesVcfFileBgzip, chrS, lim1S, lim2S,
					mixedImputeMMInfoFile, mixedImputeMMErateFile, mixedImputeMMRecFile, mixedImputeMMM3VCFFile,
					mixedImputeMMLogFile, mixedImputeFileBgzip, mixedImputeFileTbi, NO_SEX, panelIndex);

			doFilterByInfo(parsingArgs, mixedImputeMMInfoFile, mixedFilteredRsIdFile, chrS);

			doQctoolS(parsingArgs, mixedImputeFileBgzip, mixedFilteredRsIdFile, mixedFilteredFile, mixedFilteredLogFile,
					chrS);

			// }

			if (chrS.equals("23")) {
				String mixedFilteredHaplotypesVcfMalesFileBgzip = commonFilesInfo
						.getFilteredHaplotypesVcfMalesFileBgzip();
				String mixedFilteredHaplotypesVcfFemalesFileBgzip = commonFilesInfo
						.getFilteredHaplotypesVcfFemalesFileBgzip();

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

				String mixedImputeMalesFileTbi = imputationFilesInfo.getImputedMalesFileTbi(panelIndex, chrNumber, lim1,
						lim2, chunkSize);
				String mixedImputeFemalesFileTbi = imputationFilesInfo.getImputedFemalesFileTbi(panelIndex, chrNumber,
						lim1, lim2, chunkSize);

				doImputationWithMinimac(parsingArgs, refVcfFile, mixedFilteredHaplotypesVcfMalesFileBgzip, chrS,
						lim1S, lim2S, mixedImputeMMInfoMalesFile, mixedImputeMMErateMalesFile,
						mixedImputeMMRecMalesFile, mixedImputeMMM3VCFMalesFile, mixedImputeMMLogMalesFile,
						mixedImputeMalesFileBgzip, mixedImputeMalesFileTbi, SEX1, panelIndex);

				doFilterByInfo(parsingArgs, mixedImputeMMInfoMalesFile, mixedFilteredRsIdMalesFile, chrS);

				doQctoolS(parsingArgs, mixedImputeMalesFileBgzip, mixedFilteredRsIdMalesFile, mixedFilteredMalesFile,
						mixedFilteredLogMalesFile, chrS);

				doImputationWithMinimac(parsingArgs, refVcfFile, mixedFilteredHaplotypesVcfFemalesFileBgzip, chrS,
						lim1S, lim2S, mixedImputeMMInfoFemalesFile, mixedImputeMMErateFemalesFile,
						mixedImputeMMRecFemalesFile, mixedImputeMMM3VCFFemalesFile, mixedImputeMMLogFemalesFile,
						mixedImputeFemalesFileBgzip, mixedImputeFemalesFileTbi, SEX2, panelIndex);

				doFilterByInfo(parsingArgs, mixedImputeMMInfoFemalesFile, mixedFilteredRsIdFemalesFile, chrS);

				doQctoolS(parsingArgs, mixedImputeFemalesFileBgzip, mixedFilteredRsIdFemalesFile,
						mixedFilteredFemalesFile, mixedFilteredLogFemalesFile, chrS);

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
		// double infoThreshold = parsingArgs.getInfoThreshold();
		double hweCohortThreshold = parsingArgs.getHweCohortThreshold();
		double hweCasesThreshold = parsingArgs.getHweCasesThreshold();
		double hweControlsThreshold = parsingArgs.getHweControlsThreshold();

		String responseVar = parsingArgs.getResponseVar(testTypeIndex);
		String covariables = parsingArgs.getCovariables(testTypeIndex);

		String chrS = Integer.toString(chrNumber);

		String mafThresholdS = Double.toString(mafThreshold);
		// String infoThresholdS = Double.toString(infoThreshold);
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

			String mixedPhasingNewSampleMalesFile = commonFilesInfo.getPhasingNewSampleMalesFile();
			String mixedPhasingNewSampleFemalesFile = commonFilesInfo.getPhasingNewSampleFemalesFile();

			doSnptest(parsingArgs, chrS, mixedFilteredMalesFile, mixedPhasingNewSampleMalesFile, snptestOutMalesFile,
					snptestLogMalesFile, responseVar, covariables);

			doSnptest(parsingArgs, chrS, mixedFilteredFemalesFile, mixedPhasingNewSampleFemalesFile,
					snptestOutFemalesFile, snptestLogFemalesFile, responseVar, covariables);

			// This true should be erased when the problem with minimac is solved
			// if (imputationTool.equals("impute") || true) {
			if (imputationTool.equals("impute")) {
				String mixedImputeMalesFileInfo = imputationFilesInfo.getImputedInfoMalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);

				String summaryMalesFile = assocFilesInfo.getSummaryMalesFile(testTypeIndex, panelIndex, lim1, lim2,
						chunkSize);

				doCollectSummary(parsingArgs, chrS, mixedImputeMalesFileInfo, snptestOutMalesFile, summaryMalesFile,
						mafThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, SEX1);

				String mixedImputeFemalesFileInfo = imputationFilesInfo.getImputedInfoFemalesFile(panelIndex, chrNumber,
						lim1, lim2, chunkSize);

				String summaryFemalesFile = assocFilesInfo.getSummaryFemalesFile(testTypeIndex, panelIndex, lim1, lim2,
						chunkSize);

				doCollectSummary(parsingArgs, chrS, mixedImputeFemalesFileInfo, snptestOutFemalesFile,
						summaryFemalesFile, mafThresholdS, hweCohortThresholdS, hweCasesThresholdS,
						hweControlsThresholdS, SEX2);

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
						mafThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, SEX1);

				doCollectSummary(parsingArgs, chrS, mixedImputedMMInfoFemalesFile, snptestOutFemalesFile,
						summaryFemalesFile, mafThresholdS, hweCohortThresholdS, hweCasesThresholdS,
						hweControlsThresholdS, SEX2);

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

		} // else {
		String snptestOutFile = assocFilesInfo.getSnptestOutFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
				chunkSize);
		String snptestLogFile = assocFilesInfo.getSnptestLogFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
				chunkSize);
		String mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
		String mixedPhasingNewSampleFile = commonFilesInfo.getPhasingNewSampleFile(chrNumber);

		// String mixedSampleFile = commonFilesInfo.getSampleFile(chrNumber);

		// String mixedPhasingSampleFile =
		// commonFilesInfo.getPhasingSampleFile(chrNumber);

		// String cmd = JAVA_HOME + " newSample.jar " + mixedSampleFile + " " +
		// mixedPhasingSampleFile + " "
		// + mixedPhasingNewSampleFile + " " + covariables + " " + responseVar;

		doSnptest(parsingArgs, chrS, mixedFilteredFile, mixedPhasingNewSampleFile, snptestOutFile, snptestLogFile,
				responseVar, covariables);

		if (imputationTool.equals("impute")) {
			String mixedImputeFileInfo = imputationFilesInfo.getImputedInfoFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String summaryFile = assocFilesInfo.getSummaryFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
					chunkSize);

			doCollectSummary(parsingArgs, chrS, mixedImputeFileInfo, snptestOutFile, summaryFile, mafThresholdS,
					hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, NO_SEX);

			String assocFilteredByAll = assocFilesInfo.getSummaryFilteredFile(testTypeIndex, panelIndex, chrNumber,
					lim1, lim2, chunkSize);
			String assocCondensed = assocFilesInfo.getSummaryCondensedFile(testTypeIndex, panelIndex, chrNumber, lim1,
					lim2, chunkSize);

			doFilterByAll(parsingArgs, summaryFile, assocFilteredByAll, assocCondensed, NO_SEX, rpanelName);

		} else if (imputationTool.equals("minimac")) {
			String mixedImputedMMInfoFile = imputationFilesInfo.getImputedMMInfoFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String summaryFile = assocFilesInfo.getSummaryFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
					chunkSize);

			doCollectSummary(parsingArgs, chrS, mixedImputedMMInfoFile, snptestOutFile, summaryFile, mafThresholdS,
					hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, NO_SEX);

			String assocFilteredByAll = assocFilesInfo.getSummaryFilteredFile(testTypeIndex, panelIndex, chrNumber,
					lim1, lim2, chunkSize);
			String assocCondensed = assocFilesInfo.getSummaryCondensedFile(testTypeIndex, panelIndex, chrNumber, lim1,
					lim2, chunkSize);

			doFilterByAll(parsingArgs, summaryFile, assocFilteredByAll, assocCondensed, NO_SEX, rpanelName);
		}

		// }

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
	 * @throws GuidanceTaskException
	 */
	private static void makeMergeOfChunksSex(ParseCmdLine parsingArgs, int ttIndex, String sex, int rpanelIndex,
			int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo, String type)
			throws GuidanceTaskException {
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
					boolean eraseA = false;
					boolean eraseB = false;
					if (processedChunks < numberOfChunks) {
						reducedA = assocFilesInfo.getSummaryFilteredMalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedA = mergeFilesInfo.getReducedFilteredMalesFile(ttIndex, rpanelIndex, indexA);
						eraseA = true;
						indexA++;
					}
					if (processedChunks < numberOfChunks - 1) {
						reducedB = assocFilesInfo.getSummaryFilteredMalesFile(ttIndex, rpanelIndex, lim1, lim2,
								chunkSize);
						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					} else {
						reducedB = mergeFilesInfo.getReducedFilteredMalesFile(ttIndex, rpanelIndex, indexA);
						eraseB = true;
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
					if (eraseA && ERASE_FILES) {
						File f = new File(reducedA);
						f.delete();

					}
					if (eraseB && ERASE_FILES) {
						File f = new File(reducedB);
						f.delete();

					}

					// Clean intermediate files
					// File fA = new File(reducedA);
					// fA.delete();
					// File fB = new File(reducedB);
					// fB.delete();

				} // End for Chunks

			}
			/*
			 * else if (type.equals(CONDENSED)) { String condensedFile =
			 * mergeFilesInfo.getCondensedMalesFile(ttIndex, rpanelIndex); for (int
			 * processedChunks = 0; processedChunks < 2 * numberOfChunks - 2;
			 * processedChunks = processedChunks + 2) { if (processedChunks <
			 * numberOfChunks) { reducedA =
			 * assocFilesInfo.getSummaryCondensedMalesFile(ttIndex, rpanelIndex, lim1, lim2,
			 * chunkSize); lim1 = lim1 + chunkSize; lim2 = lim2 + chunkSize; } else {
			 * reducedA = mergeFilesInfo.getReducedCondensedMalesFile(ttIndex, rpanelIndex,
			 * indexA); indexA++; } if (processedChunks < numberOfChunks - 1) { reducedB =
			 * assocFilesInfo.getSummaryCondensedMalesFile(ttIndex, rpanelIndex, lim1, lim2,
			 * chunkSize); lim1 = lim1 + chunkSize; lim2 = lim2 + chunkSize; } else {
			 * reducedB = mergeFilesInfo.getReducedCondensedMalesFile(ttIndex, rpanelIndex,
			 * indexA); indexA++; }
			 * 
			 * if (processedChunks == 2 * numberOfChunks - 4) {
			 * doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB,
			 * condensedFile); indexC++; } else { reducedC =
			 * mergeFilesInfo.getReducedCondensedMalesFile(ttIndex, rpanelIndex, indexC);
			 * doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
			 * indexC++; }
			 * 
			 * // Clean intermediate files // File fA = new File(reducedA); // fA.delete();
			 * // File fB = new File(reducedB); // fB.delete();
			 * 
			 * } // End of for Chunks }
			 */
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

			}
			/*
			 * else if (type.equals(CONDENSED)) { String condensedFile =
			 * mergeFilesInfo.getCondensedFemalesFile(ttIndex, rpanelIndex); for (int
			 * processedChunks = 0; processedChunks < 2 * numberOfChunks - 2;
			 * processedChunks = processedChunks + 2) { if (processedChunks <
			 * numberOfChunks) { reducedA =
			 * assocFilesInfo.getSummaryCondensedFemalesFile(ttIndex, rpanelIndex, lim1,
			 * lim2, chunkSize); lim1 = lim1 + chunkSize; lim2 = lim2 + chunkSize; } else {
			 * reducedA = mergeFilesInfo.getReducedCondensedFemalesFile(ttIndex,
			 * rpanelIndex, indexA); indexA++; } if (processedChunks < numberOfChunks - 1) {
			 * reducedB = assocFilesInfo.getSummaryCondensedFemalesFile(ttIndex,
			 * rpanelIndex, lim1, lim2, chunkSize); lim1 = lim1 + chunkSize; lim2 = lim2 +
			 * chunkSize; } else { reducedB =
			 * mergeFilesInfo.getReducedCondensedFemalesFile(ttIndex, rpanelIndex, indexA);
			 * indexA++; }
			 * 
			 * if (processedChunks == 2 * numberOfChunks - 4) {
			 * doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB,
			 * condensedFile); indexC++; } else { reducedC =
			 * mergeFilesInfo.getReducedCondensedFemalesFile(ttIndex, rpanelIndex, indexC);
			 * doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB, reducedC);
			 * indexC++; }
			 * 
			 * // Clean intermediate files // File fA = new File(reducedA); // fA.delete();
			 * // File fB = new File(reducedB); // fB.delete();
			 * 
			 * } // End of for Chunks }
			 */
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
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
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
	 * @throws GuidanceTaskException
	 */
	private static void makeMergeOfChunksCombinedSex(ParseCmdLine parsingArgs, int ttIndex, int rpanelIndex,
			int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo, MergeFiles mergeFilesInfo, String type)
			throws GuidanceTaskException {

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
		// String condensedFileMales =
		// mergeFilesInfo.getCombinedCondensedMalesFile(ttIndex, rpanelIndex);

		if (type.equals(FILTERED)) {
			// LOGGER.info("Number of chunks for testType " + ttIndex + " | rpanel " +
			// rpanelIndex + " |chr " +
			// chr + " " + numberOfChunks);
			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
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
				if (ERASE_FILES) {
					File f = new File(reducedA);
					f.delete();
					f = new File(reducedB);
					f.delete();
				}
			} // End for Chunks

			indexA = 0;
			indexC = 0;

			lim1 = minSize;
			lim2 = lim1 + chunkSize - 1;

			for (int processedChunks = 0; processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks
					+ 2) {
				if (processedChunks < numberOfChunks) {
					reducedA = assocFilesInfo.getCombinedFilteredFemalesFile(ttIndex, rpanelIndex, lim1, lim2,
							chunkSize);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} else {
					reducedA = mergeFilesInfo.getCombinedReducedFilteredFemalesFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}
				if (processedChunks < numberOfChunks - 1) {
					reducedB = assocFilesInfo.getCombinedFilteredFemalesFile(ttIndex, rpanelIndex, lim1, lim2,
							chunkSize);
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
				if (ERASE_FILES) {
					File f = new File(reducedA);
					f.delete();
					f = new File(reducedB);
					f.delete();
				}

				// Clean intermediate files
				// File fA = new File(reducedA);
				// fA.delete();
				// File fB = new File(reducedB);
				// fB.delete();

			} // End for Chunks

		}
		/*
		 * else if (type.equals(CONDENSED)) { for (int processedChunks = 0;
		 * processedChunks < 2 * numberOfChunks - 2; processedChunks = processedChunks +
		 * 2) { if (processedChunks < numberOfChunks) { reducedA =
		 * assocFilesInfo.getCombinedCondensedMalesFile(ttIndex, rpanelIndex, lim1,
		 * lim2, chunkSize); lim1 = lim1 + chunkSize; lim2 = lim2 + chunkSize; } else {
		 * reducedA = mergeFilesInfo.getCombinedReducedCondensedMalesFile(ttIndex,
		 * rpanelIndex, indexA); indexA++; } if (processedChunks < numberOfChunks - 1) {
		 * reducedB = assocFilesInfo.getCombinedCondensedMalesFile(ttIndex, rpanelIndex,
		 * lim1, lim2, chunkSize); lim1 = lim1 + chunkSize; lim2 = lim2 + chunkSize; }
		 * else { reducedB =
		 * mergeFilesInfo.getCombinedReducedCondensedMalesFile(ttIndex, rpanelIndex,
		 * indexA); indexA++; }
		 * 
		 * if (processedChunks == 2 * numberOfChunks - 4) {
		 * doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB,
		 * condensedFileMales); indexC++; } else { reducedC =
		 * mergeFilesInfo.getCombinedReducedCondensedMalesFile(ttIndex, rpanelIndex,
		 * indexC); doMergeTwoChunksUnconditional(parsingArgs, reducedA, reducedB,
		 * reducedC); indexC++; }
		 * 
		 * // Clean intermediate files // File fA = new File(reducedA); // fA.delete();
		 * // File fB = new File(reducedB); // fB.delete();
		 * 
		 * } // End of for Chunks }
		 */
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

		if (numberOfChrs == 1 && startChr != 23) { // There is only one chr to process.
			rpanelFlag = "YES";
			filteredA = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, startChr);
			filteredB = filteredA;
			filteredC = mergeFilesInfo.getAdditionalFilteredByAllFile(ttIndex, rpanelIndex, indexC);
			doJointFilteredByAllFiles(parsingArgs, filteredA, filteredB, filteredC, rpanelName, rpanelFlag);
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

			filteredA = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, 23);
			filteredC = mergeFilesInfo.getAdditionalFilteredByAllXFile(ttIndex, rpanelIndex);

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

		// PLACES TO STORE ALL COMBINED INFORMATION
		String filteredCombineAll = combinedPanelsFilesInfo.getCombinedFilteredByAllFile(ttIndex);
		String condensedCombineAll = combinedPanelsFilesInfo.getCombinedCondensedFile(ttIndex);

		// We use same file when we do not need to compute the chromo 23 (doTopHits task
		// uses this information)
		String filteredCombineAllXMales = filteredCombineAll;
		String filteredCombineAllXFemales = filteredCombineAll;
		String filteredCombineAllX = filteredCombineAll;

		// We combine the panels per chromosome
		final int startChr = parsingArgs.getStart();
		final int endChr = parsingArgs.getEnd();
		final int chunkSize = parsingArgs.getChunkSize();

		if (endChr == 23) {
			filteredCombineAllXMales = combinedPanelsFilesInfo.getCombinedFilteredByAllXMalesFile(ttIndex);
			filteredCombineAllXFemales = combinedPanelsFilesInfo.getCombinedFilteredByAllXFemalesFile(ttIndex);
			filteredCombineAllX = combinedPanelsFilesInfo.getCombinedFilteredByAllXFile(ttIndex);
			if (startChr == 23) {
				filteredCombineAll = filteredCombineAllXMales;
			}
		}

		final boolean refPanelCombine = parsingArgs.getRefPanelCombine();
		if (refPanelCombine && parsingArgs.getStageStatus("combinePanelsComplex") == 1) {
			// We are not asked to combine the panels. End

			// INITIALIZE THE FILES TO STORE ALL COMBINED INFORMATION (ADD HEADER AND
			// COMPRESS)

			LinkedList<String> filteredCombined = new LinkedList<>();
			// LinkedList<String> condensedCombined = new LinkedList<>();

			// CHR LOOP
			for (int chr = startChr; chr <= endChr; chr++) {
				int minSize = ChromoInfo.getMinSize(chr);
				int maxSize = ChromoInfo.getMaxSize(chr);

				int lim1 = minSize;
				int lim2 = lim1 + chunkSize - 1;

				// CHUNK LOOP
				for (int j = minSize; j < maxSize; j = j + chunkSize) {
					// -- FILTERED PART --
					// Construct three queues with all the filtered panels to combine in order to
					// get the combined filtered files
					LinkedList<String> filteredPanelsToCombine = new LinkedList<>();
					LinkedList<String> filteredPanelsToCombineMales = new LinkedList<>();
					LinkedList<String> filteredPanelsToCombineFemales = new LinkedList<>();
					LinkedList<String> filteredPanelsToCombineX = new LinkedList<>();

					for (int k = 0; k < rpanelTypes.size(); ++k) {
						if (chr == 23) {
							String filteredPanelMales = assocFilesInfo.getSummaryFilteredMalesFile(ttIndex, k, lim1,
									lim2, chunkSize);
							if (DEBUG) {
								System.out.println("[Guidance] Adding filtered file to combine " + filteredPanelMales);
							}
							filteredPanelsToCombineMales.add(filteredPanelMales);

							String filteredPanelFemales = assocFilesInfo.getSummaryFilteredFemalesFile(ttIndex, k, lim1,
									lim2, chunkSize);
							if (DEBUG) {
								System.out
										.println("[Guidance] Adding filtered file to combine " + filteredPanelFemales);
							}
							filteredPanelsToCombineFemales.add(filteredPanelFemales);

							String filteredPanelX = assocFilesInfo.getSummaryFilteredFile(ttIndex, k, chr, lim1, lim2,
									chunkSize);
							if (DEBUG) {
								System.out
										.println("[Guidance] Adding filtered file to combine " + filteredPanelFemales);
							}
							filteredPanelsToCombineX.add(filteredPanelX);
						} // else {
						String filteredPanel = assocFilesInfo.getSummaryFilteredFile(ttIndex, k, chr, lim1, lim2,
								chunkSize);
						if (DEBUG) {
							LOGGER.debug("[Guidance] Adding filtered file to combine " + filteredPanel);
						}
						filteredPanelsToCombine.add(filteredPanel);
						// }
					}

					GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);

					// Combine all the filtered panels 2 by 2 until there are no remaining panels
					if (chr == 23) {
						String destFilteredPanelMales = assocFilesInfo.getCombinedFilteredMalesFile(ttIndex, 0, lim1,
								lim2, chunkSize);
						String baseFilteredPanelMales = destFilteredPanelMales.substring(0,
								destFilteredPanelMales.length() - 7);
						int counter = 0;
						while (filteredPanelsToCombineMales.size() > 1) {
							++counter;
							String filteredPanelC;
							if (filteredPanelsToCombineMales.size() == 2) {
								filteredPanelC = destFilteredPanelMales;
							} else {
								filteredPanelC = baseFilteredPanelMales + "_reduce_" + Integer.toString(counter)
										+ ".txt.tgz";
							}

							String filteredPanelA = filteredPanelsToCombineMales.poll();
							String filteredPanelB = filteredPanelsToCombineMales.poll();
							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB
										+ " to " + filteredPanelC);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelC, lim1,
									lim2);
							if (ERASE_FILES) {
								File f = new File(filteredPanelA);
								f.delete();
								f = new File(filteredPanelB);
								f.delete();
							}
							filteredPanelsToCombineMales.add(filteredPanelC);
						}

						String destFilteredPanelFemales = assocFilesInfo.getCombinedFilteredFemalesFile(ttIndex, 0,
								lim1, lim2, chunkSize);
						String baseFilteredPanelFemales = destFilteredPanelFemales.substring(0,
								destFilteredPanelFemales.length() - 7);
						counter = 0;
						while (filteredPanelsToCombineFemales.size() > 1) {
							++counter;
							String filteredPanelC;

							if (filteredPanelsToCombineFemales.size() == 2) {
								filteredPanelC = destFilteredPanelFemales;
							} else {
								filteredPanelC = baseFilteredPanelFemales + "_reduce_" + Integer.toString(counter)
										+ ".txt.tgz";
							}
							String filteredPanelA = filteredPanelsToCombineFemales.poll();
							String filteredPanelB = filteredPanelsToCombineFemales.poll();
							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB
										+ " to " + filteredPanelC);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelC, lim1,
									lim2);
							if (ERASE_FILES) {
								File f = new File(filteredPanelA);
								f.delete();
								f = new File(filteredPanelB);
								f.delete();
							}
							filteredPanelsToCombineFemales.add(filteredPanelC);
						}

						String destFilteredPanelX = assocFilesInfo.getCombinedFilteredFile(ttIndex, 0, chr, lim1, lim2,
								chunkSize);
						String baseFilteredPanelX = destFilteredPanelX.substring(0, destFilteredPanelX.length() - 7);
						counter = 0;
						while (filteredPanelsToCombineX.size() > 1) {
							++counter;
							String filteredPanelC;

							if (filteredPanelsToCombineX.size() == 2) {
								filteredPanelC = destFilteredPanelX;
							} else {
								filteredPanelC = baseFilteredPanelX + "_reduce_" + Integer.toString(counter)
										+ ".txt.tgz";
							}
							String filteredPanelA = filteredPanelsToCombineX.poll();
							String filteredPanelB = filteredPanelsToCombineX.poll();
							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB
										+ " to " + filteredPanelC);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelC, lim1,
									lim2);
							if (ERASE_FILES) {
								File f = new File(filteredPanelA);
								f.delete();
								f = new File(filteredPanelB);
								f.delete();
							}
							filteredPanelsToCombineX.add(filteredPanelC);
						}
					} else {
						String destFilteredPanel = assocFilesInfo.getCombinedFilteredFile(ttIndex, 0, chr, lim1, lim2,
								chunkSize);
						String baseFilteredPanel = destFilteredPanel.substring(0, destFilteredPanel.length() - 7);
						int counter = 0;
						while (filteredPanelsToCombine.size() > 1) {
							++counter;
							String filteredPanelC;
							if (filteredPanelsToCombine.size() == 2) {
								filteredPanelC = destFilteredPanel;
							} else {
								filteredPanelC = baseFilteredPanel + "_reduce_" + Integer.toString(counter)
										+ ".txt.tgz";
							}

							String filteredPanelA = filteredPanelsToCombine.poll();
							String filteredPanelB = filteredPanelsToCombine.poll();

							// Filtered part: combines A and B into A
							if (DEBUG) {
								LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB
										+ " to " + filteredPanelC);
							}
							doCombinePanelsComplex(parsingArgs, filteredPanelA, filteredPanelB, filteredPanelC, lim1,
									lim2);
							if (ERASE_FILES) {
								File f = new File(filteredPanelA);
								f.delete();
								f = new File(filteredPanelB);
								f.delete();
							}
							filteredPanelsToCombine.add(filteredPanelC);
						}
					}

					GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);

					// Clean partial results
					// new File(chunkResultsCondensed).delete();

					// Increase loop variables
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				} // End for chunk

				if (chr == 23) {
					makeMergeOfChunksCombinedSex(parsingArgs, ttIndex, 0, minSize, maxSize, chunkSize, assocFilesInfo,
							mergeFilesInfo, FILTERED);
				} // else {
				makeMergeOfChunksCombined(parsingArgs, ttIndex, 0, chr, minSize, maxSize, chunkSize, assocFilesInfo,
						mergeFilesInfo, FILTERED);
				// }
				if (chr == 23) {
					String filteredByAllCurrentMalesFile = mergeFilesInfo.getCombinedFilteredByAllMalesFile(ttIndex, 0);
					doCopyFile(parsingArgs, filteredByAllCurrentMalesFile, filteredCombineAllXMales);

					String filteredByAllCurrentFemalesFile = mergeFilesInfo.getCombinedFilteredByAllFemalesFile(ttIndex,
							0);
					doCopyFile(parsingArgs, filteredByAllCurrentFemalesFile, filteredCombineAllXFemales);

					String filteredByAllCurrentFile = mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, 0, chr);
					doCopyFile(parsingArgs, filteredByAllCurrentFile, filteredCombineAllX);
				} else {
					String filteredByAllCurrentFile = mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, 0, chr);
					filteredCombined.add(filteredByAllCurrentFile);
				}

			} // End for chromosomes

			if (filteredCombined.size() == 1)

			{
				String singleFilteredFile = filteredCombined.peek();
				String destinationFilteredFile = filteredCombineAll;
				doCopyFile(parsingArgs, singleFilteredFile, destinationFilteredFile);
			}

			int reduceCounter = 0;
			while (filteredCombined.size() > 1) {
				String originFilteredFileA = filteredCombined.poll();
				String originFilteredFileB = filteredCombined.poll();

				String destinationFilteredFile;
				boolean erase = false;
				if (filteredCombined.isEmpty()) {
					destinationFilteredFile = filteredCombineAll;
				} else {
					destinationFilteredFile = filteredCombineAll.substring(0, filteredCombineAll.length() - 7)
							+ "_reduce_" + Integer.toString(reduceCounter) + ".txt.gz";
					erase = true;
				}

				doMergeTwoChunksUnconditional(parsingArgs, originFilteredFileA, originFilteredFileB,
						destinationFilteredFile);
				filteredCombined.add(destinationFilteredFile);
				if (erase && ERASE_FILES) {
					File f = new File(destinationFilteredFile);
					f.delete();
				}
				reduceCounter += 1;
			}
		}

		if (refPanelCombine && parsingArgs.getStageStatus("combGenerateManhattanTop") == 1) {

			String topHitsCombinedResults = combinedPanelsFilesInfo.getTopHitsFile(ttIndex);
			String topHitsCrossRanges = combinedPanelsFilesInfo.getCrossRangesFile(ttIndex);

			doGenerateCondensedAndTopHitsFile(parsingArgs, filteredCombineAll, filteredCombineAllXMales,
					filteredCombineAllXFemales, condensedCombineAll, topHitsCombinedResults, topHitsCrossRanges);

			String combinedQqPlotPdfFile = combinedPanelsFilesInfo.getQqPlotPdfFile(ttIndex);
			String combinedQqPlotTiffFile = combinedPanelsFilesInfo.getQqPlotTiffFile(ttIndex);
			String combinedManhattanPlotPdfFile = combinedPanelsFilesInfo.getManhattanPdfFile(ttIndex);
			String combinedManhattanPlotTiffFile = combinedPanelsFilesInfo.getManhattanTiffFile(ttIndex);

			doGenerateQQManhattanPlots(parsingArgs, condensedCombineAll, combinedQqPlotPdfFile,
					combinedManhattanPlotPdfFile, combinedQqPlotTiffFile, combinedManhattanPlotTiffFile);

		}
	}

	/**
	 * Generate the condensed file from the filtered files
	 * 
	 * @param parsingArgs
	 * @param filteredFile
	 * @param filteredMalesFile
	 * @param filteredFemalesFile
	 * @param condensedFile
	 */
	private static void doGenerateCondensedAndTopHitsFile(ParseCmdLine parsingArgs, String filteredFile,
			String filteredMalesFile, String filteredFemalesFile, String condensedFile, String topHitsFile,
			String crossRanges) {

		String pvaThreshold = Double.toString(parsingArgs.getPvaThreshold());

		String cmdToStore = R_SCRIPT_BIN_DIR + "/Rscript " + R_SCRIPT_DIR + "/condensed_tophits_crossmodel.R "
				+ filteredFile + " " + filteredMalesFile + " " + filteredFemalesFile + " " + condensedFile + " "
				+ topHitsFile + " " + crossRanges + " " + pvaThreshold;

		listOfCommands.add(cmdToStore);

		try {
			GuidanceImpl.generateCondensedAndTopHitsFile(filteredFile, filteredMalesFile, filteredFemalesFile,
					condensedFile, topHitsFile, crossRanges, pvaThreshold, cmdToStore);
		} catch (GuidanceTaskException gte) {
			LOGGER.error("[Guidance] Exception trying the execution of generateCondensedFile task", gte);
		}

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
	 * @throws GuidanceTaskException
	 */
	private static void makePhenotypeAnalysis(ParseCmdLine parsingArgs, CombinedPanelsFiles combinedPanelsFilesInfo,
			ResultsFiles resultsFilesInfo, PhenomeAnalysisFiles phenomeAnalysisFilesInfo) throws GuidanceTaskException {

		if (parsingArgs.getStageStatus("phenoAnalysis") != 1) {
			return;
		}

		List<String> combinedTopHits = new ArrayList<>();

		int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();

		for (int test = 0; test < numberOfTestTypes; ++test) {
			combinedTopHits.add(combinedPanelsFilesInfo.getTopHitsFile(test));
		}

		String topHitsAllPheno = phenomeAnalysisFilesInfo.getTopHitsAllPhenos();

		String combinedTopHitsString = combinedTopHits.get(0);
		GuidanceImpl.getFile(combinedTopHits.get(0), combinedTopHits.get(0));
		for (int i = 1; i < combinedTopHits.size(); ++i) {
			combinedTopHitsString += ("," + combinedTopHits.get(i));
			GuidanceImpl.getFile(combinedTopHits.get(i), combinedTopHits.get(i));
		}
		COMPSs.barrier();
		String cmdToStore = R_SCRIPT_BIN_DIR + "/Rscript " + R_SCRIPT_DIR + "/tophits_all_phenotypes.R "
				+ combinedTopHitsString + " " + topHitsAllPheno;
		listOfCommands.add(cmdToStore);

		// This is a sequential invocation that implies bringing back all the combined
		// condensed files to the master
		GuidanceImpl.generateTopHitsAllPhenos(combinedTopHits, topHitsAllPheno);

		String condensedFile = null;
		String mergedPhenoFile = null;
		String pheno = null;
		for (int test = 0; test < numberOfTestTypes; ++test) {
			condensedFile = combinedPanelsFilesInfo.getCombinedCondensedFile(test);
			mergedPhenoFile = phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(test);
			pheno = parsingArgs.getTestTypeName(test);

			cmdToStore = R_SCRIPT_BIN_DIR + "/Rscript " + R_SCRIPT_DIR + "/merging_tophits_all_pheno.R "
					+ topHitsAllPheno + " " + condensedFile + " " + mergedPhenoFile + " " + pheno;
			listOfCommands.add(cmdToStore);

			// This is a task
			GuidanceImpl.generateMergedPhenoTopHits(topHitsAllPheno, condensedFile, mergedPhenoFile, pheno);
		}

		List<String> phenoMergedTopHits = new ArrayList<>();
		phenoMergedTopHits.add(phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(0));
		GuidanceImpl.getFile(phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(0),
				phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(0));
		String mergedTopHitsString = phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(0);
		for (int test = 1; test < numberOfTestTypes; ++test) {
			phenoMergedTopHits.add(phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(test));
			GuidanceImpl.getFile(phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(test),
					phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(test));
			mergedTopHitsString += ("," + phenomeAnalysisFilesInfo.getCrossPhenoMergedTop(test));
		}

		String pvaThreshold = Double.toString(parsingArgs.getPvaThreshold());

		String crossPhenoAll = phenomeAnalysisFilesInfo.getCrossPhenoAll();
		String crossPhenoRanges = phenomeAnalysisFilesInfo.getCrossPhenoRanges();
		String crossPhenoTopVariants = phenomeAnalysisFilesInfo.getCrossPhenoAssocTop();

		cmdToStore = R_SCRIPT_BIN_DIR + "/Rscript " + R_SCRIPT_DIR + "/crossphenotype.R " + mergedTopHitsString + " "
				+ crossPhenoAll + " " + crossPhenoRanges + " " + crossPhenoTopVariants + " " + pvaThreshold;
		listOfCommands.add(cmdToStore);
		COMPSs.barrier();
		// This is a sequential invocation that implies bringing back all the
		// phenoMerged files
		GuidanceImpl.computeCrossPheno(phenoMergedTopHits, crossPhenoAll, crossPhenoRanges, crossPhenoTopVariants,
				pvaThreshold);

		// FileUtils.getFile(crossPhenoAll);
		// FileUtils.getFile(crossPhenoRanges);
		// FileUtils.getFile(crossPhenoTopVariants);

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
			String basePath = mixedBedFile.substring(0, mixedBedFile.length() - 4);
			String cmdToStore = PLINKBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
					+ " --chr " + theChromo + " --out " + basePath + " --make-bed";
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
				myPrefix = bedChr23File.substring(0, bedChr23File.length() - 4);

				cmdToStore = PLINKBINARY + " --noweb --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
						+ " --filter-males --out " + myPrefix + " --make-bed";
			} else if (sex.equals(SEX2)) {
				myPrefix = bedChr23File.substring(0, bedChr23File.length() - 4);

				cmdToStore = PLINKBINARY + " --noweb --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
						+ " --filter-females --out " + myPrefix + " --make-bed";
			}

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.splitChr23(bedFile, bimFile, famFile, bedChr23File, bimChr23File, famChr23File, logFile,
						sex, theChromo, cmdToStore);
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
		String myPrefix = phasingSampleFile.split("\\.")[0];
		String cmd;

		if (parsingArgs.getStageStatus("phasingBed") == 1) {

			if (phasingTool.equals("shapeit")) {

				// If we process chromoso X (23) then we change the cmdToStore
				if (theChromo.equals("23")) {
					cmd = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map "
							+ gmapFile + " --chrX --output-max " + phasingHapsFile + " " + phasingSampleFile
							+ " --thread 48 --effective-size 20000 --output-log " + phasingLogFile;

				} else {
					cmd = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map "
							+ gmapFile + " --output-max " + phasingHapsFile + " " + phasingSampleFile
							+ " --thread 48 --effective-size 20000 --output-log " + phasingLogFile;
				}

				listOfCommands.add(new String(cmd));
				try {
					GuidanceImpl.phasingBed(theChromo, bedFile, bimFile, famFile, gmapFile, phasingHapsFile,
							phasingSampleFile, phasingLogFile, sex, phasingTool, cmd);
				} catch (Exception e) {
					System.err.println("[Guidance] Exception trying the execution of phasing task");
					System.err.println(e.getMessage());
				}

			} else if (phasingTool.equals("eagle")) {

				if (parsingArgs.getStageStatus("phasingBed") == 1) {
					// If we process chromoso X (23) then we change the cmdToStore
					if (theChromo.equals("23")) {

						if (sex.equals(SEX1)) {
							String baseDirOrigin = bedFile.substring(0, bedFile.length() - 4);
							String baseDirDest = phasingSampleFile.substring(0, phasingSampleFile.length() - 7);
							cmd = PLINKBINARY + " --bfile " + baseDirOrigin + " --recode vcf --out " + baseDirDest;
							listOfCommands.add(cmd);

							cmd = BCFTOOLSBINARY + " convert " + baseDirDest + ".vcf" + " --hapsample " + baseDirDest
									+ " --vcf-ids";
							listOfCommands.add(cmd);

							String generatedSample = phasingSampleFile + "s";

							cmd = "mv " + phasingSampleFile + " " + generatedSample + "; echo \"ID_1 ID_2 missing\" > "
									+ phasingSampleFile + "; echo \"0 0 0\" >> " + phasingSampleFile + "; tail -n +2 "
									+ generatedSample + " | tr \"_\" \" \" | awk '{ print $1\"_\"$2\" \"$3\" 0\" }' >> "
									+ phasingSampleFile;
						} else {
							cmd = EAGLEBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
									+ " --chrom " + theChromo + " --geneticMapFile " + gmapFile
									+ " --numThreads 48 --outPrefix " + myPrefix;
						}

					} else {
						cmd = EAGLEBINARY + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
								+ " --chrom " + theChromo + " --geneticMapFile " + gmapFile
								+ " --numThreads 48 --outPrefix " + myPrefix;
					}

					listOfCommands.add(cmd);
					try {
						GuidanceImpl.phasingBed(theChromo, bedFile, bimFile, famFile, gmapFile, phasingHapsFile,
								phasingSampleFile, phasingLogFile, sex, phasingTool, cmd);
					} catch (Exception e) {
						System.err.println("[Guidance] Exception trying the execution of phasing task");
						System.err.println(e.getMessage());
					}
				}
			} else {
				System.err.println("[Guidance] Only Eagle and Shapeit are available for phasing haplotypes");

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

			cmd = SAMTOOLSBINARY + "/bgzip " + filteredHaplotypesVcfFile;

			listOfCommands.add(new String(cmd));
			try {
				GuidanceImpl.filterHaplotypes(phasingHapsFile, phasingSampleFile, excludedSnpsFile,
						filteredHaplotypesLogFile, filteredHaplotypesVcfFileBgzip, cmd);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of filterHaplotypes task");
				System.err.println(e.getMessage());
			}

			// try {
			// GuidanceImpl.samtoolsBgzip(filteredHaplotypesVcfFile,
			// filteredHaplotypesVcfFileBgzip, cmd);
			// } catch (Exception e) {
			// System.err.println("[Guidance] Exception trying the execution of
			// samtoolsBgzip task");
			// System.err.println(e.getMessage());
			// }
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
			String imputeFileWarnings, String sex, int refpanel) {
		String cmdToStore = null;

		// TODO: remove when minimac es fixed
		if (parsingArgs.getStageStatus("imputeWithImpute") == 1) {
			// || ((parsingArgs.getStageStatus("imputeWithMinimac") == 1) &&
			// chrS.equals("23"))) {
			// Submitting the impute task per chunk

			if (chrS.equals("23")) {

				// if (sex.equals(SEX1)) {
				// cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " +
				// knownHapFile + " -l "
				// + legendFile + " -known_haps_g " + phasingHapsFile + " -sample_g " +
				// phasingSampleFile
				// + " -int " + lim1S + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile
				// + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo +
				// " -r "
				// + imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info
				// -o_gz ";
				// } else if (sex.equals(SEX2)) {
				cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
						+ legendFile + " -known_haps_g " + phasingHapsFile + " -sample_g " + phasingSampleFile
						+ " -int " + lim1S + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile
						+ " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r "
						+ imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz ";
				// }

			} else {
				cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
						+ legendFile + " -known_haps_g " + phasingHapsFile + " -int " + lim1S + " " + lim2S
						+ " -exclude_snps_g " + pairsFile + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i "
						+ imputeFileInfo + " -r " + imputeFileSummary + " -w " + imputeFileWarnings
						+ " -no_sample_qc_info -o_gz";
			}
			listOfCommands.add(cmdToStore);
			try {
				String panelMemory = parsingArgs.getRpanelMemory(refpanel);
				if (panelMemory.equals("HIGH")) {
					GuidanceImpl.imputeWithImputeHigh(gmapFile, knownHapFile, legendFile, phasingHapsFile,
							phasingSampleFile, lim1S, lim2S, pairsFile, imputeFile, imputeFileInfo, imputeFileSummary,
							imputeFileWarnings, chrS, sex, cmdToStore);
				} else if (panelMemory.equals("MEDIUM")) {
					GuidanceImpl.imputeWithImputeMedium(gmapFile, knownHapFile, legendFile, phasingHapsFile,
							phasingSampleFile, lim1S, lim2S, pairsFile, imputeFile, imputeFileInfo, imputeFileSummary,
							imputeFileWarnings, chrS, sex, cmdToStore);
				} else if (panelMemory.equals("LOW")) {
					GuidanceImpl.imputeWithImputeLow(gmapFile, knownHapFile, legendFile, phasingHapsFile,
							phasingSampleFile, lim1S, lim2S, pairsFile, imputeFile, imputeFileInfo, imputeFileSummary,
							imputeFileWarnings, chrS, sex, cmdToStore);
				} else {
					throw new GuidanceTaskException("Incorrect panel memory " + panelMemory);
				}

			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of impute task for chr " + chrS
						+ " with sex " + sex);
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
			String filteredHaplotypesVcfFileBgzip, String chrS, String lim1S, String lim2S, String imputeFileInfo,
			String imputeFileErate, String imputeFileRec, String imputeFileM3vcf, String imputeFileLog,
			String imputeFileBgzip, String imputeFileTbi, String sex, int refpanel) {

		if (parsingArgs.getStageStatus("imputeWithMinimac") == 1) {
			// Submitting the impute task per chunk
			// We don't distinguish chrS 23 since the cmdToStore is the same
			String cmdToStore = null;
			String myPrefix = imputeFileBgzip.split("\\.")[0];

			if (chrS.equals("23")) {
				// if (sex.equals(SEX1)) {
				// cmdToStore = MINIMAC_BINARY + " --refHaps " + refVcfFile + " --haps "
				// + filteredHaplotypesVcfFileBgzip + " --start " + lim1S + " --end " + lim2S
				// + " --chr X --window 500000 --prefix " + myPrefix
				// + " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
				// } else if (sex.equals(SEX2)) {
				// } else {
				cmdToStore = MINIMAC4_BINARY + " --refHaps " + refVcfFile + " --haps " + filteredHaplotypesVcfFileBgzip
						+ " --start " + lim1S + " --end " + lim2S + " --chr X --window 500000 --prefix " + myPrefix
						+ " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
				// }
			} else {
				cmdToStore = MINIMAC4_BINARY + " --refHaps " + refVcfFile + " --haps " + filteredHaplotypesVcfFileBgzip
						+ " --start " + lim1S + " --end " + lim2S + " --chr " + chrS + " --window 500000 --prefix "
						+ myPrefix + " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
			}

			listOfCommands.add(cmdToStore);
			try {
				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}

			try {
				String panelMemory = parsingArgs.getRpanelMemory(refpanel);
				if (panelMemory.equals("HIGH")) {
					GuidanceImpl.imputeWithMinimacHigh(refVcfFile, filteredHaplotypesVcfFileBgzip, imputeFileBgzip,
							imputeFileInfo, imputeFileErate, imputeFileRec, imputeFileM3vcf, imputeFileLog, chrS, lim1S,
							lim2S, myPrefix, sex, cmdToStore);
				} else if (panelMemory.equals("MEDIUM")) {
					GuidanceImpl.imputeWithMinimacMedium(refVcfFile, filteredHaplotypesVcfFileBgzip, imputeFileBgzip,
							imputeFileInfo, imputeFileErate, imputeFileRec, imputeFileM3vcf, imputeFileLog, chrS, lim1S,
							lim2S, myPrefix, sex, cmdToStore);
				} else if (panelMemory.equals("LOW")) {
					GuidanceImpl.imputeWithMinimacLow(refVcfFile, filteredHaplotypesVcfFileBgzip, imputeFileBgzip,
							imputeFileInfo, imputeFileErate, imputeFileRec, imputeFileM3vcf, imputeFileLog, chrS, lim1S,
							lim2S, myPrefix, sex, cmdToStore);
				} else {
					System.err.println("Incorrect panel memory " + panelMemory);
					throw new GuidanceTaskException("Incorrect panel memory " + panelMemory);
				}
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of impute task");
				System.err.println(e.getMessage());
			}

			/*
			 * 
			 * cmdToStore = SAMTOOLSBINARY + "/bgzip " + imputeFile + " " + imputeFileBgzip;
			 * 
			 * listOfCommands.add(cmdToStore); try { flushCommands(); } catch (IOException
			 * e) { LOGGER.error("[Guidance] Exception writing to list of commands file " +
			 * e); } try { GuidanceImpl.samtoolsBgzip(imputeFile, imputeFileBgzip,
			 * cmdToStore); } catch (Exception e) { LOGGER.
			 * error("[Guidance] Exception trying the execution of samtoolsBgzip task"); }
			 * 
			 */

			cmdToStore = SAMTOOLSBINARY + "/tabix -p vcf --force " + imputeFileBgzip + " " + imputeFileTbi;

			listOfCommands.add(cmdToStore);
			try {
				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}
			try {
				GuidanceImpl.samtoolsTabix(imputeFileBgzip, imputeFileTbi, cmdToStore);
			} catch (Exception e) {
				LOGGER.error("[Guidance] Exception trying the execution of samtoolsTabix task");
			}
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
			String chromo) {

		if (parsingArgs.getStageStatus("filterByInfo") == 1) {
			String infoThresholdS = null;
			String mafThresholdS = Double.toString(parsingArgs.getMafThreshold());

			String imputationTool = parsingArgs.getImputationTool();
			// if (chromo.equals("23")) {
			// imputationTool = "impute";
			// }

			if (imputationTool.equals("impute")) {
				infoThresholdS = Double.toString(parsingArgs.getImputeThreshold());
			} else if (imputationTool.equals("minimac")) {
				infoThresholdS = Double.toString(parsingArgs.getMinimacThreshold());
			}

			String cmdToStore = null;
			// We create the list of rsId that are greater than or equal to the
			// infoThreshold value
			cmdToStore = JAVA_HOME + "/java filterByInfo " + imputationTool + " " + imputeFileInfo + " "
					+ filteredRsIdFile + " " + infoThresholdS;
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.filterByInfo(imputationTool, imputeFileInfo, filteredRsIdFile, infoThresholdS,
						mafThresholdS, cmdToStore);
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
			String filteredFile, String filteredLogFile, String chromo) {

		double mafThreshold = parsingArgs.getMafThreshold();
		String mafThresholdS = Double.toString(mafThreshold);
		String imputationTool = parsingArgs.getImputationTool();

		if (parsingArgs.getStageStatus("qctoolS") == 1) {
			String cmdToStore = null;
			// if (imputationTool.equals("impute") || chromo.equals("23")) {
			if (imputationTool.equals("impute")) {

				imputationTool = "impute";
				cmdToStore = QCTOOL_BINARY + " -g " + imputeFile + " -og " + filteredFile + " -incl-rsids "
						+ filteredRsIdFile + " -omit-chromosome -force -log " + filteredLogFile + " -maf "
						+ mafThresholdS + " 1";

			} else if (imputationTool.equals("minimac")) {

				cmdToStore = QCTOOL_BINARY + " -g " + imputeFile + " -og " + filteredFile + " -incl-rsids "
						+ filteredRsIdFile + " -omit-chromosome -force -log " + filteredLogFile
						+ " -vcf-genotype-field GP";
			}

			listOfCommands.add(cmdToStore);
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
			//if (covariables.equals("none")) {
			//	cmdToStore = SNPTEST_BINARY + " -data " + mergedGenFile + " " + mergedSampleFile + " -o "
			//			+ snptestOutFile + " -pheno " + responseVar + " -hwe -log " + snptestLogFile;
			//} else {
			//	cmdToStore = SNPTEST_BINARY + " -data " + mergedGenFile + " " + mergedSampleFile + " -o "
			//			+ snptestOutFile + " -pheno " + responseVar + " -cov_names " + newStr + " -hwe -log "
			//			+ snptestLogFile;
			//}

			// Different parameters for chromo 23 (X) and the rest.
			//if (chrS.equals("23")) {
			//	cmdToStore = cmdToStore + " -method newml -assume_chromosome X -stratify_on sex -frequentist 1";
			//} else {
			//	cmdToStore = cmdToStore + " -method em -frequentist 1 2 3 4 5";
			//}
			
			cmdToStore = BASHSCRIPTDIR + "/snptest.sh " + SNPTEST_BINARY + " " + mergedGenFile + " " + mergedSampleFile
					+ " " + snptestOutFile + " " + snptestLogFile + " " + responseVar + " " + chrS + " " + covariables;

			listOfCommands.add(cmdToStore);
			try {
				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
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
			String snptestOutFile, String summaryFile, String mafThresholdS, String hweCohortThresholdS,
			String hweCasesThresholdS, String hweControlsThresholdS, String sex) {

		if (parsingArgs.getStageStatus("collectSummary") == 1) {
			// Submitting the collect_summary task per this chunk
			String cmdToStore = JAVA_HOME + "/java collectSummary " + chrS + " " + parsingArgs.getImputationTool() + " "
					+ imputeFileInfo + " " + snptestOutFile + " " + summaryFile + " " + mafThresholdS + " "
					+ hweCohortThresholdS + " " + hweCasesThresholdS + " " + hweControlsThresholdS;

			listOfCommands.add(cmdToStore);

			try {
				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}

			try {
				GuidanceImpl.collectSummary(chrS, parsingArgs.getImputationTool(), imputeFileInfo, snptestOutFile,
						summaryFile, mafThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, sex,
						cmdToStore);
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
	/*
	 * private static void doJointCondenseFiles(ParseCmdLine parsingArgs, String
	 * condensedA, String condensedB, String condensedC) {
	 * 
	 * if (parsingArgs.getStageStatus("jointCondensedFiles") == 1) { String
	 * cmdToStore = JAVA_HOME + "/java jointCondensedFiles " + condensedA + " " +
	 * condensedB + " " + condensedC; listOfCommands.add(cmdToStore); try {
	 * flushCommands(); } catch (IOException e) {
	 * LOGGER.error("[Guidance] Exception writing to list of commands file " + e); }
	 * try { GuidanceImpl.jointCondensedFiles(condensedA, condensedB, condensedC,
	 * cmdToStore); } catch (GuidanceTaskException gte) { LOGGER.
	 * error("[Guidance] Exception trying the execution of jointCondensedFiles task"
	 * , gte); } } }
	 */

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
				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
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

		String[] manhattans = parsingArgs.getManhattanOptions();
		String thresh = String.valueOf(parsingArgs.getPvaThreshold());
		String columnName = null;
		String manpdf = null;
		String mantiff = null;
		String qqpdf = null;
		String qqtiff = null;

		for (String option : manhattans) {
			manpdf = fileToInheritance(manhattanPlotFile, 4, option);
			mantiff = fileToInheritance(manhattanPlotTiffFile, 5, option);
			qqpdf = fileToInheritance(qqPlotFile, 4, option);
			qqtiff = fileToInheritance(qqPlotTiffFile, 5, option);
			columnName = "frequentist_" + option + "_pvalue";

			String cmdToStore = R_SCRIPT_BIN_DIR + "/Rscript " + R_SCRIPT_DIR + "/qqplot_manhattan_all_models.R "
					+ condensedFile + " " + qqpdf + " " + manpdf + " " + qqtiff + " " + mantiff + " " + columnName + " "
					+ thresh;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.generateQQManhattanPlots(condensedFile, qqpdf, manpdf, qqtiff, mantiff, columnName, thresh,
						cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of generateQQManhattanPlots task", gte);
			}

		}
	}

	private static String fileToInheritance(String filename, int length, String option) {
		return filename.substring(0, filename.length() - length) + "_" + option
				+ filename.substring(filename.length() - length, filename.length());
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
	/*
	 * private static void doMergeTwoChunks(ParseCmdLine parsingArgs, String
	 * reduceA, String reduceB, String reduceC) {
	 * 
	 * if (parsingArgs.getStageStatus("mergeTwoChunks") == 1) { // Task String
	 * cmdToStore = JAVA_HOME + "/java mergeTwoChunks " + reduceA + " " + reduceB +
	 * " " + reduceC; listOfCommands.add(cmdToStore); try { flushCommands(); } catch
	 * (IOException e) {
	 * LOGGER.error("[Guidance] Exception writing to list of commands file " + e); }
	 * try { GuidanceImpl.mergeTwoChunks(reduceA, reduceB, reduceC, cmdToStore); }
	 * catch (GuidanceTaskException gte) { LOGGER.
	 * error("[Guidance] Exception trying the execution of mergeTwoChunks task",
	 * gte); }
	 * 
	 * } }
	 */

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
			GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
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

		// Task
		if (parsingArgs.getStageStatus("filterByAll") == 1) {

			String imputationTool = parsingArgs.getImputationTool();
			double mafThreshold = parsingArgs.getMafThreshold();
			double infoThreshold = 0.0;
			double hweCohortThreshold = parsingArgs.getHweCohortThreshold();
			double hweCasesThreshold = parsingArgs.getHweCasesThreshold();
			double hweControlsThreshold = parsingArgs.getHweControlsThreshold();

			String mafThresholdS = Double.toString(mafThreshold);
			String infoThresholdS = Double.toString(infoThreshold);
			String hweCohortThresholdS = Double.toString(hweCohortThreshold);
			String hweCasesThresholdS = Double.toString(hweCasesThreshold);
			String hweControlsThresholdS = Double.toString(hweControlsThreshold);

			// Chr 23
			// if (sex.equals(SEX1) || sex.equals(SEX2)) {
			// imputationTool = "impute";
			// }

			if (imputationTool.equals("impute")) {
				infoThreshold = parsingArgs.getImputeThreshold();
			} else if (imputationTool.equals("minimac")) {
				infoThreshold = parsingArgs.getMinimacThreshold();
			}

			String cmdToStore = JAVA_HOME + "/java filterByAll " + imputationTool + " " + inputFile + " " + outputFile
					+ " " + mafThresholdS + " " + infoThresholdS + " " + hweCohortThresholdS + " " + hweCasesThresholdS
					+ " " + hweControlsThresholdS + " " + sex + " " + rpanelName;

			listOfCommands.add(cmdToStore);
			try {
				GeneralUtils.flushCommands(listOfStages, listOfCommands, FLUSH);
			} catch (IOException e) {
				LOGGER.error("[Guidance] Exception writing to list of commands file " + e);
			}

			try {
				GuidanceImpl.filterByAll(imputationTool, inputFile, outputFile, mafThresholdS, infoThresholdS,
						hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, sex, rpanelName, cmdToStore);
			} catch (Exception e) {
				System.err.println("[Guidance] Exception trying the execution of filterByAll task");
				System.err.println(e.getMessage());
			}
		}

	}

	private static void doCopyFile(ParseCmdLine parsingArgs, String originPath, String destinationPath)
			throws IOException, GuidanceTaskException {
		String cmdToStore = JAVA_HOME + "/java copyFile " + originPath + " " + destinationPath;

		listOfCommands.add(cmdToStore);

		GuidanceImpl.copyFile(originPath, destinationPath);
	}

	/**
	 * Method to verify that all the environment variables have been well defined
	 * 
	 */
	private static void verifyEnvVar() throws GuidanceEnvironmentException {
		GeneralUtils.verify(PLINK_BINARY);

		GeneralUtils.verifyEnvVarDefined(R_SCRIPT_BIN_DIR);
		GeneralUtils.verifyEnvVarDefined(R_SCRIPT_DIR);

		GeneralUtils.verify(QCTOOL_BINARY);
		GeneralUtils.verify(SHAPEIT_BINARY);
		GeneralUtils.verify(IMPUTE2_BINARY);
		GeneralUtils.verify(SNPTEST_BINARY);

		GeneralUtils.verifyEnvVarDefined(JAVA_HOME);
	}

}
