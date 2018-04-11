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
	private static final String SHAPEIT_BINARY = System.getenv("SHAPEITBINARY");
	private static final String IMPUTE2_BINARY = System.getenv("IMPUTE2BINARY");
	private static final String SNPTEST_BINARY = System.getenv("SNPTESTBINARY");
	private static final String JAVA_HOME = System.getenv("JAVA_HOME");

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
		ArrayList<String> listOfCommands = new ArrayList<>();

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
		File listOfStages = new File(listOfStagesFileName);
		if (!listOfStages.exists()) {
			LOGGER.info("\n[Guidance] File to store the tasks list: " + listOfStagesFileName);
		} else {
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
		doMixed(parsingArgs, outDir, rpanelTypes, listOfCommands);

		// Finally, we print the commands in the output file defined for this.
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(listOfStages))) {
			for (String str : listOfCommands) {
				writer.write(str);
				writer.newLine();
				writer.newLine();
			}

			// Close the file with the list of commands...
			writer.flush();
		}
		LOGGER.info("[Guidance] Everything is working with Guidance, just wait...");
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
	 */
	private static void doMixed(ParseCmdLine parsingArgs, String outDir, List<String> rpanelTypes,
			ArrayList<String> listOfCommands) throws IOException {
		// Create some general objects
		int startChr = parsingArgs.getStart();
		int endChr = parsingArgs.getEnd();
		// int endChrNormal = parsingArgs.getEndNormal();

		String exclCgatFlag = parsingArgs.getExclCgatSnp();
		String exclSVFlag = parsingArgs.getExclSVSnp();

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

		// Create all file names used in the workflow.
		// Also, define which of them will be temporal or permanent.
		setFinalStatusForCommonFiles(parsingArgs, commonFilesInfo);
		setFinalStatusForImputationFiles(parsingArgs, imputationFilesInfo, rpanelTypes);
		setFinalStatusForAssocFiles(parsingArgs, assocFilesInfo, rpanelTypes);

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

			String gmapFile = parsingArgs.getGmapDir() + "/" + parsingArgs.getGmapFileName(chr);
			String mixedShapeitHapsFile = commonFilesInfo.getShapeitHapsFile(chr);
			String mixedShapeitSampleFile = commonFilesInfo.getShapeitSampleFile(chr);
			String mixedShapeitLogFile = commonFilesInfo.getShapeitLogFile(chr);
			String mixedExcludedSnpsFile = commonFilesInfo.getExcludedSnpsFile(chr);

			String mixedFilteredHaplotypesFile = commonFilesInfo.getFilteredHaplotypesFile(chr);
			String mixedFilteredHaplotypesSampleFile = commonFilesInfo.getFilteredHaplotypesSampleFile(chr);
			String mixedFilteredHaplotypesLogFile = commonFilesInfo.getFilteredHaplotypesLogFile(chr);
			String mixedFilteredHaplotypesVcfFile = commonFilesInfo.getFilteredHaplotypesVcfFile(chr);
			String mixedFilteredListOfSnpsFile = commonFilesInfo.getListOfSnpsFile(chr);

			// Check if the input are in GEN o BED format
			inputFormat = parsingArgs.getInputFormat();
			if (inputFormat.equals("BED")) {
				bedFile = commonFilesInfo.getBedFile();
				bimFile = commonFilesInfo.getBimFile();
				famFile = commonFilesInfo.getFamFile();

				mixedBedFile = commonFilesInfo.getByChrBedFile(chr);
				mixedBimFile = commonFilesInfo.getByChrBimFile(chr);
				mixedFamFile = commonFilesInfo.getByChrFamFile(chr);
				mixedBedToBedLogFile = commonFilesInfo.getBedToBedLogFile(chr);

				doConvertFromBedToBed(parsingArgs, listOfCommands, bedFile, bimFile, famFile, mixedBedFile,
						mixedBimFile, mixedFamFile, mixedBedToBedLogFile, theChromo);

				// Create the RsId list of SNPs that are AT, TA, CG, or GC
				// In that case, because inputType is BED we pass the newBimFile
				doCreateRsIdList(parsingArgs, listOfCommands, mixedBimFile, exclCgatFlag, mixedPairsFile, inputFormat);

				doPhasingBed(parsingArgs, listOfCommands, theChromo, mixedBedFile, mixedBimFile, mixedFamFile, gmapFile,
						mixedExcludedSnpsFile, mixedShapeitHapsFile, mixedShapeitSampleFile, mixedShapeitLogFile,
						mixedFilteredHaplotypesFile, mixedFilteredHaplotypesSampleFile, mixedFilteredHaplotypesLogFile,
						mixedFilteredHaplotypesVcfFile, mixedFilteredListOfSnpsFile, exclCgatFlag, exclSVFlag);

			} else { // The inputFormat is GEN
				mixedGenFile = commonFilesInfo.getGenFile(chr);
				// Task
				doCreateRsIdList(parsingArgs, listOfCommands, mixedGenFile, exclCgatFlag, mixedPairsFile, inputFormat);

				doPhasing(parsingArgs, listOfCommands, theChromo, mixedGenFile, mixedSampleFile, gmapFile,
						mixedExcludedSnpsFile, mixedShapeitHapsFile, mixedShapeitSampleFile, mixedShapeitLogFile,
						mixedFilteredHaplotypesFile, mixedFilteredHaplotypesSampleFile, mixedFilteredHaplotypesLogFile,
						mixedFilteredHaplotypesVcfFile, mixedFilteredListOfSnpsFile, exclCgatFlag, exclSVFlag);

			} // End of inputFormat GEN

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
					makeImputationPerChunk(parsingArgs, chr, lim1, lim2, panel, gmapFile, imputationFilesInfo,
							commonFilesInfo, listOfCommands);
					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				}
			} // End for panel types

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
							String panelName = rpanelTypes.get(panel);
							LOGGER.debug("Making association for TEST " + testName + ", PANEL " + panelName
									+ ", CHROMO " + chr + ", CHUNK " + lim1 + " - " + lim2);
						}
						makeAssociationPerChunk(parsingArgs, test, panel, chr, lim1, lim2, imputationFilesInfo,
								commonFilesInfo, assocFilesInfo, listOfCommands);

						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}

					// Now we perform the merge of chunks for each chromosome
					makeMergeOfChunks(parsingArgs, listOfCommands, test, panel, chr, minSize, maxSize, chunkSize,
							assocFilesInfo, mergeFilesInfo, FILTERED);
					makeMergeOfChunks(parsingArgs, listOfCommands, test, panel, chr, minSize, maxSize, chunkSize,
							assocFilesInfo, mergeFilesInfo, CONDENSED);
				} // End for Chromo

				// Now we have to joint the condensedFiles of each chromosome. There is not
				// problem if chr23 is being
				// processed, because
				// the format of condensedFiles is the same for all chromosome.
				makeJointCondensedFiles(parsingArgs, listOfCommands, test, panel, startChr, endChr, mergeFilesInfo);

				// Now we have to joint the filteredByAllFiles of each chromosome. Here there is
				// an additional
				// complexity due to the chr23,
				// because the file format changes if chr23 is being processed. This situation
				// is taken into account
				// inside the next function.
				String rpanelName = rpanelTypes.get(panel);
				makeJointFilteredByAllFiles(parsingArgs, listOfCommands, test, panel, rpanelName, startChr, endChr,
						mergeFilesInfo);

				// Generate TopHits
				String lastCondensedFile = mergeFilesInfo.getFinalCondensedFile(test, panel);
				String lastFilteredByAllFile = mergeFilesInfo.getFinalFilteredByAllFile(test, panel);
				String topHitsResults = resultsFilesInfo.getTopHitsFile(test, panel);
				String filteredByAllXFile;
				if (endChr == 23) {
					filteredByAllXFile = mergeFilesInfo.getAdditionalFilteredByAllXFile(test, panel, 0);
				} else {
					filteredByAllXFile = lastFilteredByAllFile;
				}

				doGenerateTopHits(parsingArgs, listOfCommands, lastFilteredByAllFile, filteredByAllXFile,
						topHitsResults, PVA_THRESHOLD_STR);

				// Generate QQManhattan Plots
				String qqPlotPdfFile = resultsFilesInfo.getQqPlotPdfFile(test, panel);
				String qqPlotTiffFile = resultsFilesInfo.getQqPlotTiffFile(test, panel);
				String manhattanPlotPdfFile = resultsFilesInfo.getManhattanPdfFile(test, panel);
				String manhattanPlotTiffFile = resultsFilesInfo.getManhattanTiffFile(test, panel);
				// String correctedPvaluesFile = resultsFilesInfo.getCorrectedPvaluesFile(test,
				// panel);

				doGenerateQQManhattanPlots(parsingArgs, listOfCommands, lastCondensedFile, qqPlotPdfFile,
						manhattanPlotPdfFile, qqPlotTiffFile, manhattanPlotTiffFile);

			} // End for refPanels

			// Now we continue with the combining of the results of the different reference
			// panels.
			// It is done if the refPanelCombine flag is true.
			makeCombinePanels(parsingArgs, assocFilesInfo, mergeFilesInfo, combinedPanelsFilesInfo, rpanelTypes, test,
					listOfCommands);

		} // End for test types

		if (1 < numberOfTestTypes) {
			makePhenotypeAnalysis(parsingArgs, mergeFilesInfo, resultsFilesInfo, phenomeAnalysisFilesInfo, rpanelTypes,
					listOfCommands);
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
			int panelIndex, String gmapFile, ImputationFiles imputationFilesInfo, CommonFiles commonFilesInfo,
			ArrayList<String> listOfCommands) {

		int chunkSize = parsingArgs.getChunkSize();
		double infoThreshold = parsingArgs.getInfoThreshold();

		String knownHapFileName = parsingArgs.getRpanelHapFileName(panelIndex, chrNumber);
		String rpanelDir = parsingArgs.getRpanelDir(panelIndex);
		String knownHapFile = rpanelDir + File.separator + knownHapFileName;

		String lim1S = Integer.toString(lim1);
		String lim2S = Integer.toString(lim2);
		String chrS = Integer.toString(chrNumber);
		String infoThresholdS = Double.toString(infoThreshold);
		String imputationTool = parsingArgs.getImputationTool();

		if (imputationTool.equals("impute")) {
			String legendFileName = parsingArgs.getRpanelLegFileName(panelIndex, chrNumber);
			String legendFile = rpanelDir + File.separator + legendFileName;

			// String mixedSampleFile = commonFilesInfo.getSampleFile(chrNumber);
			String mixedShapeitHapsFile = commonFilesInfo.getShapeitHapsFile(chrNumber);
			String mixedShapeitSampleFile = commonFilesInfo.getShapeitSampleFile(chrNumber);
			String mixedPairsFile = commonFilesInfo.getPairsFile(chrNumber);
			String mixedImputeFile = imputationFilesInfo.getImputedFile(panelIndex, chrNumber, lim1, lim2, chunkSize);
			String mixedImputeFileInfo = imputationFilesInfo.getImputedInfoFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeFileSummary = imputationFilesInfo.getImputedSummaryFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputeFileWarnings = imputationFilesInfo.getImputedWarningsFile(panelIndex, chrNumber, lim1,
					lim2, chunkSize);
			// String mixedImputeLogFile = imputationFilesInfo.getImputedLogFile(panelIndex,
			// chrNumber, lim1, lim2,
			// chunkSize);

			String mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedFilteredLogFile = imputationFilesInfo.getFilteredLogFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);

			// We create the list of rsId that are greater than or equal to the
			// infoThreshold value
			String mixedFilteredRsIdFile = imputationFilesInfo.getFilteredRsIdFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);

			doImputationWithImpute(parsingArgs, listOfCommands, chrS, gmapFile, knownHapFile, legendFile,
					mixedShapeitHapsFile, mixedShapeitSampleFile, lim1S, lim2S, mixedPairsFile, mixedImputeFile,
					mixedImputeFileInfo, mixedImputeFileSummary, mixedImputeFileWarnings);

			doFilterByInfo(parsingArgs, listOfCommands, mixedImputeFileInfo, mixedFilteredRsIdFile, infoThresholdS);

			doQctoolS(parsingArgs, listOfCommands, mixedImputeFile, mixedFilteredRsIdFile, mixedFilteredFile,
					mixedFilteredLogFile);
		} else if (imputationTool.equals("minimac")) {
			String mixedFilteredHaplotypesFile = commonFilesInfo.getFilteredHaplotypesFile(chrNumber);
			String mixedFilteredHaplotypesSampleFile = commonFilesInfo.getFilteredHaplotypesSampleFile(chrNumber);
			String mixedFilteredListOfSnpsFile = commonFilesInfo.getListOfSnpsFile(chrNumber);

			String mixedImputedMMFileName = imputationFilesInfo.getImputedMMFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputedMMInfoFile = imputationFilesInfo.getImputedMMInfoFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			// String mixedImputedMMDraftFile =
			// imputationFilesInfo.getImputedMMDraftFile(panelIndex, chrNumber, lim1,
			// lim2, chunkSize);
			String mixedImputedMMErateFile = imputationFilesInfo.getImputedMMErateFile(panelIndex, chrNumber, lim1,
					lim2, chunkSize);
			String mixedImputedMMRecFile = imputationFilesInfo.getImputedMMRecFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputedMMDoseFile = imputationFilesInfo.getImputedMMDoseFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);
			String mixedImputedMMLogFile = imputationFilesInfo.getImputedMMLogFile(panelIndex, chrNumber, lim1, lim2,
					chunkSize);

			doImputationWithMinimac(parsingArgs, listOfCommands, knownHapFile, mixedFilteredHaplotypesFile,
					mixedFilteredHaplotypesSampleFile, mixedFilteredListOfSnpsFile, mixedImputedMMFileName,
					mixedImputedMMInfoFile, mixedImputedMMErateFile, mixedImputedMMRecFile, mixedImputedMMDoseFile,
					mixedImputedMMLogFile, chrS, lim1S, lim2S);
		} else {
			LOGGER.fatal(
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
			int chrNumber, int lim1, int lim2, ImputationFiles imputationFilesInfo, CommonFiles commonFilesInfo,
			AssocFiles assocFilesInfo, ArrayList<String> listOfCommands) {

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

		String snptestOutFile = assocFilesInfo.getSnptestOutFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
				chunkSize);
		String snptestLogFile = assocFilesInfo.getSnptestLogFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
				chunkSize);

		String mixedSampleFile = commonFilesInfo.getSampleFile(chrNumber);
		String mixedImputeFileInfo = imputationFilesInfo.getImputedInfoFile(panelIndex, chrNumber, lim1, lim2,
				chunkSize);

		String mixedFilteredFile = imputationFilesInfo.getFilteredFile(panelIndex, chrNumber, lim1, lim2, chunkSize);

		// String inputFormat = parsingArgs.getInputFormat();
		// String exclCgatFlag = parsingArgs.getExclCgatSnp();

		doSnptest(parsingArgs, listOfCommands, chrS, mixedFilteredFile, mixedSampleFile, snptestOutFile, snptestLogFile,
				responseVar, covariables);

		String summaryFile = assocFilesInfo.getSummaryFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2, chunkSize);

		doCollectSummary(parsingArgs, listOfCommands, chrS, mixedImputeFileInfo, snptestOutFile, summaryFile,
				mafThresholdS, infoThresholdS, hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);

		String assocFilteredByAll = assocFilesInfo.getSummaryFilteredFile(testTypeIndex, panelIndex, chrNumber, lim1,
				lim2, chunkSize);
		String assocCondensed = assocFilesInfo.getSummaryCondensedFile(testTypeIndex, panelIndex, chrNumber, lim1, lim2,
				chunkSize);

		doFilterByAll(parsingArgs, listOfCommands, summaryFile, assocFilteredByAll, assocCondensed);
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
	private static void makeMergeOfChunks(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, int ttIndex,
			int rpanelIndex, int chr, int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo,
			MergeFiles mergeFilesInfo, String type) {

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

		String filteredByAllFile = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, chr);
		String condensedFile = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, chr);

		if (type.equals(FILTERED)) {
			// LOGGER.info("Number of chunks for testType " + ttIndex + " | rpanel " +
			// rpanelIndex + " |chr " +
			// chr + " " + numberOfChunks);
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
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, filteredByAllFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getReducedFilteredFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, reducedC);
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
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, condensedFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getReducedCondensedFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, reducedC);
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
	private static void makeMergeOfChunksCombined(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			int ttIndex, int rpanelIndex, int chr, int minSize, int maxSize, int chunkSize, AssocFiles assocFilesInfo,
			MergeFiles mergeFilesInfo, String type) {

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

		String filteredByAllFile = mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, rpanelIndex, chr);
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
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, filteredByAllFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedFilteredFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, reducedC);
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
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, condensedFile);
					indexC++;
				} else {
					reducedC = mergeFilesInfo.getCombinedReducedCondensedFile(ttIndex, rpanelIndex, chr, indexC);
					doMergeTwoChunks(parsingArgs, listOfCommands, reducedA, reducedB, reducedC);
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
	 */
	private static void makeJointCondensedFiles(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, int ttIndex,
			int rpanelIndex, int startChr, int endChr, MergeFiles mergeFilesInfo) {

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
			condensedB = condensedA;

			condensedC = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexC);
			doJointCondenseFiles(parsingArgs, listOfCommands, condensedA, condensedB, condensedC);
		} else {
			for (int processedCondensed = 0; processedCondensed < 2 * numberOfChrs
					- 2; processedCondensed = processedCondensed + 2) {
				if (processedCondensed < numberOfChrs) {
					i = startChr + processedCondensed;
					condensedA = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, i);
				} else {
					condensedA = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}

				if (processedCondensed < numberOfChrs - 1) {
					i = startChr + processedCondensed + 1;
					condensedB = mergeFilesInfo.getCondensedFile(ttIndex, rpanelIndex, i);
				} else {
					condensedB = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexA);
					indexA++;
				}
				condensedC = mergeFilesInfo.getAdditionalCondensedFile(ttIndex, rpanelIndex, indexC);

				doJointCondenseFiles(parsingArgs, listOfCommands, condensedA, condensedB, condensedC);
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
	private static void makeJointFilteredByAllFiles(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			int ttIndex, int rpanelIndex, String rpanelName, int startChr, int endChr, MergeFiles mergeFilesInfo) {

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
				doJointFilteredByAllFiles(parsingArgs, listOfCommands, filteredA, filteredB, filteredC, rpanelName,
						rpanelFlag);
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

				doJointFilteredByAllFiles(parsingArgs, listOfCommands, filteredA, filteredB, filteredC, rpanelName,
						rpanelFlag);
				indexC++;
			} // End for fitlered chromosomes
				// processedFiltered +2)
		}

		// Now we process the chr 23 if this is defined in this execution
		if (endChr == 23) {
			rpanelFlag = "YES";
			filteredA = mergeFilesInfo.getFilteredByAllFile(ttIndex, rpanelIndex, endChr);
			filteredC = mergeFilesInfo.getAdditionalFilteredByAllXFile(ttIndex, rpanelIndex, 0);

			doJointFilteredByAllFiles(parsingArgs, listOfCommands, filteredA, filteredA, filteredC, rpanelName,
					rpanelFlag);
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
			CombinedPanelsFiles combinedPanelsFilesInfo, List<String> rpanelTypes, int ttIndex,
			ArrayList<String> listOfCommands) throws IOException {

		final boolean refPanelCombine = parsingArgs.getRefPanelCombine();
		if (!refPanelCombine) {
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
		String tmpFCAX = filteredCombineAll;
		if (endChr == 23) {
			tmpFCAX = combinedPanelsFilesInfo.getCombinedFilteredByAllXFile(ttIndex);
		}
		final String filteredCombineAllX = tmpFCAX;

		final String condensedCombineAll = combinedPanelsFilesInfo.getCombinedCondensedFile(ttIndex);

		// INITIALIZE THE FILES TO STORE ALL COMBINED INFORMATION (ADD HEADER AND
		// COMPRESS)
		final String filteredHeader = Headers.constructHeader();
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

		String filteredXHeader = null;
		if (endChr == 23) {
			filteredXHeader = Headers.constructHeaderX();
			final String plainfilteredCombineAllX = filteredCombineAllX.substring(0, filteredCombineAllX.length() - 3);
			try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainfilteredCombineAllX))) {
				writer.write(filteredXHeader);
				writer.newLine();
				writer.flush();
			} catch (IOException ioe) {
				LOGGER.error("[Guidance] Exception when initializing makeCombinePanel filteredX ALL file", ioe);
			}
			FileUtils.gzipFile(plainfilteredCombineAllX, filteredCombineAllX);
			new File(plainfilteredCombineAllX).delete();
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

		LinkedList<String> filteredCombined = new LinkedList<>();
		LinkedList<String> condensedCombined = new LinkedList<>();

		// CHR LOOP
		for (int chr = startChr; chr <= endChr; chr++) {
			int minSize = ChromoInfo.getMinSize(chr);
			int maxSize = ChromoInfo.getMaxSize(chr);

			int lim1 = minSize;
			int lim2 = lim1 + chunkSize - 1;
			int indexFC = 0;
			int indexCC = 0;
			int indexXFC = 0;

			// CHUNK LOOP
			for (int j = minSize; j < maxSize; j = j + chunkSize) {
				// -- FILTERED PART --
				// Construct a queue with all the filtered panels to combine
				LinkedList<String> filteredPanelsToCombine = new LinkedList<>();
				for (int k = 0; k < rpanelTypes.size(); ++k) {
					String filteredPanel = assocFilesInfo.getSummaryFilteredFile(ttIndex, k, chr, lim1, lim2,
							chunkSize);
					if (DEBUG) {
						LOGGER.debug("[Guidance] Adding filtered file to combine " + filteredPanel);
					}
					filteredPanelsToCombine.add(filteredPanel);
				}
				// Combine all the filtered panels 2 by 2 until there are no remaining panels
				while (!filteredPanelsToCombine.isEmpty()) {
					String filteredPanelA = filteredPanelsToCombine.poll();
					if (!filteredPanelsToCombine.isEmpty()) {
						String filteredPanelB = filteredPanelsToCombine.poll();
						// Filtered part: combines A and B into A
						if (DEBUG) {
							LOGGER.debug("[Guidance] Combining " + filteredPanelA + " and " + filteredPanelB);
						}
						doCombinePanelsComplex(parsingArgs, listOfCommands, filteredPanelA, filteredPanelB, lim1, lim2);
						// Adds A to the queue again
						filteredPanelsToCombine.add(filteredPanelA);

						// Deletes B since it is no longer needed
						// new File(filteredPanelB).delete();
					} else {
						String destFilteredPanel = assocFilesInfo.getCombinedFilteredFile(ttIndex, 0, chr, lim1, lim2,
								chunkSize);
						GuidanceImpl.copyFile(filteredPanelA, destFilteredPanel);
					}
				}

				// -- CONDENSED PART --
				// Construct a queue with all the condensed panels to combine
				LinkedList<String> condensedPanelsToCombine = new LinkedList<>();
				for (int k = 0; k < rpanelTypes.size(); ++k) {
					String condensedPanel = assocFilesInfo.getSummaryCondensedFile(ttIndex, k, chr, lim1, lim2,
							chunkSize);
					if (DEBUG) {
						LOGGER.debug("[Guidance] Adding condensed file to combine " + condensedPanel);
					}
					condensedPanelsToCombine.add(condensedPanel);
				}
				// Combine all the filtered panels 2 by 2 until there are no remaining panels
				while (!condensedPanelsToCombine.isEmpty()) {
					String condensedPanelA = condensedPanelsToCombine.poll();
					if (!condensedPanelsToCombine.isEmpty()) {
						String condensedPanelB = condensedPanelsToCombine.poll();
						// Filtered part: combines A and B into A
						if (DEBUG) {
							LOGGER.debug("[Guidance] Combining " + condensedPanelA + " and " + condensedPanelB);
						}
						doCombinePanelsComplex(parsingArgs, listOfCommands, condensedPanelA, condensedPanelB, lim1,
								lim2);
						// Adds A to the queue again
						condensedPanelsToCombine.add(condensedPanelA);

						// Deletes B since it is no longer needed
						// new File(condensedPanelB).delete();
					} else {
						String destCondensedPanel = assocFilesInfo.getCombinedCondensedFile(ttIndex, 0, chr, lim1, lim2,
								chunkSize);
						GuidanceImpl.copyFile(condensedPanelA, destCondensedPanel);
					}
				}

				// -- COMBINE FILTERED TO FINAL CHUNK AND ALL FILES
				// Files are merged to the first or the last file depending on the panel Size
				int mergeIndex = (rpanelTypes.size() % 2 == 0) ? 0 : rpanelTypes.size() - 1;
				String chunkResultsFiltered = assocFilesInfo.getSummaryFilteredFile(ttIndex, mergeIndex, chr, lim1,
						lim2, chunkSize);
				if (chr != 23) {
					// We merge in the chunk files
					String filteredCombinePerChunk = combinedPanelsFilesInfo.getCombinedFilteredByAllChromoFile(ttIndex,
							indexFC);
					if (!new File(filteredCombinePerChunk).exists()) {
						// Initialize filtered per chunk file to default header and compress it
						String plainfilteredCombinePerChunk = filteredCombinePerChunk.substring(0,
								filteredCombinePerChunk.length() - 3);
						try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainfilteredCombinePerChunk))) {
							writer.write(filteredHeader);
							writer.newLine();
							writer.flush();
						} catch (IOException ioe) {
							LOGGER.error(
									"[Guidance] Exception when initializing makeCombinePanel filteredCombinePerChunk file",
									ioe);
						}
						FileUtils.gzipFile(plainfilteredCombinePerChunk, filteredCombinePerChunk);
						new File(plainfilteredCombinePerChunk).delete();
					}
					/*
					 * if (DEBUG) { LOGGER.debug( "[Guidance] Merging filtered into chunk-file: " +
					 * filteredCombinePerChunk + " and " + chunkResultsFiltered); }
					 * doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands,
					 * filteredCombinePerChunk, chunkResultsFiltered, Integer.toString(chr),
					 * FILTERED);
					 */
					// We merge in the all files
					if (DEBUG) {
						LOGGER.debug("[Guidance] Merging filtered into ALL: " + filteredCombineAll + " and "
								+ chunkResultsFiltered);
					}
					/*
					 * doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands, filteredCombineAll,
					 * chunkResultsFiltered, Integer.toString(chr), FILTERED);
					 */
					++indexFC;
				} else {
					// We merge in the chunk files
					String filteredCombineXPerChunk = combinedPanelsFilesInfo
							.getCombinedFilteredByAllXChromoFile(ttIndex, indexXFC);
					if (!new File(filteredCombineXPerChunk).exists()) {
						// Initialize filtered per chunk file to default header and compress it
						String plainfilteredCombineXPerChunk = filteredCombineXPerChunk.substring(0,
								filteredCombineXPerChunk.length() - 3);
						try (BufferedWriter writer = new BufferedWriter(
								new FileWriter(plainfilteredCombineXPerChunk))) {
							writer.write(filteredXHeader);
							writer.newLine();
							writer.flush();
						} catch (IOException ioe) {
							LOGGER.error(
									"[Guidance] Exception when initializing makeCombinePanel filteredCombineXPerChunk file",
									ioe);
						}
						FileUtils.gzipFile(plainfilteredCombineXPerChunk, filteredCombineXPerChunk);
						new File(plainfilteredCombineXPerChunk).delete();
					}

					// We merge in the all files
					/*
					 * doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands, filteredCombineAllX,
					 * chunkResultsFiltered, Integer.toString(chr), FILTERED);
					 */
					++indexXFC;
				}

				// Clean partial results
				// new File(chunkResultsFiltered).delete();

				// -- COMBINE CONDENSED TO FINAL CHUNK AND ALL FILES
				String chunkResultsCondensed = assocFilesInfo.getSummaryCondensedFile(ttIndex, mergeIndex, chr, lim1,
						lim2, chunkSize);
				// We merge in the chunk files
				String condensedCombinePerChunk = combinedPanelsFilesInfo.getCombinedCondensedChromoFile(ttIndex,
						indexCC);
				if (!new File(condensedCombinePerChunk).exists()) {
					// Initialize filtered per chunk file to default header and compress it
					String plainCondensedCombinePerChunk = condensedCombinePerChunk.substring(0,
							condensedCombinePerChunk.length() - 3);
					try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainCondensedCombinePerChunk))) {
						writer.write(condensedHeader);
						writer.newLine();
						writer.flush();
					} catch (IOException ioe) {
						LOGGER.error(
								"[Guidance] Exception when initializing makeCombinePanel plainCondensedCombinePerChunk file",
								ioe);
					}
					FileUtils.gzipFile(plainCondensedCombinePerChunk, condensedCombinePerChunk);
					new File(plainCondensedCombinePerChunk).delete();
				}
				/*
				 * if (DEBUG) { LOGGER.debug( "[Guidance] Merging combined into chunk-file: " +
				 * condensedCombinePerChunk + " and " + chunkResultsCondensed); }
				 * doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands,
				 * condensedCombinePerChunk, chunkResultsCondensed, Integer.toString(chr),
				 * CONDENSED);
				 */

				// We merge in the all files
				if (DEBUG) {
					LOGGER.debug("[Guidance] Merging combined into ALL: " + condensedCombineAll + " and "
							+ chunkResultsCondensed);
				}
				/*
				 * doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands, condensedCombineAll,
				 * chunkResultsCondensed, Integer.toString(chr), CONDENSED);
				 */
				++indexCC;

				// Clean partial results
				// new File(chunkResultsCondensed).delete();

				// Increase loop variables
				lim1 = lim1 + chunkSize;
				lim2 = lim2 + chunkSize;
			} // End for chunk
			makeMergeOfChunksCombined(parsingArgs, listOfCommands, ttIndex, 0, chr, minSize, maxSize, chunkSize,
					assocFilesInfo, mergeFilesInfo, FILTERED);

			makeMergeOfChunksCombined(parsingArgs, listOfCommands, ttIndex, 0, chr, minSize, maxSize, chunkSize,
					assocFilesInfo, mergeFilesInfo, CONDENSED);

			String filteredByAllCurrentFile = mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, 0, chr);
			String condensedCurrentFile = mergeFilesInfo.getCombinedCondensedFile(ttIndex, 0, chr);

			filteredCombined.add(filteredByAllCurrentFile);
			condensedCombined.add(condensedCurrentFile);

		} // End for chromosomes

		int reduceCounter = 0;

		while (filteredCombined.size() > 1) {
			String originFilteredFileA = filteredCombined.poll();
			String originFilteredFileB = filteredCombined.poll();
			String originCondensedFileA = condensedCombined.poll();
			String originCondensedFileB = condensedCombined.poll();

			String destinationFilteredFile;
			String destinationCondensedFile;

			if (condensedCombined.isEmpty()) {
				destinationFilteredFile = filteredCombineAll;
				destinationCondensedFile = condensedCombineAll;
			} else {
				destinationFilteredFile = filteredCombineAll.substring(0, filteredCombineAll.length() - 7) + "_reduce_"
						+ Integer.toString(reduceCounter) + ".txt.gz";
				destinationCondensedFile = condensedCombineAll.substring(0, condensedCombineAll.length() - 7)
						+ "_reduce_" + Integer.toString(reduceCounter) + ".txt.gz";
			}

			doMergeTwoChunks(parsingArgs, listOfCommands, originFilteredFileA, originFilteredFileB,
					destinationFilteredFile);
			filteredCombined.add(destinationFilteredFile);

			doMergeTwoChunks(parsingArgs, listOfCommands, originCondensedFileA, originCondensedFileB,
					destinationCondensedFile);
			condensedCombined.add(destinationCondensedFile);

			reduceCounter += 1;
		}
		/*
		 * for (int chr = startChr; chr <= endChr; chr++) { String filteredFile =
		 * mergeFilesInfo.getCombinedFilteredByAllFile(ttIndex, 0, chr);
		 * doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands, filteredCombineAll,
		 * filteredFile, Integer.toString(chr), FILTERED);
		 * 
		 * String condensedFile = mergeFilesInfo.getCombinedCondensedFile(ttIndex, 0,
		 * chr); doMergeTwoChunksInTheFirst(parsingArgs, listOfCommands,
		 * condensedCombineAll, condensedFile, Integer.toString(chr), CONDENSED); }
		 */

		// Finally, we create topHits from filteredCombined, and QQ and Manhattan plots
		// from condensedCombined
		String topHitsCombinedResults = combinedPanelsFilesInfo.getTopHitsFile(ttIndex);
		doGenerateTopHits(parsingArgs, listOfCommands, filteredCombineAll, filteredCombineAllX, topHitsCombinedResults,
				PVA_THRESHOLD_STR);

		String combinedQqPlotPdfFile = combinedPanelsFilesInfo.getQqPlotPdfFile(ttIndex);
		String combinedQqPlotTiffFile = combinedPanelsFilesInfo.getQqPlotTiffFile(ttIndex);
		String combinedManhattanPlotPdfFile = combinedPanelsFilesInfo.getManhattanPdfFile(ttIndex);
		String combinedManhattanPlotTiffFile = combinedPanelsFilesInfo.getManhattanTiffFile(ttIndex);
		// String combinedCorrectedPvaluesFile =
		// combinedPanelsFilesInfo.getCorrectedPvaluesFile(ttIndex);

		doGenerateQQManhattanPlots(parsingArgs, listOfCommands, condensedCombineAll, combinedQqPlotPdfFile,
				combinedManhattanPlotPdfFile, combinedQqPlotTiffFile, combinedManhattanPlotTiffFile);

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
			ResultsFiles resultsFilesInfo, PhenomeAnalysisFiles phenomeAnalysisFilesInfo, List<String> rpanelTypes,
			ArrayList<String> listOfCommands) {

		int endChr = parsingArgs.getEnd();
		String endChrS = Integer.toString(endChr);

		String topHitsFile = null;

		String filteredByAllFile = null;
		String filteredByAllXFile = null;

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

		doInitPhenoMatrix(parsingArgs, listOfCommands, topHitsFile, ttName, rpName, phenomeFileA);

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

				doAddToPhenoMatrix(parsingArgs, listOfCommands, phenomeFileA, topHitsFile, ttName, rpName,
						phenomeFileB);

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
					filteredByAllXFile = mergeFilesInfo.getAdditionalFilteredByAllXFile(ttIndex, rpIndex, 0);
				} else {
					filteredByAllXFile = mergeFilesInfo.getFinalFilteredByAllFile(ttIndex, rpIndex);
				}

				phenomeFileB = phenomeAnalysisFilesInfo.getPhenotypeFile(phenoIndex);

				doFilloutPhenoMatrix(parsingArgs, listOfCommands, phenomeFileA, filteredByAllFile, filteredByAllXFile,
						endChrS, ttName, rpName, phenomeFileB);

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

				doFinalizePhenoMatrix(parsingArgs, listOfCommands, phenomeFileA, phenomeFileB, ttName, rpName,
						phenomeFileC);

				phenomeFileA = phenomeFileC;
				phenoIndex++;
			}
		}

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
	private static void doConvertFromBedToBed(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String bedFile, String bimFile, String famFile, String mixedBedFile, String mixedBimFile,
			String mixedFamFile, String mixedBedToBedLogFile, String theChromo) {

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
	private static void doCreateRsIdList(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String mixedBimOrGenFile, String exclCgatFlag, String mixedPairsFile, String inputFormat) {

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
	private static void doPhasingBed(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String theChromo,
			String bedFile, String bimFile, String famFile, String gmapFile, String excludedSnpsFile,
			String shapeitHapsFile, String shapeitSampleFile, String shapeitLogFile, String filteredHaplotypesFile,
			String filteredHaplotypesSampleFile, String filteredHaplotypesLogFile, String filteredHaplotypesVcfFile,
			String listOfSnpsFile, String exclCgatFlag, String exclSVFlag) {

		if (parsingArgs.getStageStatus("phasingBed") == 1) {
			// If we process chromoso X (23) then we change the cmdToStore
			String cmdToStore;
			if (theChromo.equals("23")) {
				cmdToStore = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile
						+ " --input-map " + gmapFile + " --chrX --output-max " + shapeitHapsFile + " "
						+ shapeitSampleFile + " --thread 47 --effective-size 20000 --output-log " + shapeitLogFile;
			} else {
				cmdToStore = SHAPEIT_BINARY + " --input-bed " + bedFile + " " + bimFile + " " + famFile
						+ " --input-map " + gmapFile + " --output-max " + shapeitHapsFile + " " + shapeitSampleFile
						+ " --thread 47 --effective-size 20000 --output-log " + shapeitLogFile;
			}

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.phasingBed(theChromo, bedFile, bimFile, famFile, gmapFile, shapeitHapsFile,
						shapeitSampleFile, shapeitLogFile, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of phasing task", gte);
			}
		}

		if (parsingArgs.getStageStatus("createListOfExcludedSnps") == 1) {
			String cmdToStore = JAVA_HOME + "/java createListOfExcludedSnps.jar " + shapeitHapsFile + " "
					+ excludedSnpsFile + " " + exclCgatFlag + " " + exclSVFlag;

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.createListOfExcludedSnps(shapeitHapsFile, excludedSnpsFile, exclCgatFlag, exclSVFlag,
						cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of createListOfExcludedSnps task", gte);
			}
		}

		if (parsingArgs.getStageStatus("filterHaplotypes") == 1) {
			String cmdToStore = JAVA_HOME + "/java filterHaplotypes.jar " + shapeitHapsFile + " " + shapeitSampleFile
					+ " " + excludedSnpsFile + " " + filteredHaplotypesFile + " " + filteredHaplotypesSampleFile + " "
					+ filteredHaplotypesLogFile + " " + filteredHaplotypesVcfFile + " " + listOfSnpsFile;

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.filterHaplotypes(shapeitHapsFile, shapeitSampleFile, excludedSnpsFile,
						filteredHaplotypesFile, filteredHaplotypesSampleFile, filteredHaplotypesLogFile,
						filteredHaplotypesVcfFile, listOfSnpsFile, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of createListOfExcludedSnps task", gte);
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
	private static void doPhasing(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String theChromo,
			String genFile, String sampleFile, String gmapFile, String excludedSnpsFile, String shapeitHapsFile,
			String shapeitSampleFile, String shapeitLogFile, String filteredHaplotypesFile,
			String filteredHaplotypesSampleFile, String filteredHaplotypesLogFile, String filteredHaplotypesVcfFile,
			String listOfSnpsFile, String exclCgatFlag, String exclSVFlag) {

		if (parsingArgs.getStageStatus("pashing") == 1) {
			String cmdToStore;
			if (theChromo.equals("23")) {
				// If we process chromoso X (23) then we change the cmdToStore
				cmdToStore = SHAPEIT_BINARY + " --input-gen " + genFile + " " + sampleFile + " --input-map " + gmapFile
						+ " --chrX --output-max " + shapeitHapsFile + " " + shapeitSampleFile
						+ " --thread 15 --effective-size 20000 --output-log " + shapeitLogFile;
			} else {
				cmdToStore = SHAPEIT_BINARY + " --input-gen " + genFile + " " + sampleFile + " --input-map " + gmapFile
						+ " --output-max " + shapeitHapsFile + " " + shapeitSampleFile
						+ " --thread 15 --effective-size 20000 --output-log " + shapeitLogFile;
			}
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.phasing(theChromo, genFile, sampleFile, gmapFile, shapeitHapsFile, shapeitSampleFile,
						shapeitLogFile, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of phasing task", gte);
			}
		}

		if (parsingArgs.getStageStatus("createListOfExcludedSnps") == 1) {
			String cmdToStore = JAVA_HOME + "/java createListOfExcludedSnps.jar " + shapeitHapsFile + " "
					+ excludedSnpsFile + " " + exclCgatFlag + " " + exclSVFlag;

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.createListOfExcludedSnps(shapeitHapsFile, excludedSnpsFile, exclCgatFlag, exclSVFlag,
						cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of createListOfExcludedSnps task", gte);
			}
		}

		if (parsingArgs.getStageStatus("filterHaplotypes") == 1) {
			String cmdToStore = JAVA_HOME + "/java filterHaplotypes.jar " + shapeitHapsFile + " " + shapeitSampleFile
					+ " " + excludedSnpsFile + " " + filteredHaplotypesFile + " " + filteredHaplotypesSampleFile + " "
					+ filteredHaplotypesLogFile + " " + filteredHaplotypesVcfFile + " " + listOfSnpsFile;

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.filterHaplotypes(shapeitHapsFile, shapeitSampleFile, excludedSnpsFile,
						filteredHaplotypesFile, filteredHaplotypesSampleFile, filteredHaplotypesLogFile,
						filteredHaplotypesVcfFile, listOfSnpsFile, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of createListOfExcludedSnps task", gte);
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
	private static void doImputationWithImpute(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String chrS,
			String gmapFile, String knownHapFile, String legendFile, String shapeitHapsFile, String shapeitSampleFile,
			String lim1S, String lim2S, String pairsFile, String imputeFile, String imputeFileInfo,
			String imputeFileSummary, String imputeFileWarnings) {

		if (parsingArgs.getStageStatus("imputeWithImpute") == 1) {
			// Submitting the impute task per chunk
			String cmdToStore;
			if (chrS.equals("23")) {
				cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
						+ legendFile + " -known_haps_g " + shapeitHapsFile + " -sample_g " + shapeitSampleFile
						+ " -int " + lim1S + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile
						+ " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r "
						+ imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz ";
			} else {
				cmdToStore = IMPUTE2_BINARY + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l "
						+ legendFile + " -known_haps_g " + shapeitHapsFile + " -int " + lim1S + " " + lim2S
						+ " -exclude_snps_g " + pairsFile + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i "
						+ imputeFileInfo + " -r " + imputeFileSummary + " -w " + imputeFileWarnings
						+ " -no_sample_qc_info -o_gz";
			}
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.imputeWithImpute(gmapFile, knownHapFile, legendFile, shapeitHapsFile, shapeitSampleFile,
						lim1S, lim2S, pairsFile, imputeFile, imputeFileInfo, imputeFileSummary, imputeFileWarnings,
						chrS, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of impute task", gte);
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
	private static void doImputationWithMinimac(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String knownHapFile, String filteredHapsFile, String filteredSampleFile, String filteredListOfSnpsFile,
			String imputedMMFileName, String imputedMMInfoFile, String imputedMMErateFile, String imputedMMRecFile,
			String imputedMMDoseFile, String imputedMMLogFile, String chrS, String lim1S, String lim2S) {

		if (parsingArgs.getStageStatus("imputeWithMinimac") == 1) {
			// Submitting the impute task per chunk
			// We don't distinguish chrS 23 since the cmdToStore is the same
			String cmdToStore = JAVA_HOME + "/java imputationWithMinimac --vcfReference --refHaps " + knownHapFile
					+ " --shape_haps " + filteredHapsFile + " --sample " + filteredSampleFile + " --snps "
					+ filteredListOfSnpsFile + " --vcfstart " + lim1S + " --vcfend " + lim2S + " --chr " + chrS
					+ " --vcfwindow --rounds 5 --states 200 --outInfo " + imputedMMInfoFile + " --outErate "
					+ imputedMMErateFile + " --outRec " + imputedMMRecFile + " --outDose " + imputedMMDoseFile
					+ " --outLog " + imputedMMLogFile;

			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.imputeWithMinimac(knownHapFile, filteredHapsFile, filteredSampleFile,
						filteredListOfSnpsFile, imputedMMFileName, imputedMMInfoFile, imputedMMErateFile,
						imputedMMRecFile, imputedMMDoseFile, imputedMMLogFile, chrS, lim1S, lim2S, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of imputationWithMinimac task", gte);
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
	private static void doFilterByInfo(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String imputeFileInfo, String filteredRsIdFile, String infoThresholdS) {

		if (parsingArgs.getStageStatus("filterByInfo") == 1) {
			// We create the list of rsId that are greater than or equal to the
			// infoThreshold value
			String cmdToStore = JAVA_HOME + "/java filterByInfo " + imputeFileInfo + " " + filteredRsIdFile + " "
					+ infoThresholdS;
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.filterByInfo(imputeFileInfo, filteredRsIdFile, infoThresholdS, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of filterByInfo tasks for controls", gte);
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
	private static void doQctoolS(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String imputeFile,
			String filteredRsIdFile, String filteredFile, String filteredLogFile) {

		double mafThreshold = parsingArgs.getMafThreshold();
		String mafThresholdS = Double.toString(mafThreshold);

		if (parsingArgs.getStageStatus("qctoolS") == 1) {
			String cmdToStore = QCTOOL_BINARY + " -g " + imputeFile + " -og " + filteredFile + " -incl-rsids "
					+ filteredRsIdFile + " -omit-chromosome -force -log " + filteredLogFile + " -maf " + mafThresholdS
					+ " 1";
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.qctoolS(imputeFile, filteredRsIdFile, mafThresholdS, filteredFile, filteredLogFile,
						cmdToStore);
			} catch (GuidanceTaskException gte) {
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
	private static void doSnptest(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String chrS,
			String mergedGenFile, String mergedSampleFile, String snptestOutFile, String snptestLogFile,
			String responseVar, String covariables) {

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
	private static void doCollectSummary(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String chrS,
			String imputeFileInfo, String snptestOutFile, String summaryFile, String mafThresholdS,
			String infoThresholdS, String hweCohortThresholdS, String hweCasesThresholdS,
			String hweControlsThresholdS) {

		if (parsingArgs.getStageStatus("collectSummary") == 1) {
			// Submitting the collect_summary task per this chunk
			String cmdToStore = JAVA_HOME + "/java collectSummary " + chrS + " " + imputeFileInfo + " " + snptestOutFile
					+ " " + summaryFile + " " + mafThresholdS + " " + infoThresholdS + " " + hweCohortThresholdS + " "
					+ hweCasesThresholdS + " " + hweControlsThresholdS;

			listOfCommands.add(cmdToStore);

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
	private static void doJointCondenseFiles(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String condensedA, String condensedB, String condensedC) {

		if (parsingArgs.getStageStatus("jointCondensedFiles") == 1) {
			String cmdToStore = JAVA_HOME + "/java jointCondensedFiles " + condensedA + " " + condensedB + " "
					+ condensedC;
			listOfCommands.add(cmdToStore);

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
	private static void doJointFilteredByAllFiles(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String filteredByAllA, String filteredByAllB, String filteredByAllC, String rpanelName, String rpanelFlag) {

		if (parsingArgs.getStageStatus("jointFilteredByAllFiles") == 1) {
			String cmdToStore = JAVA_HOME + "/java jointFilteredByAllFiles " + filteredByAllA + " " + filteredByAllB
					+ " " + filteredByAllC + " " + rpanelName + " " + rpanelFlag;
			listOfCommands.add(cmdToStore);
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
	private static void doGenerateTopHits(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String filteredFile, String filteredXFile, String topHitsResults, String pvaThrS) {

		if (parsingArgs.getStageStatus("generateTopHits") == 1) {
			String cmdToStore = JAVA_HOME + "/java generateTopHits " + filteredFile + " " + filteredXFile + " "
					+ topHitsResults + " " + pvaThrS;
			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.generateTopHitsAll(filteredFile, filteredXFile, topHitsResults, pvaThrS, cmdToStore);
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
	private static void doGenerateQQManhattanPlots(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String condensedFile, String qqPlotFile, String manhattanPlotFile, String qqPlotTiffFile,
			String manhattanPlotTiffFile) {

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
	private static void doCombinePanelsComplex(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String resultsPanelA, String resultsPanelB, int lim1, int lim2) {

		if (parsingArgs.getStageStatus("combinePanelsComplex") == 1) {
			String cmdToStore = JAVA_HOME + "/java combinePanelsComplex " + resultsPanelA + " " + resultsPanelB + " "
					+ lim1 + " " + lim2;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.combinePanelsComplex(resultsPanelA, resultsPanelB, lim1, lim2, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of combinePanelsComplex task", gte);
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
	private static void doMergeTwoChunks(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String reduceA,
			String reduceB, String reduceC) {

		if (parsingArgs.getStageStatus("mergeTwoChunks") == 1) {
			// Task
			String cmdToStore = JAVA_HOME + "/java mergeTwoChunks " + reduceA + " " + reduceB + " " + reduceC;
			listOfCommands.add(cmdToStore);

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
	private static void doMergeTwoChunksInTheFirst(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String reduceA, String reduceB, String theChromo, String type) {

		if (parsingArgs.getStageStatus("mergeTwoChunks") == 1) {
			// Task
			String cmdToStore = JAVA_HOME + "/java mergeTwoChunksInTheFirst " + reduceA + " " + reduceB + " "
					+ theChromo;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.mergeTwoChunksInTheFirst(reduceA, reduceB, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of mergeTwoChunks task", gte);
			}

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
	private static void doFilterByAll(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands, String inputFile,
			String outputFile, String outputCondensedFile) {

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
			String cmdToStore = JAVA_HOME + "/java filterByAll " + inputFile + " " + outputFile + " "
					+ outputCondensedFile + " " + mafThresholdS + " " + infoThresholdS + " " + hweCohortThresholdS + " "
					+ hweCasesThresholdS + " " + hweControlsThresholdS;

			listOfCommands.add(cmdToStore);
			try {
				GuidanceImpl.filterByAll(inputFile, outputFile, outputCondensedFile, mafThresholdS, infoThresholdS,
						hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, cmdToStore);
			} catch (GuidanceTaskException gte) {
				LOGGER.error("[Guidance] Exception trying the execution of filterByAll task", gte);
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
	private static void doInitPhenoMatrix(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String topHitsFile, String ttName, String rpName, String phenomeFile) {

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
	private static void doAddToPhenoMatrix(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String phenomeFileA, String topHitsFile, String ttName, String rpName, String phenomeFileB) {

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
	private static void doFilloutPhenoMatrix(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String phenomeFileA, String filteredByAllFile, String filteredByAllXFile, String endChrS, String ttName,
			String rpName, String phenomeFileB) {

		if (parsingArgs.getStageStatus("filloutPhenoMatrix") == 1) {
			String cmdToStore = JAVA_HOME + "/java filloutPhenoMatrix " + phenomeFileA + " " + filteredByAllFile + " "
					+ filteredByAllXFile + " " + endChrS + " " + ttName + " " + rpName + " " + phenomeFileB;
			listOfCommands.add(cmdToStore);

			try {
				GuidanceImpl.filloutPhenoMatrix(phenomeFileA, filteredByAllFile, filteredByAllXFile, endChrS, ttName,
						rpName, phenomeFileB, cmdToStore);
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
	private static void doFinalizePhenoMatrix(ParseCmdLine parsingArgs, ArrayList<String> listOfCommands,
			String phenomeFileA, String phenomeFileB, String ttName, String rpName, String phenomeFileC) {

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
			commonFilesInfo.setShapeitHapsFileFinalStatus(i, finalStatus);
			commonFilesInfo.setShapeitSampleFileFinalStatus(i, finalStatus);
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
						imputationFilesInfo.setImputedFileFinalStatus(j, chromo, lim1, lim2, chunkSize, finalStatus);
						imputationFilesInfo.setFilteredFileFinalStatus(j, chromo, lim1, lim2, chunkSize, finalStatus);
						imputationFilesInfo.setImputedInfoFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
								finalStatus);

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
						imputationFilesInfo.setImputedMMInfoFileFinalStatus(j, chromo, lim1, lim2, chunkSize,
								finalStatus);
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
						assocFilesInfo.setSnptestOutFileFinalStatus(tt, j, chromo, lim1, lim2, chunkSize, finalStatus);
						assocFilesInfo.setSummaryFileFinalStatus(tt, j, chromo, lim1, lim2, chunkSize, finalStatus);

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
