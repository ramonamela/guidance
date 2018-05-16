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

package guidance.files;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import guidance.utils.ChromoInfo;
import guidance.utils.ParseCmdLine;

import java.io.File;
import java.util.ArrayList;

public class ImputationFiles {

	// Logger
	private static final Logger LOGGER = LogManager.getLogger("Console");

	private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;
	private static final int MIN_LIMIT = 1;
	private static final int MAX_LIMIT = 252_000_000;

	private int startChr;
	private int endChr;

	// Two dimensional array for the names of chunks:
	private ArrayList<ArrayList<String>> imputedOutDir = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedInfoFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedInfoMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedInfoFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedSummaryFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedSummaryMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedSummaryFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedWarningsFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedWarningsMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedWarningsFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedLogFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedLogMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedLogFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFileLogFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredLogMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredLogFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFileRsIdFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredRsIdMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredRsIdFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<String>> imputedMMOutDir = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMDoseVCFFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMDoseVCFMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMDoseVCFFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedFileTbi = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMalesFileTbi = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedFemalesFileTbi = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedFileBgzip = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMalesFileBgzip = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedFemalesFileBgzip = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMInfoFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMInfoMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMInfoFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMErateFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMErateMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMErateFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMRecFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMRecMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMRecFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFFileBgzip = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFMalesFileBgzip = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFFemalesFileBgzip = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFFileTbi = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFMalesFileTbi = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMM3VCFFemalesFileTbi = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMLogFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMLogMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMLogFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

	/**
	 * A constructor for the class
	 * 
	 * @param parsingArgs
	 * @param generalChromoInfo
	 * @param myOutDir
	 * @param refPanels
	 */
	public ImputationFiles(ParseCmdLine parsingArgs, String myOutDir, List<String> refPanels) {
		this.startChr = parsingArgs.getStart();
		this.endChr = parsingArgs.getEnd();

		String imputationTool = parsingArgs.getImputationTool();
		if (imputationTool.equals("impute")) {
			initializeForImpute(parsingArgs, myOutDir, refPanels);
		} else if (imputationTool.equals("minimac")) {
			initializeForMinimac(parsingArgs, myOutDir, refPanels);
		} else {
			System.err.println(
					"[ImputationFiles] Error, this imputation tool (" + imputationTool + ") is not supported yet!.");
			System.exit(1);
		}
	}

	/**
	 * Initialize file structure for impute
	 * 
	 * @param parsingArgs
	 * @param generalChromoInfo
	 * @param myOutDir
	 * @param refPanels
	 */
	private void initializeForImpute(ParseCmdLine parsingArgs, String myOutDir, List<String> refPanels) {
		int chunkSize = parsingArgs.getChunkSize();

		// We create the first directory name: the cohort directory.
		String cohort = parsingArgs.getCohort();
		String tmpOutDir = myOutDir + File.separator + cohort;

		// Now I create the directories for imputedOutDir
		for (int j = 0; j < refPanels.size(); j++) {
			String rPanel = refPanels.get(j);
			String tmpPanelDir = tmpOutDir + File.separator + rPanel;
			// Next level: Create directories.
			String mixOutDir = tmpPanelDir + File.separator + "mixed";

			ArrayList<String> chromoListImputedOutDir = new ArrayList<String>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedInfoFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedSummaryFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedWarningsFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedLogFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredLogFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdFile = new ArrayList<ArrayList<GenericFile>>();

			ArrayList<ArrayList<GenericFile>> chromoListImputedMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedInfoMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedSummaryMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedWarningsMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedLogMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredLogMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdMalesFile = new ArrayList<ArrayList<GenericFile>>();

			ArrayList<ArrayList<GenericFile>> chromoListImputedFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedInfoFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedSummaryFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedWarningsFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedLogFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredLogFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdFemalesFile = new ArrayList<ArrayList<GenericFile>>();

			// int maxSize = chromoInformation.getMaxSize(i);
			for (int i = this.startChr; i <= this.endChr; i++) {
				int chromo = i;
				int lim1 = 1;
				int lim2 = lim1 + chunkSize - 1;
				String tmpChrDir = mixOutDir + File.separator + "Chr_" + chromo;
				chromoListImputedOutDir.add(tmpChrDir);

				ArrayList<GenericFile> chunkListImputedFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedInfoFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedSummaryFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedWarningsFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedLogFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredLogFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredRsIdFile = new ArrayList<GenericFile>();

				ArrayList<GenericFile> chunkListImputedMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedInfoMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedSummaryMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedWarningsMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedLogMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredLogMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredRsIdMalesFile = new ArrayList<GenericFile>();

				ArrayList<GenericFile> chunkListImputedFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedInfoFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedSummaryFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedWarningsFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedLogFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredLogFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredRsIdFemalesFile = new ArrayList<GenericFile>();

				int numberOfChunks = ChromoInfo.getMaxSize(chromo) / chunkSize;
				int module = ChromoInfo.getMaxSize(chromo) % chunkSize;
				if (module != 0)
					numberOfChunks++;

				for (int k = 0; k < numberOfChunks; k++) {
					if (chromo < 23) {

						String imputedFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ ".impute";
						GenericFile myChunkListImputedFile = new GenericFile(tmpChrDir, imputedFileName + ".gz",
								"compressed", "none");
						chunkListImputedFile.add(myChunkListImputedFile);

						GenericFile myChunkListImputedInfoFile = new GenericFile(tmpChrDir, imputedFileName + "_info",
								"compressed", "none");
						chunkListImputedInfoFile.add(myChunkListImputedInfoFile);

						GenericFile myChunkListImputedSummaryFile = new GenericFile(tmpChrDir,
								imputedFileName + "_summary", "decompressed", "none");
						chunkListImputedSummaryFile.add(myChunkListImputedSummaryFile);

						GenericFile myChunkListImputedWarningsFile = new GenericFile(tmpChrDir,
								imputedFileName + "_warnings", "decompressed", "none");
						chunkListImputedWarningsFile.add(myChunkListImputedWarningsFile);

						GenericFile myChunkListImputedLogFile = new GenericFile(tmpChrDir, imputedFileName + ".log",
								"decompressed", "none");
						chunkListImputedLogFile.add(myChunkListImputedLogFile);

						String filteredFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_filtered.impute.gz";
						GenericFile myChunkListFilteredFile = new GenericFile(tmpChrDir, filteredFileName, "compressed",
								"none");
						chunkListFilteredFile.add(myChunkListFilteredFile);

						String filteredLogFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_filtered.impute.log";
						GenericFile myChunkListFilteredLogFile = new GenericFile(tmpChrDir, filteredLogFileName,
								"decompressed", "none");
						chunkListFilteredLogFile.add(myChunkListFilteredLogFile);

						String filteredRsIdFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_filtered_rsid.txt";
						GenericFile myChunkListFilteredRsIdFile = new GenericFile(tmpChrDir, filteredRsIdFileName,
								"compressed", "none");
						chunkListFilteredRsIdFile.add(myChunkListFilteredRsIdFile);

					} else if (i == 23) {
						String imputedFileMalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_males.impute";
						GenericFile myChunkListImputedMalesFile = new GenericFile(tmpChrDir,
								imputedFileMalesName + ".gz", "compressed", "none");
						chunkListImputedMalesFile.add(myChunkListImputedMalesFile);

						GenericFile myChunkListImputedInfoMalesFile = new GenericFile(tmpChrDir,
								imputedFileMalesName + "_info", "compressed", "none");
						chunkListImputedInfoMalesFile.add(myChunkListImputedInfoMalesFile);

						GenericFile myChunkListImputedSummaryMalesFile = new GenericFile(tmpChrDir,
								imputedFileMalesName + "_summary", "decompressed", "none");
						chunkListImputedSummaryMalesFile.add(myChunkListImputedSummaryMalesFile);

						GenericFile myChunkListImputedWarningsMalesFile = new GenericFile(tmpChrDir,
								imputedFileMalesName + "_warnings", "decompressed", "none");
						chunkListImputedWarningsMalesFile.add(myChunkListImputedWarningsMalesFile);

						GenericFile myChunkListImputedLogMalesFile = new GenericFile(tmpChrDir,
								imputedFileMalesName + ".log", "decompressed", "none");
						chunkListImputedLogMalesFile.add(myChunkListImputedLogMalesFile);

						String filteredFileMalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_males_filtered.impute.gz";
						GenericFile myChunkListFilteredMalesFile = new GenericFile(tmpChrDir, filteredFileMalesName,
								"compressed", "none");
						chunkListFilteredMalesFile.add(myChunkListFilteredMalesFile);

						String filteredLogFileMalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_males_filtered.impute.log";
						GenericFile myChunkListFilteredLogMalesFile = new GenericFile(tmpChrDir,
								filteredLogFileMalesName, "decompressed", "none");
						chunkListFilteredLogMalesFile.add(myChunkListFilteredLogMalesFile);

						String filteredRsIdFileMalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_"
								+ lim2 + "_males_filtered_rsid.txt";
						GenericFile myChunkListFilteredRsIdMalesFile = new GenericFile(tmpChrDir,
								filteredRsIdFileMalesName, "compressed", "none");
						chunkListFilteredRsIdMalesFile.add(myChunkListFilteredRsIdMalesFile);

						String imputedFileFemalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_females.impute";
						GenericFile myChunkListImputedFemalesFile = new GenericFile(tmpChrDir,
								imputedFileFemalesName + ".gz", "compressed", "none");
						chunkListImputedFemalesFile.add(myChunkListImputedFemalesFile);

						GenericFile myChunkListImputedInfoFemalesFile = new GenericFile(tmpChrDir,
								imputedFileFemalesName + "_info", "compressed", "none");
						chunkListImputedInfoFemalesFile.add(myChunkListImputedInfoFemalesFile);

						GenericFile myChunkListImputedSummaryFemalesFile = new GenericFile(tmpChrDir,
								imputedFileFemalesName + "_summary", "decompressed", "none");
						chunkListImputedSummaryFemalesFile.add(myChunkListImputedSummaryFemalesFile);

						GenericFile myChunkListImputedWarningsFemalesFile = new GenericFile(tmpChrDir,
								imputedFileFemalesName + "_warnings", "decompressed", "none");
						chunkListImputedWarningsFemalesFile.add(myChunkListImputedWarningsFemalesFile);

						GenericFile myChunkListImputedLogFemalesFile = new GenericFile(tmpChrDir,
								imputedFileFemalesName + ".log", "decompressed", "none");
						chunkListImputedLogFemalesFile.add(myChunkListImputedLogFemalesFile);

						String filteredFileFemalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_females_filtered.impute.gz";
						GenericFile myChunkListFilteredFemalesFile = new GenericFile(tmpChrDir, filteredFileFemalesName,
								"compressed", "none");
						chunkListFilteredFemalesFile.add(myChunkListFilteredFemalesFile);

						String filteredLogFileFemalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_"
								+ lim2 + "_females_filtered.impute.log";
						GenericFile myChunkListFilteredLogFemalesFile = new GenericFile(tmpChrDir,
								filteredLogFileFemalesName, "decompressed", "none");
						chunkListFilteredLogFemalesFile.add(myChunkListFilteredLogFemalesFile);

						String filteredRsIdFileFemalesName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_"
								+ lim2 + "_females_filtered_rsid.txt";
						GenericFile myChunkListFilteredRsIdFemalesFile = new GenericFile(tmpChrDir,
								filteredRsIdFileFemalesName, "compressed", "none");
						chunkListFilteredRsIdFemalesFile.add(myChunkListFilteredRsIdFemalesFile);

					}

					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				}

				if (chromo < 23) {
					chromoListImputedFile.add(chunkListImputedFile);
					chromoListImputedInfoFile.add(chunkListImputedInfoFile);
					chromoListImputedSummaryFile.add(chunkListImputedSummaryFile);
					chromoListImputedWarningsFile.add(chunkListImputedWarningsFile);
					chromoListImputedLogFile.add(chunkListImputedLogFile);
					chromoListFilteredFile.add(chunkListFilteredFile);
					chromoListFilteredLogFile.add(chunkListFilteredLogFile);
					chromoListFilteredRsIdFile.add(chunkListFilteredRsIdFile);
				} else if (chromo == 23) {
					chromoListImputedMalesFile.add(chunkListImputedMalesFile);
					chromoListImputedInfoMalesFile.add(chunkListImputedInfoMalesFile);
					chromoListImputedSummaryMalesFile.add(chunkListImputedSummaryMalesFile);
					chromoListImputedWarningsMalesFile.add(chunkListImputedWarningsMalesFile);
					chromoListImputedLogMalesFile.add(chunkListImputedLogMalesFile);
					chromoListFilteredMalesFile.add(chunkListFilteredMalesFile);
					chromoListFilteredLogMalesFile.add(chunkListFilteredLogMalesFile);
					chromoListFilteredRsIdMalesFile.add(chunkListFilteredRsIdMalesFile);

					chromoListImputedFemalesFile.add(chunkListImputedFemalesFile);
					chromoListImputedInfoFemalesFile.add(chunkListImputedInfoFemalesFile);
					chromoListImputedSummaryFemalesFile.add(chunkListImputedSummaryFemalesFile);
					chromoListImputedWarningsFemalesFile.add(chunkListImputedWarningsFemalesFile);
					chromoListImputedLogFemalesFile.add(chunkListImputedLogFemalesFile);
					chromoListFilteredFemalesFile.add(chunkListFilteredFemalesFile);
					chromoListFilteredLogFemalesFile.add(chunkListFilteredLogFemalesFile);
					chromoListFilteredRsIdFemalesFile.add(chunkListFilteredRsIdFemalesFile);
				}

				this.imputedOutDir.add(chromoListImputedOutDir);
				if (chromo < 23) {
					imputedFile.add(chromoListImputedFile);

					imputedInfoFile.add(chromoListImputedInfoFile);
					imputedSummaryFile.add(chromoListImputedSummaryFile);
					imputedWarningsFile.add(chromoListImputedWarningsFile);
					imputedLogFile.add(chromoListImputedLogFile);
					filteredFile.add(chromoListFilteredFile);
					filteredFileLogFile.add(chromoListFilteredLogFile);
					filteredFileRsIdFile.add(chromoListFilteredRsIdFile);
				} else if (chromo == 23) {
					imputedMalesFile.add(chromoListImputedMalesFile);

					imputedInfoMalesFile.add(chromoListImputedInfoMalesFile);
					imputedSummaryMalesFile.add(chromoListImputedSummaryMalesFile);
					imputedWarningsMalesFile.add(chromoListImputedWarningsMalesFile);
					imputedLogMalesFile.add(chromoListImputedLogMalesFile);
					filteredMalesFile.add(chromoListFilteredMalesFile);
					filteredLogMalesFile.add(chromoListFilteredLogMalesFile);
					filteredRsIdMalesFile.add(chromoListFilteredRsIdMalesFile);

					imputedFemalesFile.add(chromoListImputedFemalesFile);

					imputedInfoFemalesFile.add(chromoListImputedInfoFemalesFile);
					imputedSummaryFemalesFile.add(chromoListImputedSummaryFemalesFile);
					imputedWarningsFemalesFile.add(chromoListImputedWarningsFemalesFile);
					imputedLogFemalesFile.add(chromoListImputedLogFemalesFile);
					filteredFemalesFile.add(chromoListFilteredFemalesFile);
					filteredLogFemalesFile.add(chromoListFilteredLogFemalesFile);
					filteredRsIdFemalesFile.add(chromoListFilteredRsIdFemalesFile);
				}
			}
		}
	}

	/**
	 * Initialize file structure for Minimac
	 * 
	 * @param parsingArgs
	 * @param generalChromoInfo
	 * @param myOutDir
	 * @param refPanels
	 */
	private void initializeForMinimac(ParseCmdLine parsingArgs, String myOutDir, List<String> refPanels) {
		int chunkSize = parsingArgs.getChunkSize();

		// We create the first directory name: the cohort directory.
		String cohort = parsingArgs.getCohort();
		String tmpOutDir = myOutDir + File.separator + cohort;

		// Now we create the directories for imputedMMOutDir
		for (int j = 0; j < refPanels.size(); j++) {
			String rPanel = refPanels.get(j);
			String tmpPanelDir = tmpOutDir + File.separator + rPanel;
			// Next level: Create directories.
			String mixOutDir = tmpPanelDir + File.separator + "mixed";

			ArrayList<String> chromoListImputedMMOutDir = new ArrayList<String>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMDoseVCFFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedFileBgzip = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedFileTbi = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMInfoFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMErateFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMRecFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFFileBgzip = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFFileTbi = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMLogFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredLogFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdFile = new ArrayList<ArrayList<GenericFile>>();

			ArrayList<ArrayList<GenericFile>> chromoListImputedMMDoseVCFMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMalesFileBgzip = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMalesFileTbi = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMInfoMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMErateMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMRecMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFMalesFileBgzip = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFMalesFileTbi = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMLogMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredLogMalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdMalesFile = new ArrayList<ArrayList<GenericFile>>();

			ArrayList<ArrayList<GenericFile>> chromoListImputedMMDoseVCFFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedFemalesFileBgzip = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedFemalesFileTbi = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMInfoFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMErateFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMRecFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFFemalesFileBgzip = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMM3VCFFemalesFileTbi = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListImputedMMLogFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredLogFemalesFile = new ArrayList<ArrayList<GenericFile>>();
			ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdFemalesFile = new ArrayList<ArrayList<GenericFile>>();

			for (int i = this.startChr; i <= this.endChr; i++) {
				int chromo = i;
				int lim1 = 1;
				int lim2 = lim1 + chunkSize - 1;
				String tmpChrDir = mixOutDir + File.separator + "Chr_" + chromo;
				chromoListImputedMMOutDir.add(tmpChrDir);

				ArrayList<GenericFile> chunkListImputedMMDoseVCFFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedFileBgzip = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedFileTbi = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMInfoFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMErateFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMRecFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFFileBgzip = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFFileTbi = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMLogFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredLogFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredRsIdFile = new ArrayList<GenericFile>();

				ArrayList<GenericFile> chunkListImputedMMDoseVCFMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMalesFileBgzip = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMalesFileTbi = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMInfoMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMErateMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMRecMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFMalesFileBgzip = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFMalesFileTbi = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMLogMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredLogMalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredRsIdMalesFile = new ArrayList<GenericFile>();

				ArrayList<GenericFile> chunkListImputedMMDoseVCFFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedFemalesFileBgzip = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedFemalesFileTbi = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMInfoFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMErateFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMRecFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFFemalesFileBgzip = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMM3VCFFemalesFileTbi = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListImputedMMLogFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredLogFemalesFile = new ArrayList<GenericFile>();
				ArrayList<GenericFile> chunkListFilteredRsIdFemalesFile = new ArrayList<GenericFile>();

				int numberOfChunks = ChromoInfo.getMaxSize(chromo) / chunkSize;
				int module = ChromoInfo.getMaxSize(chromo) % chunkSize;
				if (module != 0)
					numberOfChunks++;

				for (int k = 0; k < numberOfChunks; k++) {

					if (chromo < 23) {

						String imputedMMFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_minimac";
						GenericFile myChunkListImputedMMDoseVCFFile = new GenericFile(tmpChrDir,
								imputedMMFileName + ".dose.vcf.gz", "decompressed", "none");
						chunkListImputedMMDoseVCFFile.add(myChunkListImputedMMDoseVCFFile);

						GenericFile myChunkListImputedFileBgzip = new GenericFile(tmpChrDir,
								imputedMMFileName + ".dose.vcf.gz", "decompressed", "none");
						chunkListImputedFileBgzip.add(myChunkListImputedFileBgzip);

						GenericFile myChunkListImputedFileTbi = new GenericFile(tmpChrDir,
								imputedMMFileName + ".dose.vcf.gz.tbi", "decompressed", "none");
						chunkListImputedFileTbi.add(myChunkListImputedFileTbi);

						GenericFile myChunkListImputedMMInfoFile = new GenericFile(tmpChrDir,
								imputedMMFileName + ".info", "decompressed", "none");
						chunkListImputedMMInfoFile.add(myChunkListImputedMMInfoFile);

						GenericFile myChunkListImputedMMErateFile = new GenericFile(tmpChrDir,
								imputedMMFileName + ".erate", "decompressed", "none");
						chunkListImputedMMErateFile.add(myChunkListImputedMMErateFile);

						GenericFile myChunkListImputedMMRecFile = new GenericFile(tmpChrDir, imputedMMFileName + ".rec",
								"decompressed", "none");
						chunkListImputedMMRecFile.add(myChunkListImputedMMRecFile);

						GenericFile myChunkListImputedMMM3VCFFile = new GenericFile(tmpChrDir,
								imputedMMFileName + ".m3vcf.gz", "decompressed", "none");
						chunkListImputedMMM3VCFFile.add(myChunkListImputedMMM3VCFFile);

						GenericFile myChunkListImputedMMM3VCFFileBgzip = new GenericFile(tmpChrDir,
								imputedMMFileName + ".m3vcf.gz", "decompressed", "none");
						chunkListImputedMMM3VCFFileBgzip.add(myChunkListImputedMMM3VCFFileBgzip);

						GenericFile myChunkListImputedMMM3VCFFileTbi = new GenericFile(tmpChrDir,
								imputedMMFileName + ".m3vcf.gz.tbi", "decompressed", "none");
						chunkListImputedMMM3VCFFileTbi.add(myChunkListImputedMMM3VCFFileTbi);

						GenericFile myChunkListImputedMMLogFile = new GenericFile(tmpChrDir,
								imputedMMFileName + ".logfile", "decompressed", "none");
						chunkListImputedMMLogFile.add(myChunkListImputedMMLogFile);

						String filteredFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_filtered.minimac.gz";
						GenericFile myChunkListFilteredFile = new GenericFile(tmpChrDir, filteredFileName, "compressed",
								"none");
						chunkListFilteredFile.add(myChunkListFilteredFile);

						String filteredLogFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_filtered.minimac.log";
						GenericFile myChunkListFilteredLogFile = new GenericFile(tmpChrDir, filteredLogFileName,
								"decompressed", "none");
						chunkListFilteredLogFile.add(myChunkListFilteredLogFile);

						// //String filteredRsIdFileName = parsingArgs.getGenFileName(chromo) + "_" +
						// rPanel + "_" + lim1 + "_" + lim2 +"_filtered_rsid.txt";
						String filteredRsIdFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_filtered_rsid.txt";
						GenericFile myChunkListFilteredRsIdFile = new GenericFile(tmpChrDir, filteredRsIdFileName,
								"compressed", "none");
						chunkListFilteredRsIdFile.add(myChunkListFilteredRsIdFile);
					} else if (chromo == 23) {

						String imputedMMMalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_males_minimac";
						GenericFile myChunkListImputedMMDoseVCFMalesFile = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".dose.vcf.gz", "decompressed", "none");
						chunkListImputedMMDoseVCFMalesFile.add(myChunkListImputedMMDoseVCFMalesFile);

						GenericFile myChunkListImputedMalesFileBgzip = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".dose.vcf.gz", "decompressed", "none");
						chunkListImputedMalesFileBgzip.add(myChunkListImputedMalesFileBgzip);

						GenericFile myChunkListImputedMalesFileTbi = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".dose.vcf.gz.tbi", "decompressed", "none");
						chunkListImputedMalesFileTbi.add(myChunkListImputedMalesFileTbi);

						GenericFile myChunkListImputedMMInfoMalesFile = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".info", "decompressed", "none");
						chunkListImputedMMInfoMalesFile.add(myChunkListImputedMMInfoMalesFile);

						GenericFile myChunkListImputedMMErateMalesFile = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".erate", "decompressed", "none");
						chunkListImputedMMErateMalesFile.add(myChunkListImputedMMErateMalesFile);

						GenericFile myChunkListImputedMMRecMalesFile = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".rec", "decompressed", "none");
						chunkListImputedMMRecMalesFile.add(myChunkListImputedMMRecMalesFile);

						GenericFile myChunkListImputedMMM3VCFMalesFile = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".m3vcf.gz", "decompressed", "none");
						chunkListImputedMMM3VCFMalesFile.add(myChunkListImputedMMM3VCFMalesFile);

						GenericFile myChunkListImputedMMM3VCFMalesFileBgzip = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".m3vcf.gz", "decompressed", "none");
						chunkListImputedMMM3VCFMalesFileBgzip.add(myChunkListImputedMMM3VCFMalesFileBgzip);

						GenericFile myChunkListImputedMMM3VCFMalesFileTbi = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".m3vcf.gz.tbi", "decompressed", "none");
						chunkListImputedMMM3VCFMalesFileTbi.add(myChunkListImputedMMM3VCFMalesFileTbi);

						GenericFile myChunkListImputedMMLogMalesFile = new GenericFile(tmpChrDir,
								imputedMMMalesFileName + ".logfile", "decompressed", "none");
						chunkListImputedMMLogMalesFile.add(myChunkListImputedMMLogMalesFile);

						String filteredMalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_males_filtered.minimac.gz";
						GenericFile myChunkListFilteredMalesFile = new GenericFile(tmpChrDir, filteredMalesFileName,
								"compressed", "none");
						chunkListFilteredMalesFile.add(myChunkListFilteredMalesFile);

						String filteredLogMalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_males_filtered.minimac.log";
						GenericFile myChunkListFilteredLogMalesFile = new GenericFile(tmpChrDir,
								filteredLogMalesFileName, "decompressed", "none");
						chunkListFilteredLogMalesFile.add(myChunkListFilteredLogMalesFile);
						String filteredRsIdMalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_"
								+ lim2 + "_males_filtered_rsid.txt";
						GenericFile myChunkListFilteredRsIdMalesFile = new GenericFile(tmpChrDir,
								filteredRsIdMalesFileName, "compressed", "none");
						chunkListFilteredRsIdMalesFile.add(myChunkListFilteredRsIdMalesFile);

						String imputedMMFemalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_females_minimac";
						GenericFile myChunkListImputedMMDoseVCFFemalesFile = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".dose.vcf.gz", "decompressed", "none");
						chunkListImputedMMDoseVCFFemalesFile.add(myChunkListImputedMMDoseVCFFemalesFile);

						GenericFile myChunkListImputedFemalesFileBgzip = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".dose.vcf.gz", "decompressed", "none");
						chunkListImputedFemalesFileBgzip.add(myChunkListImputedFemalesFileBgzip);

						GenericFile myChunkListImputedFemalesFileTbi = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".dose.vcf.gz.tbi", "decompressed", "none");
						chunkListImputedFemalesFileTbi.add(myChunkListImputedFemalesFileTbi);

						GenericFile myChunkListImputedMMInfoFemalesFile = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".info", "decompressed", "none");
						chunkListImputedMMInfoFemalesFile.add(myChunkListImputedMMInfoFemalesFile);

						GenericFile myChunkListImputedMMErateFemalesFile = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".erate", "decompressed", "none");
						chunkListImputedMMErateFemalesFile.add(myChunkListImputedMMErateFemalesFile);

						GenericFile myChunkListImputedMMRecFemalesFile = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".rec", "decompressed", "none");
						chunkListImputedMMRecFemalesFile.add(myChunkListImputedMMRecFemalesFile);

						GenericFile myChunkListImputedMMM3VCFFemalesFile = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".m3vcf.gz", "decompressed", "none");
						chunkListImputedMMM3VCFFemalesFile.add(myChunkListImputedMMM3VCFFemalesFile);

						GenericFile myChunkListImputedMMM3VCFFemalesFileBgzip = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".m3vcf.gz", "decompressed", "none");
						chunkListImputedMMM3VCFFemalesFileBgzip.add(myChunkListImputedMMM3VCFFemalesFileBgzip);

						GenericFile myChunkListImputedMMM3VCFFemalesFileTbi = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".m3vcf.gz.tbi", "decompressed", "none");
						chunkListImputedMMM3VCFFemalesFileTbi.add(myChunkListImputedMMM3VCFFemalesFileTbi);

						GenericFile myChunkListImputedMMLogFemalesFile = new GenericFile(tmpChrDir,
								imputedMMFemalesFileName + ".logfile", "decompressed", "none");
						chunkListImputedMMLogFemalesFile.add(myChunkListImputedMMLogFemalesFile);

						String filteredFemalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
								+ "_females_filtered.minimac.gz";
						GenericFile myChunkListFilteredFemalesFile = new GenericFile(tmpChrDir, filteredFemalesFileName,
								"compressed", "none");
						chunkListFilteredFemalesFile.add(myChunkListFilteredFemalesFile);

						String filteredLogFemalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_"
								+ lim2 + "_females_filtered.minimac.log";
						GenericFile myChunkListFilteredLogFemalesFile = new GenericFile(tmpChrDir,
								filteredLogFemalesFileName, "decompressed", "none");
						chunkListFilteredLogFemalesFile.add(myChunkListFilteredLogFemalesFile);

						String filteredRsIdFemalesFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_"
								+ lim2 + "_females_filtered_rsid.txt";
						GenericFile myChunkListFilteredRsIdFemalesFile = new GenericFile(tmpChrDir,
								filteredRsIdFemalesFileName, "compressed", "none");
						chunkListFilteredRsIdFemalesFile.add(myChunkListFilteredRsIdFemalesFile);
					}

					lim1 = lim1 + chunkSize;
					lim2 = lim2 + chunkSize;
				}
				if (chromo < 23) {
					chromoListImputedMMDoseVCFFile.add(chunkListImputedMMDoseVCFFile);
					chromoListImputedFileBgzip.add(chunkListImputedFileBgzip);
					chromoListImputedFileTbi.add(chunkListImputedFileTbi);
					chromoListImputedMMInfoFile.add(chunkListImputedMMInfoFile);
					chromoListImputedMMErateFile.add(chunkListImputedMMErateFile);
					chromoListImputedMMRecFile.add(chunkListImputedMMRecFile);
					chromoListImputedMMM3VCFFile.add(chunkListImputedMMM3VCFFile);
					chromoListImputedMMM3VCFFileBgzip.add(chunkListImputedMMM3VCFFileBgzip);
					chromoListImputedMMM3VCFFileTbi.add(chunkListImputedMMM3VCFFileTbi);
					chromoListImputedMMLogFile.add(chunkListImputedMMLogFile);
					chromoListFilteredFile.add(chunkListFilteredFile);
					chromoListFilteredLogFile.add(chunkListFilteredLogFile);
					chromoListFilteredRsIdFile.add(chunkListFilteredRsIdFile);
				} else if (chromo == 23) {

					chromoListImputedMMDoseVCFMalesFile.add(chunkListImputedMMDoseVCFMalesFile);
					chromoListImputedMalesFileBgzip.add(chunkListImputedMalesFileBgzip);
					chromoListImputedMalesFileTbi.add(chunkListImputedMalesFileTbi);
					chromoListImputedMMInfoMalesFile.add(chunkListImputedMMInfoMalesFile);
					chromoListImputedMMErateMalesFile.add(chunkListImputedMMErateMalesFile);
					chromoListImputedMMRecMalesFile.add(chunkListImputedMMRecMalesFile);
					chromoListImputedMMM3VCFMalesFile.add(chunkListImputedMMM3VCFMalesFile);
					chromoListImputedMMM3VCFMalesFileBgzip.add(chunkListImputedMMM3VCFMalesFileBgzip);
					chromoListImputedMMM3VCFMalesFileTbi.add(chunkListImputedMMM3VCFMalesFileTbi);
					chromoListImputedMMLogMalesFile.add(chunkListImputedMMLogMalesFile);
					chromoListFilteredMalesFile.add(chunkListFilteredMalesFile);
					chromoListFilteredLogMalesFile.add(chunkListFilteredLogMalesFile);
					chromoListFilteredRsIdMalesFile.add(chunkListFilteredRsIdMalesFile);

					chromoListImputedMMDoseVCFFemalesFile.add(chunkListImputedMMDoseVCFFemalesFile);
					chromoListImputedFemalesFileBgzip.add(chunkListImputedFemalesFileBgzip);
					chromoListImputedFemalesFileTbi.add(chunkListImputedFemalesFileTbi);
					chromoListImputedMMInfoFemalesFile.add(chunkListImputedMMInfoFemalesFile);
					chromoListImputedMMErateFemalesFile.add(chunkListImputedMMErateFemalesFile);
					chromoListImputedMMRecFemalesFile.add(chunkListImputedMMRecFemalesFile);
					chromoListImputedMMM3VCFFemalesFile.add(chunkListImputedMMM3VCFFemalesFile);
					chromoListImputedMMM3VCFFemalesFileBgzip.add(chunkListImputedMMM3VCFFemalesFileBgzip);
					chromoListImputedMMM3VCFFemalesFileTbi.add(chunkListImputedMMM3VCFFemalesFileTbi);
					chromoListImputedMMLogFemalesFile.add(chunkListImputedMMLogFemalesFile);
					chromoListFilteredFemalesFile.add(chunkListFilteredFemalesFile);
					chromoListFilteredLogFemalesFile.add(chunkListFilteredLogFemalesFile);
					chromoListFilteredRsIdFemalesFile.add(chunkListFilteredRsIdFemalesFile);
				}

				if (chromo < 23) {

					imputedMMOutDir.add(chromoListImputedMMOutDir);
					imputedMMDoseVCFFile.add(chromoListImputedMMDoseVCFFile);
					imputedFileBgzip.add(chromoListImputedFileBgzip);
					imputedFileTbi.add(chromoListImputedFileTbi);
					imputedMMInfoFile.add(chromoListImputedMMInfoFile);
					imputedMMErateFile.add(chromoListImputedMMErateFile);
					imputedMMRecFile.add(chromoListImputedMMRecFile);
					imputedMMM3VCFFile.add(chromoListImputedMMM3VCFFile);
					imputedMMM3VCFFileBgzip.add(chromoListImputedMMM3VCFFileBgzip);
					imputedMMM3VCFFileTbi.add(chromoListImputedMMM3VCFFileTbi);
					imputedMMLogFile.add(chromoListImputedMMLogFile);
					filteredFile.add(chromoListFilteredFile);
					filteredFileLogFile.add(chromoListFilteredLogFile);
					filteredFileRsIdFile.add(chromoListFilteredRsIdFile);

					System.out.println("mon :: INIT i :: " + i + "\n");
					System.out.println("mon :: i :: " + i + " size :: " + imputedMMInfoFile.size());
					System.out.println("mon :: startChr :: " + startChr + " endChr :: " + endChr);
					System.out.println("");

				} else if (chromo == 23) {

					imputedMMDoseVCFMalesFile.add(chromoListImputedMMDoseVCFMalesFile);
					imputedMalesFileBgzip.add(chromoListImputedMalesFileBgzip);
					imputedMalesFileTbi.add(chromoListImputedMalesFileTbi);
					imputedMMInfoMalesFile.add(chromoListImputedMMInfoMalesFile);
					imputedMMErateMalesFile.add(chromoListImputedMMErateMalesFile);
					imputedMMRecMalesFile.add(chromoListImputedMMRecMalesFile);
					imputedMMM3VCFMalesFile.add(chromoListImputedMMM3VCFMalesFile);
					imputedMMM3VCFMalesFileBgzip.add(chromoListImputedMMM3VCFMalesFileBgzip);
					imputedMMM3VCFMalesFileTbi.add(chromoListImputedMMM3VCFMalesFileTbi);
					imputedMMLogMalesFile.add(chromoListImputedMMLogMalesFile);
					filteredMalesFile.add(chromoListFilteredMalesFile);
					filteredLogMalesFile.add(chromoListFilteredLogMalesFile);
					filteredRsIdMalesFile.add(chromoListFilteredRsIdMalesFile);

					imputedMMDoseVCFFemalesFile.add(chromoListImputedMMDoseVCFFemalesFile);
					imputedFemalesFileBgzip.add(chromoListImputedFemalesFileBgzip);
					imputedFemalesFileTbi.add(chromoListImputedFemalesFileTbi);
					imputedMMInfoFemalesFile.add(chromoListImputedMMInfoFemalesFile);
					imputedMMErateFemalesFile.add(chromoListImputedMMErateFemalesFile);
					imputedMMRecFemalesFile.add(chromoListImputedMMRecFemalesFile);
					imputedMMM3VCFFemalesFile.add(chromoListImputedMMM3VCFFemalesFile);
					imputedMMM3VCFFemalesFileBgzip.add(chromoListImputedMMM3VCFFemalesFileBgzip);
					imputedMMM3VCFFemalesFileTbi.add(chromoListImputedMMM3VCFFemalesFileTbi);
					imputedMMLogFemalesFile.add(chromoListImputedMMLogFemalesFile);
					filteredFemalesFile.add(chromoListFilteredFemalesFile);
					filteredLogFemalesFile.add(chromoListFilteredLogFemalesFile);
					filteredRsIdFemalesFile.add(chromoListFilteredRsIdFemalesFile);
				}
			}
		}
	}

	/**
	 * Method to access the outputDir for the imputation process
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getOutputDir(int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int indexChr = chromo - this.startChr;
		return this.imputedOutDir.get(rPanelIndex).get(indexChr);
	}

	/**
	 * Method to access the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to get the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to get the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to get the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.imputedFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.imputedInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMalesInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedInfoMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFemalesInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedInfoFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to get the finalStatus of the imputeInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to get the finalStatus of the imputeInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedInfoMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedInfoMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to get the finalStatus of the imputeInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedInfoFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedInfoFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to access the imputedInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedInfoFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedInfoMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedInfoMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedInfoFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedInfoFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedInfoFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedSummaryFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedSummaryFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedSummaryFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedSummaryFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedSummaryMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedSummaryMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedSummaryFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedSummaryFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedSummaryFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedWarningsFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedWarningsFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedWarningsFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedWarningsFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedWarningsMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedWarningsMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedWarningsFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedWarningsFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedWarningsFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedLogFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedLogMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedLogMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedLogFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedLogFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMDoseVCFFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedMMDoseVCFFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMDoseVCFMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMDoseVCFMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMDoseVCFFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMDoseVCFFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMDoseVCFFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.imputedMMDoseVCFFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMDoseVCFMalesFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMDoseVCFMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMDoseVCFFemalesFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMDoseVCFFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFileBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return imputedFileBgzip.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMalesFileBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return imputedMalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFemalesFileBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return imputedFemalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFileBgzipFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		imputedFileBgzip.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMalesFileBgzipFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFemalesFileBgzipFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedFemalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFileTbi(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return imputedFileTbi.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMalesFileTbi(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return imputedMalesFileTbi.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedFemalesFileTbi(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return imputedFemalesFileTbi.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFileTbiFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		imputedFileTbi.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMalesFileTbiFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMalesFileTbi.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedFemalesFileTbiFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedFemalesFileTbi.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the imputedMMFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMInfoFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedMMInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMInfoMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMInfoMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMInfoFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMInfoFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - startChr;
		int indexChunk = lim1 / chunkSize;
		imputedMMInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMInfoMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMMInfoMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMInfoFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMMInfoFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the imputedMMErateFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMErateFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedMMErateFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMErateFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMErateMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMErateMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMErateFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMErateFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMErateFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to set the finalStatus of the imputedMMErateFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMErateFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.imputedMMErateFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedMMErateFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMErateMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMErateMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedMMErateFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMErateFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMErateFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the imputedMMRecFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMRecFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedMMRecFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMRecFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMRecMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMRecMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMRecFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMRecFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMRecFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to set the finalStatus of the imputeMMRecFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMRecFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.imputedMMRecFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputeMMRecFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMRecMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMRecMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputeMMRecFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMRecFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMRecFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the imputedMMM3VCFFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedMMM3VCFFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMM3VCFFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMM3VCFMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the imputedMMM3VCFFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMM3VCFFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to set the finalStatus of the imputedMMDoseFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - startChr;
		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedMMDoseFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set the finalStatus of the imputedMMDoseFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFFileBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - startChr;
		int indexChunk = lim1 / chunkSize;
		return imputedMMM3VCFFileBgzip.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFMalesFileBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return imputedMMM3VCFMalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFFemalesFileBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return imputedMMM3VCFFemalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFFileFinalStatusBgzip(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - startChr;
		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFFileBgzip.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFMalesFileFinalStatusBgzip(int rPanelIndex, int chromo, int lim1, int lim2,
			int chunkSize, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFMalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFFemalesFileFinalStatusBgzip(int rPanelIndex, int chromo, int lim1, int lim2,
			int chunkSize, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFFemalesFileBgzip.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMM3VCFFileTbi(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - startChr;
		int indexChunk = lim1 / chunkSize;
		return imputedMMM3VCFFileTbi.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}

	/**
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMM3VCFFileFinalStatusTbi(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - startChr;
		int indexChunk = lim1 / chunkSize;
		imputedMMM3VCFFileTbi.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the imputedMMLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMLogFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.imputedMMLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}
	
	/**
	 * Method to access the imputedMMLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMLogMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMLogMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}
	
	/**
	 * Method to access the imputedMMLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getImputedMMLogFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.imputedMMLogFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to set the finalStatus of the imputedMMLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMLogFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.imputedMMLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}
	
	/**
	 * Method to set the finalStatus of the imputedMMLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMLogMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMLogMalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}
	
	/**
	 * Method to set the finalStatus of the imputedMMLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setImputedMMLogFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.imputedMMLogFemalesFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.filteredFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName() + ".gz";
	}
	
	/**
	 * Method to access the filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName() + ".gz";
	}
	
	/**
	 * Method to access the filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName() + ".gz";
	}

	/**
	 * Method to set the finalStatus of filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setFilteredFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		this.filteredFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
	}
	
	/**
	 * Method to set the finalStatus of filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setFilteredMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.filteredFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}
	
	/**
	 * Method to set the finalStatus of filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @param finalStatus
	 */
	public void setFilteredFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		this.filteredFile.get(rPanelIndex).get(0).get(indexChunk).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.filteredFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFinalStatus();
	}
	
	/**
	 * Method to access the finalStatus of filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredMalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFinalStatus();
	}
	
	/**
	 * Method to access the finalStatus of filteredFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredFemalesFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFinalStatus();
	}

	/**
	 * Method to access the filteredLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredLogFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.filteredFileLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}
	
	/**
	 * Method to access the filteredLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredLogMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredLogMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}
	
	/**
	 * Method to access the filteredLogFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredLogFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredLogFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to access the filteredRsIdFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredRsIdFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexChunk = lim1 / chunkSize;
		return this.filteredFileRsIdFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
	}
	
	/**
	 * Method to access the filteredRsIdFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredRsIdMalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredRsIdMalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}
	
	/**
	 * Method to access the filteredRsIdFile
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getFilteredRsIdFemalesFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChunk = lim1 / chunkSize;
		return this.filteredRsIdFemalesFile.get(rPanelIndex).get(0).get(indexChunk).getFullName();
	}

	/**
	 * Method to print all the impute files that will be created
	 * 
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 */
	public void printImputationFiles(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int indexChr = chromo - this.startChr;
		int indexLow = lim1 / chunkSize;
		int indexHigh = lim2 / chunkSize;
		LOGGER.info("-------------------------------------------------");
		LOGGER.info("Files information for the chromosome " + chromo);
		LOGGER.info(
				"OutDir[" + rPanelIndex + "][Chr_" + chromo + "]=" + this.imputedOutDir.get(rPanelIndex).get(indexChr));

		for (int j = indexLow; j < indexHigh; j++) {
			LOGGER.info("        ImputedFile[" + rPanelIndex + "][Chr_" + chromo + "][" + lim1 + "-" + lim2 + "]="
					+ imputedFile.get(rPanelIndex).get(indexChr).get(indexLow).getFullName());
			LOGGER.info("    ImputedInfoFile[" + rPanelIndex + "][Chr_" + chromo + "][" + lim1 + "-" + lim2 + "]="
					+ imputedInfoFile.get(rPanelIndex).get(indexChr).get(indexLow).getFullName());
			LOGGER.info(" ImputedSummaryFile[" + rPanelIndex + "][Chr_" + chromo + "][" + lim1 + "-" + lim2 + "]="
					+ imputedSummaryFile.get(rPanelIndex).get(indexChr).get(indexLow).getFullName());
			LOGGER.info("ImputedWarningsFile[" + rPanelIndex + "][Chr_" + chromo + "][" + lim1 + "-" + lim2 + "]="
					+ imputedWarningsFile.get(rPanelIndex).get(indexChr).get(indexLow).getFullName());
			LOGGER.info("     ImputedLogFile[" + rPanelIndex + "][Chr_" + chromo + "][" + lim1 + "-" + lim2 + "]="
					+ imputedLogFile.get(rPanelIndex).get(indexChr).get(indexLow).getFullName());
		}

		LOGGER.info("-------------------------------------------------");
	}

	private void checkChromoIndex(int chromo) {
		if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
			LOGGER.fatal("[ImputationFiles] Error, chromosomesome " + chromo + "does not exist");
			System.exit(1);
		}
	}

	private void checkLimits(int chromo, int lim1, int lim2) {
		if ((lim1 < MIN_LIMIT) || (lim2 > MAX_LIMIT)) {
			LOGGER.fatal(
					"[ImputationFiles] Error, Chunk " + lim1 + "_" + lim2 + "does not exist for chromosome " + chromo);
			System.exit(1);
		}
	}

}
