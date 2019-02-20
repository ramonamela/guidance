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

public class AssocFiles {

	// Logger
	private static final Logger LOGGER = LogManager.getLogger("Console");

	// Constants
	private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;
	private static final int MIN_LIMIT = 1;
	private static final int MAX_LIMIT = 252_000_000;
	private static final String UNCOMPRESSED = "uncompressed";

	private ArrayList<ArrayList<ArrayList<String>>> outDir = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestOutFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestOutMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestOutFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestLogFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestLogMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestLogFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFilteredFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFilteredMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFilteredFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryCondensedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryCondensedMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryCondensedFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedFilteredFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedFilteredMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedFilteredFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedCondensedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedCondensedMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedCondensedFemalesFile = new ArrayList<>();

	private int startChr = 0;
	private int endChr = 0;

	/**
	 * New AssocFiles instance
	 * 
	 * @param parsingArgs
	 * @param generalChromoInfo
	 * @param baseOutDir
	 * @param refPanels
	 */
	public AssocFiles(ParseCmdLine parsingArgs, String baseOutDir, List<String> refPanels) {
		this.startChr = parsingArgs.getStart();
		this.endChr = parsingArgs.getEnd();

		int chunkSize = parsingArgs.getChunkSize();
		int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();

		// We create the first directory name: the cohort directory
		String mixedCohort = parsingArgs.getCohort();

		for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
			String testTypeName = parsingArgs.getTestTypeName(tt);

			String testTypeOutDir = baseOutDir + File.separator + "associations" + File.separator + testTypeName;

			ArrayList<ArrayList<String>> rpanelListOutDir = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestOutFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestLogFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestOutMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestLogMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryMalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestOutFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestLogFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFemalesFile = new ArrayList<ArrayList<ArrayList<GenericFile>>>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFilteredFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFilteredMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFilteredFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryCondensedMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryCondensedFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedFilteredFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedFilteredMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedFilteredFemalesFile = new ArrayList<>();
			
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedCondensedMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedCondensedFemalesFile = new ArrayList<>();

			String rPanel = null;
			String testTypeOutDir2 = baseOutDir + File.separator + "associations" + File.separator + testTypeName
					+ File.separator + mixedCohort + "_combined_panels";
			for (int j = 0; j < refPanels.size(); j++) {
				rPanel = refPanels.get(j);
				testTypeOutDir2 = testTypeOutDir2 + "_" + rPanel;
			}

			rPanel = refPanels.get(0);
			String prefixFilteredName = "filteredByAll_results_" + testTypeName + "_" + mixedCohort + "_" + rPanel;
			String prefixCondensedName = "condensed_results_" + testTypeName + "_" + mixedCohort + "_" + rPanel;

			for (int j = 1; j < refPanels.size(); j++) {
				rPanel = refPanels.get(j);

				prefixFilteredName = prefixFilteredName + "_" + rPanel;
				prefixCondensedName = prefixCondensedName + "_" + rPanel;

				ArrayList<ArrayList<GenericFile>> chromoListCombinedFilteredFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedFilteredMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedFilteredFemalesFile = new ArrayList<>();
				
				ArrayList<ArrayList<GenericFile>> chromoListCombinedCondensedFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedCondensedMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedCondensedFemalesFile = new ArrayList<>();

				for (int chromo = this.startChr; chromo <= this.endChr; chromo++) {
					int maxSize = ChromoInfo.getMaxSize(chromo);
					int totalChunks = maxSize / chunkSize;
					int module = maxSize % chunkSize;
					if (module != 0) {
						totalChunks++;
					}
					int lim1 = 1;
					int lim2 = lim1 + chunkSize - 1;

					ArrayList<GenericFile> chunkListCombinedFilteredFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListCombinedFilteredMalesFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListCombinedFilteredFemalesFile = new ArrayList<>();

					ArrayList<GenericFile> chunkListCombinedCondensedFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListCombinedCondensedMalesFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListCombinedCondensedFemalesFile = new ArrayList<>();

					for (int k = 0; k < totalChunks; k++) {
						if (chromo == 23) {
							String tmpCombinedFilteredMalesFileName = prefixFilteredName + "_chr_" + chromo + "_" + lim1
									+ "_" + lim2 + "_combined_males.txt.gz";
							GenericFile myChunkListCombinedFilteredMalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedFilteredMalesFileName, UNCOMPRESSED, "none");
							chunkListCombinedFilteredMalesFile.add(myChunkListCombinedFilteredMalesFile);

							String tmpCombinedCondensedMalesFileName = prefixCondensedName + "_chr_" + chromo + "_"
									+ lim1 + "_" + lim2 + "_combined_males.txt.gz";
							GenericFile myChunkListCombinedCondensedMalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedCondensedMalesFileName, UNCOMPRESSED, "none");
							chunkListCombinedCondensedMalesFile.add(myChunkListCombinedCondensedMalesFile);

							String tmpCombinedFilteredFemalesFileName = prefixFilteredName + "_chr_" + chromo + "_"
									+ lim1 + "_" + lim2 + "_combined_females.txt.gz";
							GenericFile myChunkListCombinedFilteredFemalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedFilteredFemalesFileName, UNCOMPRESSED, "none");
							chunkListCombinedFilteredFemalesFile.add(myChunkListCombinedFilteredFemalesFile);

							String tmpCombinedCondensedFemalesFileName = prefixCondensedName + "_chr_" + chromo + "_"
									+ lim1 + "_" + lim2 + "_combined_females.txt.gz";
							GenericFile myChunkListCombinedCondensedFemalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedCondensedFemalesFileName, UNCOMPRESSED, "none");
							chunkListCombinedCondensedFemalesFile.add(myChunkListCombinedCondensedFemalesFile);
						} else {
							String tmpCombinedFilteredFileName = prefixFilteredName + "_chr_" + chromo + "_" + lim1
									+ "_" + lim2 + "_combined.txt.gz";
							GenericFile myChunkListCombinedFilteredFile = new GenericFile(testTypeOutDir2,
									tmpCombinedFilteredFileName, UNCOMPRESSED, "none");
							chunkListCombinedFilteredFile.add(myChunkListCombinedFilteredFile);

							String tmpCombinedCondensedFileName = prefixCondensedName + "_chr_" + chromo + "_" + lim1
									+ "_" + lim2 + "_combined.txt.gz";
							GenericFile myChunkListCombinedCondensedFile = new GenericFile(testTypeOutDir2,
									tmpCombinedCondensedFileName, UNCOMPRESSED, "none");
							chunkListCombinedCondensedFile.add(myChunkListCombinedCondensedFile);
						}

						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}

					if (chromo == 23) {
						chromoListCombinedFilteredMalesFile.add(chunkListCombinedFilteredMalesFile);
						chromoListCombinedFilteredFemalesFile.add(chunkListCombinedFilteredFemalesFile);
						
						chromoListCombinedCondensedMalesFile.add(chunkListCombinedCondensedMalesFile);
						chromoListCombinedCondensedFemalesFile.add(chunkListCombinedCondensedFemalesFile);
					} else {
						chromoListCombinedFilteredFile.add(chunkListCombinedFilteredFile);
						chromoListCombinedCondensedFile.add(chunkListCombinedCondensedFile);
					}

					

				} // End for Chromo

				rpanelListCombinedFilteredFile.add(chromoListCombinedFilteredFile);
				rpanelListCombinedFilteredMalesFile.add(chromoListCombinedFilteredMalesFile);
				rpanelListCombinedFilteredFemalesFile.add(chromoListCombinedFilteredFemalesFile);

				rpanelListCombinedCondensedFile.add(chromoListCombinedCondensedFile);
				rpanelListCombinedCondensedMalesFile.add(chromoListCombinedCondensedMalesFile);
				rpanelListCombinedCondensedFemalesFile.add(chromoListCombinedCondensedFemalesFile);

			} // End of for refPanels

			this.combinedFilteredFile.add(rpanelListCombinedFilteredFile);
			this.combinedFilteredMalesFile.add(rpanelListCombinedFilteredMalesFile);
			this.combinedFilteredFemalesFile.add(rpanelListCombinedFilteredFemalesFile);

			this.combinedCondensedFile.add(rpanelListCombinedCondensedFile);
			this.combinedCondensedMalesFile.add(rpanelListCombinedCondensedMalesFile);
			this.combinedCondensedFemalesFile.add(rpanelListCombinedCondensedFemalesFile);

			for (int j = 0; j < refPanels.size(); j++) {
				// String rpanel = refPanels.get(j);
				rPanel = refPanels.get(j);
				String rpanelOutDir = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel;

				ArrayList<String> chromoListOutDir = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListSnptestOutFile = new ArrayList<ArrayList<GenericFile>>();
				ArrayList<ArrayList<GenericFile>> chromoListSnptestLogFile = new ArrayList<ArrayList<GenericFile>>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryFile = new ArrayList<ArrayList<GenericFile>>();

				ArrayList<ArrayList<GenericFile>> chromoListSnptestOutMalesFile = new ArrayList<ArrayList<GenericFile>>();
				ArrayList<ArrayList<GenericFile>> chromoListSnptestLogMalesFile = new ArrayList<ArrayList<GenericFile>>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryMalesFile = new ArrayList<ArrayList<GenericFile>>();

				ArrayList<ArrayList<GenericFile>> chromoListSnptestOutFemalesFile = new ArrayList<ArrayList<GenericFile>>();
				ArrayList<ArrayList<GenericFile>> chromoListSnptestLogFemalesFile = new ArrayList<ArrayList<GenericFile>>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryFemalesFile = new ArrayList<ArrayList<GenericFile>>();

				ArrayList<ArrayList<GenericFile>> chromoListSummaryFilteredFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryFilteredMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryFilteredFemalesFile = new ArrayList<>();

				ArrayList<ArrayList<GenericFile>> chromoListSummaryCondensedFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryCondensedMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListSummaryCondensedFemalesFile = new ArrayList<>();

				for (int i = this.startChr; i <= this.endChr; i++) {
					int chromo = i;

					String tmpChrDir = rpanelOutDir + File.separator + "Chr_" + chromo;
					chromoListOutDir.add(tmpChrDir);

					int maxSize = ChromoInfo.getMaxSize(chromo);
					int totalChunks = maxSize / chunkSize;
					int module = maxSize % chunkSize;
					if (module != 0)
						totalChunks++;
					
					int lim1 = 1;
					int lim2 = lim1 + chunkSize - 1;

					ArrayList<GenericFile> chunkListSnptestOutFile = new ArrayList<GenericFile>();
					ArrayList<GenericFile> chunkListSnptestLogFile = new ArrayList<GenericFile>();
					ArrayList<GenericFile> chunkListSummaryFile = new ArrayList<GenericFile>();

					ArrayList<GenericFile> chunkListSnptestOutMalesFile = new ArrayList<GenericFile>();
					ArrayList<GenericFile> chunkListSnptestLogMalesFile = new ArrayList<GenericFile>();
					ArrayList<GenericFile> chunkListSummaryMalesFile = new ArrayList<GenericFile>();

					ArrayList<GenericFile> chunkListSnptestOutFemalesFile = new ArrayList<GenericFile>();
					ArrayList<GenericFile> chunkListSnptestLogFemalesFile = new ArrayList<GenericFile>();
					ArrayList<GenericFile> chunkListSummaryFemalesFile = new ArrayList<GenericFile>();

					ArrayList<GenericFile> chunkListSummaryFilteredFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListSummaryFilteredMalesFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListSummaryFilteredFemalesFile = new ArrayList<>();

					ArrayList<GenericFile> chunkListSummaryCondensedFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListSummaryCondensedMalesFile = new ArrayList<>();
					ArrayList<GenericFile> chunkListSummaryCondensedFemalesFile = new ArrayList<>();

					for (int k = 0; k < totalChunks; k++) {
						
						if (chromo == 23) {
							// Now we have to create the impute files for this iteration
							String tmpSnptestOutMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_males_snptest.out.gz";

							GenericFile myChunkListSnptestOutMalesFile = new GenericFile(tmpChrDir,
									tmpSnptestOutMalesFileName, UNCOMPRESSED, "none");

							chunkListSnptestOutMalesFile.add(myChunkListSnptestOutMalesFile);

							String tmpSnptestLogMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_males_snptest.log";

							GenericFile myChunkListSnptestLogMalesFile = new GenericFile(tmpChrDir,
									tmpSnptestLogMalesFileName, UNCOMPRESSED, "none");

							chunkListSnptestLogMalesFile.add(myChunkListSnptestLogMalesFile);

							String tmpSummaryMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_"
									+ lim1 + "_" + lim2 + "_males_summary.txt.gz";

							GenericFile myChunkListSummaryMalesFile = new GenericFile(tmpChrDir,
									tmpSummaryMalesFileName, UNCOMPRESSED, "none");

							chunkListSummaryMalesFile.add(myChunkListSummaryMalesFile);

							// Now we have to create the impute files for this iteration
							String tmpSnptestOutFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_females_snptest.out.gz";

							GenericFile myChunkListSnptestOutFemalesFile = new GenericFile(tmpChrDir,
									tmpSnptestOutFemalesFileName, UNCOMPRESSED, "none");

							chunkListSnptestOutFemalesFile.add(myChunkListSnptestOutFemalesFile);

							String tmpSnptestLogFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_females_snptest.log";

							GenericFile myChunkListSnptestLogFemalesFile = new GenericFile(tmpChrDir,
									tmpSnptestLogFemalesFileName, UNCOMPRESSED, "none");

							chunkListSnptestLogFemalesFile.add(myChunkListSnptestLogFemalesFile);

							String tmpSummaryFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_"
									+ lim1 + "_" + lim2 + "_females_summary.txt.gz";

							GenericFile myChunkListSummaryFemalesFile = new GenericFile(tmpChrDir,
									tmpSummaryFemalesFileName, UNCOMPRESSED, "none");

							chunkListSummaryFemalesFile.add(myChunkListSummaryFemalesFile);

							String tmpSummaryFilteredMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_males_summary_filtered.txt.gz";
							GenericFile myChunkListSummaryFilteredMalesFile = new GenericFile(tmpChrDir,
									tmpSummaryFilteredMalesFileName, UNCOMPRESSED, "none");
							chunkListSummaryFilteredMalesFile.add(myChunkListSummaryFilteredMalesFile);

							String tmpSummaryFilteredFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_"
									+ rPanel + "_" + lim1 + "_" + lim2 + "_females_summary_filtered.txt.gz";
							GenericFile myChunkListSummaryFilteredFemalesFile = new GenericFile(tmpChrDir,
									tmpSummaryFilteredFemalesFileName, UNCOMPRESSED, "none");
							chunkListSummaryFilteredFemalesFile.add(myChunkListSummaryFilteredFemalesFile);

							String tmpSummaryCondensedMalesFileName = "chr_" + chromo + "_" + testTypeName + "_"
									+ rPanel + "_" + lim1 + "_" + lim2 + "_males_summary_condensed.txt.gz";
							GenericFile myChunkListSummaryCondensedMalesFile = new GenericFile(tmpChrDir,
									tmpSummaryCondensedMalesFileName, UNCOMPRESSED, "none");
							chunkListSummaryCondensedMalesFile.add(myChunkListSummaryCondensedMalesFile);

							String tmpSummaryCondensedFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_"
									+ rPanel + "_" + lim1 + "_" + lim2 + "_females_summary_condensed.txt.gz";
							GenericFile myChunkListSummaryCondensedFemalesFile = new GenericFile(tmpChrDir,
									tmpSummaryCondensedFemalesFileName, UNCOMPRESSED, "none");
							chunkListSummaryCondensedFemalesFile.add(myChunkListSummaryCondensedFemalesFile);

						} else {
							// Now we have to create the impute files for this iteration
							String tmpSnptestOutFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_"
									+ lim1 + "_" + lim2 + "_snptest.out.gz";
							GenericFile myChunkListSnptestOutFile = new GenericFile(tmpChrDir, tmpSnptestOutFileName,
									UNCOMPRESSED, "none");
							chunkListSnptestOutFile.add(myChunkListSnptestOutFile);

							String tmpSnptestLogFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_"
									+ lim1 + "_" + lim2 + "_snptest.log";
							GenericFile myChunkListSnptestLogFile = new GenericFile(tmpChrDir, tmpSnptestLogFileName,
									UNCOMPRESSED, "none");
							chunkListSnptestLogFile.add(myChunkListSnptestLogFile);

							String tmpSummaryFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_" + lim1
									+ "_" + lim2 + "_summary.txt.gz";
							GenericFile myChunkListSummaryFile = new GenericFile(tmpChrDir, tmpSummaryFileName,
									UNCOMPRESSED, "none");
							chunkListSummaryFile.add(myChunkListSummaryFile);

							String tmpSummaryFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_summary_filtered.txt.gz";
							GenericFile myChunkListSummaryFilteredFile = new GenericFile(tmpChrDir,
									tmpSummaryFilteredFileName, UNCOMPRESSED, "none");
							chunkListSummaryFilteredFile.add(myChunkListSummaryFilteredFile);

							String tmpSummaryCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ "_" + lim1 + "_" + lim2 + "_summary_condensed.txt.gz";
							GenericFile myChunkListSummaryCondensedFile = new GenericFile(tmpChrDir,
									tmpSummaryCondensedFileName, UNCOMPRESSED, "none");
							chunkListSummaryCondensedFile.add(myChunkListSummaryCondensedFile);
						}

						lim1 = lim1 + chunkSize;
						lim2 = lim2 + chunkSize;
					}

					if (chromo == 23) {
						chromoListSnptestOutMalesFile.add(chunkListSnptestOutMalesFile);
						chromoListSnptestLogMalesFile.add(chunkListSnptestLogMalesFile);
						chromoListSummaryMalesFile.add(chunkListSummaryMalesFile);
						chromoListSummaryCondensedMalesFile.add(chunkListSummaryCondensedMalesFile);
						chromoListSummaryFilteredMalesFile.add(chunkListSummaryFilteredMalesFile);

						chromoListSnptestOutFemalesFile.add(chunkListSnptestOutFemalesFile);
						chromoListSnptestLogFemalesFile.add(chunkListSnptestLogFemalesFile);
						chromoListSummaryFemalesFile.add(chunkListSummaryFemalesFile);
						chromoListSummaryCondensedFemalesFile.add(chunkListSummaryCondensedFemalesFile);
						chromoListSummaryFilteredFemalesFile.add(chunkListSummaryFilteredFemalesFile);
					} else {
						chromoListSnptestOutFile.add(chunkListSnptestOutFile);
						chromoListSnptestLogFile.add(chunkListSnptestLogFile);
						chromoListSummaryFile.add(chunkListSummaryFile);
						chromoListSummaryCondensedFile.add(chunkListSummaryCondensedFile);
						chromoListSummaryFilteredFile.add(chunkListSummaryFilteredFile);
					}

				} // End for chromo

				rpanelListOutDir.add(chromoListOutDir);

				rpanelListSnptestOutFile.add(chromoListSnptestOutFile);
				rpanelListSnptestLogFile.add(chromoListSnptestLogFile);
				rpanelListSummaryFile.add(chromoListSummaryFile);

				rpanelListSnptestOutMalesFile.add(chromoListSnptestOutMalesFile);
				rpanelListSnptestLogMalesFile.add(chromoListSnptestLogMalesFile);
				rpanelListSummaryMalesFile.add(chromoListSummaryMalesFile);

				rpanelListSnptestOutFemalesFile.add(chromoListSnptestOutFemalesFile);
				rpanelListSnptestLogFemalesFile.add(chromoListSnptestLogFemalesFile);
				rpanelListSummaryFemalesFile.add(chromoListSummaryFemalesFile);

				rpanelListSummaryFilteredFile.add(chromoListSummaryFilteredFile);
				rpanelListSummaryFilteredMalesFile.add(chromoListSummaryFilteredMalesFile);
				rpanelListSummaryFilteredFemalesFile.add(chromoListSummaryFilteredFemalesFile);

				rpanelListSummaryCondensedFile.add(chromoListSummaryCondensedFile);
				rpanelListSummaryCondensedMalesFile.add(chromoListSummaryCondensedMalesFile);
				rpanelListSummaryCondensedFemalesFile.add(chromoListSummaryCondensedFemalesFile);

			} // End of for refPanels

			this.snptestOutFile.add(rpanelListSnptestOutFile);
			this.snptestLogFile.add(rpanelListSnptestLogFile);
			this.summaryFile.add(rpanelListSummaryFile);

			this.snptestOutMalesFile.add(rpanelListSnptestOutMalesFile);
			this.snptestLogMalesFile.add(rpanelListSnptestLogMalesFile);
			this.summaryMalesFile.add(rpanelListSummaryMalesFile);

			this.snptestOutFemalesFile.add(rpanelListSnptestOutFemalesFile);
			this.snptestLogFemalesFile.add(rpanelListSnptestLogFemalesFile);
			this.summaryFemalesFile.add(rpanelListSummaryFemalesFile);

			this.summaryFilteredFile.add(rpanelListSummaryFilteredFile);
			this.summaryFilteredMalesFile.add(rpanelListSummaryFilteredMalesFile);
			this.summaryFilteredFemalesFile.add(rpanelListSummaryFilteredFemalesFile);

			this.summaryCondensedFile.add(rpanelListSummaryCondensedFile);
			this.summaryCondensedMalesFile.add(rpanelListSummaryCondensedMalesFile);
			this.summaryCondensedFemalesFile.add(rpanelListSummaryCondensedFemalesFile);

			this.outDir.add(rpanelListOutDir);
		} // End of for test types
	}

	/**
	 * Method to access outDir information
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getAssocOutDir(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;
		return this.outDir.get(testTypeIndex).get(rPanelIndex).get(i);
	}

	/**
	 * Method to access snptestOutFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSnptestOutFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access snptestOutFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSnptestOutMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		// The offset is because the array start in position 0
		int index = lim1 / chunkSize;
		return this.snptestOutMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access snptestOutFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSnptestOutFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		// The offset is because the array start in position 0
		int index = lim1 / chunkSize;
		return this.snptestOutFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/* Filtered -------------------------------------- */
	/**
	 * Method to access summaryFilteredFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryFilteredFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2,
			int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.summaryFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access summaryFilteredFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryFilteredMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		int index = lim1 / chunkSize;
		return this.summaryFilteredMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access summaryFilteredFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryFilteredFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		int index = lim1 / chunkSize;
		return this.summaryFilteredFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access combinedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getCombinedFilteredFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2,
			int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.combinedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access combinedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getCombinedFilteredMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		int index = lim1 / chunkSize;
		return this.combinedFilteredMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access combinedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getCombinedFilteredFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2,
			int chunkSize) {
		int index = lim1 / chunkSize;
		return this.combinedFilteredFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/* Condensed ------------------------------------- */

	/**
	 * Method to access summaryCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryCondensedFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2,
			int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.summaryCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access summaryCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryCondensedMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		int index = lim1 / chunkSize;
		return this.summaryCondensedMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access summaryCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryCondensedFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2,
			int chunkSize) {
		int index = lim1 / chunkSize;
		return this.summaryCondensedFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access combinedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getCombinedCondensedFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2,
			int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.combinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}
	
	/**
	 * Method to access combinedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getCombinedCondensedMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2,
			int chunkSize) {
		int index = lim1 / chunkSize;
		return this.combinedCondensedMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}
	
	/**
	 * Method to access combinedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getCombinedCondensedFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2,
			int chunkSize) {
		int index = lim1 / chunkSize;
		return this.combinedCondensedFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access snptestLogFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSnptestLogFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.snptestLogFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access snptestLogFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSnptestLogMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		int index = lim1 / chunkSize;
		return this.snptestLogMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access snptestLogFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSnptestLogFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2, int chunkSize) {
		int index = lim1 / chunkSize;
		return this.snptestLogFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access snptestLogFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		// The offset is because the array start in position 0
		int i = chromo - this.startChr;
		int index = lim1 / chunkSize;
		return this.summaryFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access snptestLogFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryMalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2,
			int chunkSize) {
		int index = lim1 / chunkSize;
		return this.summaryMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to access snptestLogFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 * @return
	 */
	public String getSummaryFemalesFile(int testTypeIndex, int rPanelIndex, int lim1, int lim2,
			int chunkSize) {
		int index = lim1 / chunkSize;
		return this.summaryFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index).getFullName();
	}

	/**
	 * Method to print Association files information
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param lim1
	 * @param lim2
	 * @param chunkSize
	 */
	public void printAssocFiles(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		// Check limits are within the bounds
		checkLimits(chromo, lim1, lim2);

		int i = chromo - 1;
		int indexLow = lim1 / chunkSize;
		int indexHigh = lim2 / chunkSize;

		for (int j = indexLow; j < indexHigh; j++) {
			LOGGER.info("-------------------------------------------------");
			LOGGER.info("Assoc files information for the chromosome " + chromo);
			LOGGER.info("outDir                  : " + outDir.get(testTypeIndex).get(rPanelIndex).get(i));
			LOGGER.info(
					"snptestOutFile    : " + snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(indexLow));
			LOGGER.info(
					"snptestLogFile    : " + snptestLogFile.get(testTypeIndex).get(rPanelIndex).get(i).get(indexLow));
			LOGGER.info("-------------------------------------------------");
		}
	}

	private void checkChromoIndex(int chromo) {
		if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
			LOGGER.fatal("[AssocFiles] Error, chromosome " + chromo + "does not exist");
			System.exit(1);
		}
	}

	private void checkLimits(int chromo, int lim1, int lim2) {
		if ((lim1 < MIN_LIMIT) || (lim2 > MAX_LIMIT)) {
			LOGGER.fatal("[AssocFiles] Error, Chunk " + lim1 + "_" + lim2 + "does not exist for chromosome " + chromo);
			System.exit(1);
		}
	}

}
