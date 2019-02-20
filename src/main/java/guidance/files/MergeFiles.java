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

public class MergeFiles {

	// Logger
	private static final Logger LOGGER = LogManager.getLogger("Console");

	private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;

	// File names and extensions
	private static final String EXT_TXT_GZ = ".txt.gz";
	private static final String REDUCE_FILTERED_FILENAME = "_reduce_filtered_file_";
	private static final String REDUCE_CONDENSED_FILENAME = "_reduce_condensed_file_";
	private static final String FILTER_BY_ALL_FILENAME = "_filteredByAll_chr_";
	private static final String UNCOMPRESSED_FILE = "uncompressed";

	private ArrayList<ArrayList<ArrayList<String>>> testTypeReducedOutDir = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedFilteredFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedFilteredMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedFilteredFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedCondensedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedCondensedMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedCondensedFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeFilteredByAllFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeFilteredByAllMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeFilteredByAllFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCondensedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCondensedMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCondensedFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalCondensedFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalFilteredByAllFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalFilteredByAllMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalFilteredByAllFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<Integer>> testTypeAdditionalCondensedIndex = new ArrayList<>();
	private ArrayList<ArrayList<Integer>> testTypeAdditionalFilteredByAllIndex = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedFilteredFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedFilteredMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedFilteredFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedCondensedFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedFilteredByAllFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedFilteredByAllMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedFilteredByAllFemalesFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedCondensedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedCondensedMalesFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedCondensedFemalesFile = new ArrayList<>();

	private int startChr = 0;
	private int endChr = 0;

	/**
	 * New MergeFiles instance
	 * 
	 * @param parsingArgs
	 * @param generalChromoInfo
	 * @param baseOutDir
	 * @param refPanels
	 */
	public MergeFiles(ParseCmdLine parsingArgs, String baseOutDir, List<String> refPanels) {
		this.startChr = parsingArgs.getStart();
		this.endChr = parsingArgs.getEnd();

		int chunkSize = parsingArgs.getChunkSize();
		int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();

		// We create the first directory name: the cohort directory
		String mixedCohort = parsingArgs.getCohort();
		for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
			String testTypeName = parsingArgs.getTestTypeName(tt);
			String testTypeOutDir = baseOutDir + File.separator + "associations" + File.separator + testTypeName;
			String testTypeOutDir2 = baseOutDir + File.separator + "associations" + File.separator + testTypeName
					+ File.separator + mixedCohort + "_combined_panels";

			String rPanel = null;
			for (int j = 0; j < refPanels.size(); j++) {
				rPanel = refPanels.get(j);
				testTypeOutDir2 = testTypeOutDir2 + "_" + rPanel;
			}

			ArrayList<ArrayList<String>> rpanelListOutDir = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedFile = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedFilteredFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedFilteredMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedFilteredFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedCondensedMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedCondensedFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelFilteredByAllFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelFilteredByAllMalesFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelFilteredByAllFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCondensedMalesFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCondensedFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelAdditionalCondensedFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelAdditionalFilteredByAllFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelAdditionalFilteredByAllMalesFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelAdditionalFilteredByAllFemalesFile = new ArrayList<>();

			ArrayList<Integer> rpanelAdditionalCondensedIndex = new ArrayList<>();
			ArrayList<Integer> rpanelAdditionalFilteredByAllIndex = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedFilteredFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedFilteredMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedFilteredFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedCondensedMalesFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedCondensedFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelCombinedFilteredByAllFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCombinedFilteredByAllMalesFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCombinedFilteredByAllFemalesFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelCombinedCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCombinedCondensedMalesFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCombinedCondensedFemalesFile = new ArrayList<>();

			for (int j = 0; j < refPanels.size(); j++) {
				rPanel = refPanels.get(j);
				String rpanelOutDir = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel;
				String rpanelOutDirSummary = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel
						+ File.separator + "summary";

				ArrayList<String> chromoListOutDir = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedFile = new ArrayList<>();

				ArrayList<ArrayList<GenericFile>> chromoListReducedFilteredFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedFilteredMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedFilteredFemalesFile = new ArrayList<>();

				ArrayList<ArrayList<GenericFile>> chromoListReducedCondensedFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedCondensedMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedCondensedFemalesFile = new ArrayList<>();

				ArrayList<GenericFile> chromoFilteredByAllFile = new ArrayList<>();
				ArrayList<GenericFile> chromoFilteredByAllMalesFile = new ArrayList<>();
				ArrayList<GenericFile> chromoFilteredByAllFemalesFile = new ArrayList<>();

				ArrayList<GenericFile> chromoCondensedFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCondensedMalesFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCondensedFemalesFile = new ArrayList<>();

				ArrayList<GenericFile> additionalCondensedFile = new ArrayList<>();

				ArrayList<GenericFile> additionalFilteredByAllFile = new ArrayList<>();
				ArrayList<GenericFile> additionalFilteredByAllMalesFile = new ArrayList<>();
				ArrayList<GenericFile> additionalFilteredByAllFemalesFile = new ArrayList<>();

				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedFilteredFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedFilteredMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedFilteredFemalesFile = new ArrayList<>();

				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedCondensedFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedCondensedMalesFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedCondensedFemalesFile = new ArrayList<>();

				ArrayList<GenericFile> chromoCombinedFilteredByAllFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCombinedFilteredByAllMalesFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCombinedFilteredByAllFemalesFile = new ArrayList<>();

				ArrayList<GenericFile> chromoCombinedCondensedFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCombinedCondensedMalesFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCombinedCondensedFemalesFile = new ArrayList<>();

				boolean addExtraCount = false;
				for (int i = this.startChr; i <= this.endChr; i++) {
					int chromo = i;

					if (i != this.startChr) {
						addExtraCount = true;
					}

					String tmpChrDir = rpanelOutDir + File.separator + "Chr_" + chromo;
					chromoListOutDir.add(tmpChrDir);

					int maxSize = ChromoInfo.getMaxSize(chromo);
					int total_chunks = maxSize / chunkSize;
					if (maxSize % chunkSize != 0) {
						total_chunks++;
					}

					// ArrayList<GenericFile> listReducedFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedFilteredFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedFilteredMalesFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedFilteredFemalesFile = new ArrayList<>();

					ArrayList<GenericFile> listReducedCondensedFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedCondensedMalesFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedCondensedFemalesFile = new ArrayList<>();

					ArrayList<GenericFile> listCombinedReducedFilteredFile = new ArrayList<>();
					ArrayList<GenericFile> listCombinedReducedFilteredMalesFile = new ArrayList<>();
					ArrayList<GenericFile> listCombinedReducedFilteredFemalesFile = new ArrayList<>();

					ArrayList<GenericFile> listCombinedReducedCondensedFile = new ArrayList<>();
					ArrayList<GenericFile> listCombinedReducedCondensedMalesFile = new ArrayList<>();
					ArrayList<GenericFile> listCombinedReducedCondensedFemalesFile = new ArrayList<>();

					// Now we have to create the rest of file names that will be used to reduce
					// files
					int counter = 0;
					for (int deep = 0; deep < total_chunks - 1; deep++) {
						if (chromo == 23) {

							String tmpReducedFilteredMalesFileName = "chr_" + chromo + "_" + testTypeName + "_males_"
									+ rPanel + REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
							String tmpReducedCondensedMalesFileName = "chr_" + chromo + "_" + testTypeName + "_males_"
									+ rPanel + REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;

							GenericFile myReducedFilteredMalesFile = new GenericFile(tmpChrDir,
									tmpReducedFilteredMalesFileName, UNCOMPRESSED_FILE, "none");
							GenericFile myReducedCondensedMalesFile = new GenericFile(tmpChrDir,
									tmpReducedCondensedMalesFileName, UNCOMPRESSED_FILE, "none");

							String tmpReducedFilteredFemalesFileName = "chr_" + chromo + "_" + testTypeName
									+ "_females_" + rPanel + REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
							String tmpReducedCondensedFemalesFileName = "chr_" + chromo + "_" + testTypeName
									+ "_females_" + rPanel + REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;

							GenericFile myReducedFilteredFemalesFile = new GenericFile(tmpChrDir,
									tmpReducedFilteredFemalesFileName, UNCOMPRESSED_FILE, "none");
							GenericFile myReducedCondensedFemalesFile = new GenericFile(tmpChrDir,
									tmpReducedCondensedFemalesFileName, UNCOMPRESSED_FILE, "none");

							listReducedFilteredMalesFile.add(myReducedFilteredMalesFile);
							listReducedCondensedMalesFile.add(myReducedCondensedMalesFile);

							listReducedFilteredFemalesFile.add(myReducedFilteredFemalesFile);
							listReducedCondensedFemalesFile.add(myReducedCondensedFemalesFile);
						} else {
							String tmpReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
							String tmpReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
									+ REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;

							GenericFile myReducedFilteredFile = new GenericFile(tmpChrDir, tmpReducedFilteredFileName,
									UNCOMPRESSED_FILE, "none");
							GenericFile myReducedCondensedFile = new GenericFile(tmpChrDir, tmpReducedCondensedFileName,
									UNCOMPRESSED_FILE, "none");

							listReducedFilteredFile.add(myReducedFilteredFile);
							listReducedCondensedFile.add(myReducedCondensedFile);
						}

						// listReducedFile.add(myReducedFile);

						if (j == 0) {
							if (chromo == 23) {
								String tmpCombinedReducedFilteredMalesFileName = "chr_" + chromo + "_males_"
										+ testTypeName + REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
								GenericFile myCombinedReducedFilteredMalesFile = new GenericFile(testTypeOutDir2,
										tmpCombinedReducedFilteredMalesFileName, UNCOMPRESSED_FILE, "none");

								listCombinedReducedFilteredMalesFile.add(myCombinedReducedFilteredMalesFile);

								String tmpCombinedReducedCondensedMalesFileName = "chr_" + chromo + "_males_"
										+ testTypeName + REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;
								GenericFile myCombinedReducedCondensedMalesFile = new GenericFile(testTypeOutDir2,
										tmpCombinedReducedCondensedMalesFileName, UNCOMPRESSED_FILE, "none");
								listCombinedReducedCondensedMalesFile.add(myCombinedReducedCondensedMalesFile);

								String tmpCombinedReducedFilteredFemalesFileName = "chr_" + chromo + "_females_"
										+ testTypeName + REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
								GenericFile myCombinedReducedFilteredFemalesFile = new GenericFile(testTypeOutDir2,
										tmpCombinedReducedFilteredFemalesFileName, UNCOMPRESSED_FILE, "none");

								listCombinedReducedFilteredFemalesFile.add(myCombinedReducedFilteredFemalesFile);

								String tmpCombinedReducedCondensedFemalesFileName = "chr_" + chromo + "_females_"
										+ testTypeName + REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;
								GenericFile myCombinedReducedCondensedFemalesFile = new GenericFile(testTypeOutDir2,
										tmpCombinedReducedCondensedFemalesFileName, UNCOMPRESSED_FILE, "none");
								listCombinedReducedCondensedFemalesFile.add(myCombinedReducedCondensedFemalesFile);
							} else {
								String tmpCombinedReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName
										+ REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
								GenericFile myCombinedReducedFilteredFile = new GenericFile(testTypeOutDir2,
										tmpCombinedReducedFilteredFileName, UNCOMPRESSED_FILE, "none");

								listCombinedReducedFilteredFile.add(myCombinedReducedFilteredFile);

								String tmpCombinedReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName
										+ REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;
								GenericFile myCombinedReducedCondensedFile = new GenericFile(testTypeOutDir2,
										tmpCombinedReducedCondensedFileName, UNCOMPRESSED_FILE, "none");
								listCombinedReducedCondensedFile.add(myCombinedReducedCondensedFile);
							}
						}
						counter++;
					}

					if (addExtraCount) {
						if (chromo == 23) {
							String tmpCombinedReducedFilteredMalesFileName = "chr_" + chromo + "_males_" + testTypeName
									+ REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
							String tmpCombinedReducedCondensedMalesFileName = "chr_" + chromo + "_males_" + testTypeName
									+ REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;
							GenericFile myCombinedReducedFilteredMalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedReducedFilteredMalesFileName, UNCOMPRESSED_FILE, "none");
							GenericFile myCombinedReducedCondensedMalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedReducedCondensedMalesFileName, UNCOMPRESSED_FILE, "none");

							listCombinedReducedFilteredMalesFile.add(myCombinedReducedFilteredMalesFile);
							listCombinedReducedCondensedFemalesFile.add(myCombinedReducedCondensedMalesFile);

							String tmpCombinedReducedFilteredFemalesFileName = "chr_" + chromo + "_females_"
									+ testTypeName + REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
							String tmpCombinedReducedCondensedFemalesFileName = "chr_" + chromo + "_females_"
									+ testTypeName + REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;
							GenericFile myCombinedReducedFilteredFemalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedReducedFilteredFemalesFileName, UNCOMPRESSED_FILE, "none");
							GenericFile myCombinedReducedCondensedFemalesFile = new GenericFile(testTypeOutDir2,
									tmpCombinedReducedCondensedFemalesFileName, UNCOMPRESSED_FILE, "none");

							listCombinedReducedFilteredFemalesFile.add(myCombinedReducedFilteredFemalesFile);
							listCombinedReducedCondensedFemalesFile.add(myCombinedReducedCondensedFemalesFile);

						} else {
							String tmpCombinedReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName
									+ REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
							String tmpCombinedReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName
									+ REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;
							GenericFile myCombinedReducedFilteredFile = new GenericFile(testTypeOutDir2,
									tmpCombinedReducedFilteredFileName, UNCOMPRESSED_FILE, "none");
							GenericFile myCombinedReducedCondensedFile = new GenericFile(testTypeOutDir2,
									tmpCombinedReducedCondensedFileName, UNCOMPRESSED_FILE, "none");

							listCombinedReducedFilteredFile.add(myCombinedReducedFilteredFile);
							listCombinedReducedCondensedFile.add(myCombinedReducedCondensedFile);

						}

						addExtraCount = false;
					}

					if (chromo == 23) {
						String tmpFilteredByAllMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_filtered_by_maf_info_hwe_males.txt.gz";
						GenericFile myFilteredByAllMalesFile = new GenericFile(tmpChrDir, tmpFilteredByAllMalesFileName,
								UNCOMPRESSED_FILE, "none");
						chromoFilteredByAllMalesFile.add(myFilteredByAllMalesFile);

						String tmpFilteredByAllFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_filtered_by_maf_info_hwe_females.txt.gz";
						GenericFile myFilteredByAllFemalesFile = new GenericFile(tmpChrDir,
								tmpFilteredByAllFemalesFileName, UNCOMPRESSED_FILE, "none");
						chromoFilteredByAllFemalesFile.add(myFilteredByAllFemalesFile);

						String tmpCombinedCondensedMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_combined_condensed_males.txt.gz";
						GenericFile myCombinedCondensedMalesFile = new GenericFile(testTypeOutDir2,
								tmpCombinedCondensedMalesFileName, UNCOMPRESSED_FILE, "none");
						chromoCombinedCondensedMalesFile.add(myCombinedCondensedMalesFile);

						String tmpCombinedCondensedFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_combined_condensed_females.txt.gz";
						GenericFile myCombinedCondensedFemalesFile = new GenericFile(testTypeOutDir2,
								tmpCombinedCondensedFemalesFileName, UNCOMPRESSED_FILE, "none");
						chromoCombinedCondensedFemalesFile.add(myCombinedCondensedFemalesFile);

						chromoListReducedFilteredMalesFile.add(listReducedFilteredMalesFile);
						chromoListReducedFilteredFemalesFile.add(listReducedFilteredFemalesFile);

						chromoListReducedCondensedMalesFile.add(listReducedCondensedMalesFile);
						chromoListReducedCondensedFemalesFile.add(listReducedCondensedFemalesFile);

						chromoListCombinedReducedFilteredMalesFile.add(listCombinedReducedFilteredMalesFile);
						chromoListCombinedReducedFilteredFemalesFile.add(listCombinedReducedFilteredFemalesFile);

						chromoListCombinedReducedCondensedMalesFile.add(listCombinedReducedCondensedMalesFile);
						chromoListCombinedReducedCondensedFemalesFile.add(listCombinedReducedCondensedFemalesFile);

						String tmpCombinedFilteredByAllMalesFileName = "chr_" + chromo + "_" + testTypeName + "_"
								+ rPanel + "_males_combined_filtered_by_maf_info_hwe.txt.gz";

						GenericFile myCombinedFilteredByAllMalesFile = new GenericFile(testTypeOutDir2,
								tmpCombinedFilteredByAllMalesFileName, UNCOMPRESSED_FILE, "none");

						chromoCombinedFilteredByAllMalesFile.add(myCombinedFilteredByAllMalesFile);

						String tmpCombinedFilteredByAllFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_"
								+ rPanel + "_females_combined_filtered_by_maf_info_hwe.txt.gz";

						GenericFile myCombinedFilteredByAllFemalesFile = new GenericFile(testTypeOutDir2,
								tmpCombinedFilteredByAllFemalesFileName, UNCOMPRESSED_FILE, "none");

						chromoCombinedFilteredByAllFemalesFile.add(myCombinedFilteredByAllFemalesFile);

						String tmpCondensedMalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_condensed_males.txt.gz";
						GenericFile myCondensedMalesFile = new GenericFile(tmpChrDir, tmpCondensedMalesFileName,
								UNCOMPRESSED_FILE, "none");
						chromoCondensedMalesFile.add(myCondensedMalesFile);

						String tmpCondensedFemalesFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_condensed_females.txt.gz";
						GenericFile myCondensedFemalesFile = new GenericFile(tmpChrDir, tmpCondensedFemalesFileName,
								UNCOMPRESSED_FILE, "none");
						chromoCondensedFemalesFile.add(myCondensedFemalesFile);
						
						chromoCondensedFile.add(myCondensedMalesFile);
						chromoCondensedFile.add(myCondensedFemalesFile);

					} else {
						String tmpFilteredByAllFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_filtered_by_maf_info_hwe.txt.gz";
						GenericFile myFilteredByAllFile = new GenericFile(tmpChrDir, tmpFilteredByAllFileName,
								UNCOMPRESSED_FILE, "none");
						chromoFilteredByAllFile.add(myFilteredByAllFile);

						String tmpCombinedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_combined_condensed.txt.gz";
						GenericFile myCombinedCondensedFile = new GenericFile(testTypeOutDir2,
								tmpCombinedCondensedFileName, UNCOMPRESSED_FILE, "none");
						chromoCombinedCondensedFile.add(myCombinedCondensedFile);

						chromoListReducedFilteredFile.add(listReducedFilteredFile);
						chromoListReducedCondensedFile.add(listReducedCondensedFile);
						chromoListCombinedReducedFilteredFile.add(listCombinedReducedFilteredFile);
						chromoListCombinedReducedCondensedFile.add(listCombinedReducedCondensedFile);

						String tmpCombinedFilteredByAllFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_combined_filtered_by_maf_info_hwe.txt.gz";

						GenericFile myCombinedFilteredByAllFile = new GenericFile(testTypeOutDir2,
								tmpCombinedFilteredByAllFileName, UNCOMPRESSED_FILE, "none");

						chromoCombinedFilteredByAllFile.add(myCombinedFilteredByAllFile);

						String tmpCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_condensed.txt.gz";
						GenericFile myCondensedFile = new GenericFile(tmpChrDir, tmpCondensedFileName,
								UNCOMPRESSED_FILE, "none");
						chromoCondensedFile.add(myCondensedFile);
					}

				} // End for Chromo
				rpanelListOutDir.add(chromoListOutDir);
				rpanelReducedFile.add(chromoListReducedFile);

				rpanelReducedFilteredFile.add(chromoListReducedFilteredFile);
				rpanelReducedFilteredMalesFile.add(chromoListReducedFilteredMalesFile);
				rpanelReducedFilteredFemalesFile.add(chromoListReducedFilteredFemalesFile);

				rpanelReducedCondensedFile.add(chromoListReducedCondensedFile);
				rpanelReducedCondensedMalesFile.add(chromoListReducedCondensedMalesFile);
				rpanelReducedCondensedFemalesFile.add(chromoListReducedCondensedFemalesFile);

				rpanelFilteredByAllFile.add(chromoFilteredByAllFile);
				rpanelFilteredByAllMalesFile.add(chromoFilteredByAllMalesFile);
				rpanelFilteredByAllFemalesFile.add(chromoFilteredByAllFemalesFile);

				rpanelCondensedFile.add(chromoCondensedFile);
				rpanelCondensedMalesFile.add(chromoCondensedMalesFile);
				rpanelCondensedFemalesFile.add(chromoCondensedFemalesFile);

				// Only for all the panels combined
				if (j == 0) {
					rpanelCombinedReducedFilteredFile.add(chromoListCombinedReducedFilteredFile);
					rpanelCombinedReducedFilteredMalesFile.add(chromoListCombinedReducedFilteredMalesFile);
					rpanelCombinedReducedFilteredFemalesFile.add(chromoListCombinedReducedFilteredFemalesFile);

					rpanelCombinedReducedCondensedFile.add(chromoListCombinedReducedCondensedFile);
					rpanelCombinedReducedCondensedMalesFile.add(chromoListCombinedReducedCondensedMalesFile);
					rpanelCombinedReducedCondensedFemalesFile.add(chromoListCombinedReducedCondensedFemalesFile);

					rpanelCombinedFilteredByAllFile.add(chromoCombinedFilteredByAllFile);
					rpanelCombinedFilteredByAllMalesFile.add(chromoCombinedFilteredByAllMalesFile);
					rpanelCombinedFilteredByAllFemalesFile.add(chromoCombinedFilteredByAllFemalesFile);

					rpanelCombinedCondensedFile.add(chromoCombinedCondensedFile);
					rpanelCombinedCondensedMalesFile.add(chromoCombinedCondensedMalesFile);
					rpanelCombinedCondensedFemalesFile.add(chromoCombinedCondensedFemalesFile);
				}

				// Here we have to create an additional list of condensed files that will be
				// used when we execute
				// jointCondensedFiles Task, for all chromosomes.
				// The number of additional files is the number of chromosomes minus 1.
				int addCondensed = 0;
				int limit = 0;
				if(this.endChr == 23) {
					limit = this.endChr + 1;
				} else {
					limit = this.endChr;
				}
				
				for (int deep = this.startChr; deep < limit; deep++) {
					String tmpAdditionalCondensedFileName = null;
					if (deep == (limit - 1)) {
						tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_chr_" + this.startChr
								+ "_to_" + this.endChr + EXT_TXT_GZ;
					} else {
						tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_" + addCondensed
								+ EXT_TXT_GZ;
					}
					GenericFile myAdditionalCondensedFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalCondensedFileName, UNCOMPRESSED_FILE, "none");
					additionalCondensedFile.add(myAdditionalCondensedFile);
					// LOGGER.info("\t[MergeFiles.java] " + tmpAdditionalCondensedFile);

					addCondensed++;
				}

				if (this.startChr == this.endChr) {
					String tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_chr_"
							+ this.startChr + EXT_TXT_GZ;
					GenericFile myAdditionalCondensedFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalCondensedFileName, UNCOMPRESSED_FILE, "none");
					additionalCondensedFile.add(myAdditionalCondensedFile);
					// LOGGER.info("\t[MergeFiles.java] only " + tmpAdditionalCondensedFile);
					addCondensed++;
				}

				rpanelAdditionalCondensedIndex.add(addCondensed);

				// Here we have to create an additional list of filteredByAll files that will be
				// used when we execute
				// jointFilteredByAllFile task for all chromosomes.
				// Unlike the previous case with condensed files, we can not include chromosome
				// 23. (Chr 23 format of filteredByAllFile is different to the rest of
				// chromosomes (because of snptest).

				// The number of additional files is the number of chromosomes minus 1.
				int addFiltered = 0;
				int endChrNormal = this.endChr;
				if (this.startChr < this.endChr) {
					if (this.endChr != 23) {
						endChrNormal = this.endChr;
					} else {
						endChrNormal = this.endChr - 1;
					}
				}

				for (int deep = this.startChr; deep < endChrNormal; deep++) {
					String tmpAdditionalFilteredByAllFileName = null;
					if (deep == (endChrNormal - 1)) {
						tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + FILTER_BY_ALL_FILENAME
								+ this.startChr + "_to_" + endChrNormal + EXT_TXT_GZ;
					} else {
						tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + "_filteredByAll_"
								+ addFiltered + EXT_TXT_GZ;
					}
					GenericFile myAdditionalFilteredByAllFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalFilteredByAllFileName, UNCOMPRESSED_FILE, "none");
					additionalFilteredByAllFile.add(myAdditionalFilteredByAllFile);
					// LOGGER.info("\t[MergeFiles.java] " + tmpAdditionalFilteredByAllFile);

					addFiltered++;
				}

				if (this.startChr == endChrNormal) {
					String tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + FILTER_BY_ALL_FILENAME
							+ this.startChr + EXT_TXT_GZ;
					GenericFile myAdditionalFilteredByAllFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalFilteredByAllFileName, UNCOMPRESSED_FILE, "none");
					additionalFilteredByAllFile.add(myAdditionalFilteredByAllFile);
					addFiltered++;
				}

				rpanelAdditionalFilteredByAllIndex.add(addFiltered);
				rpanelAdditionalCondensedFile.add(additionalCondensedFile);
				rpanelAdditionalFilteredByAllFile.add(additionalFilteredByAllFile);

				// If there is chr 23:
				if (this.endChr == 23) {
					String tmpAdditionalFilteredByAllMalesFileName = testTypeName + "_" + rPanel
							+ FILTER_BY_ALL_FILENAME + this.endChr + "_males" + EXT_TXT_GZ;
					GenericFile myAdditionalFilteredByAllMalesFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalFilteredByAllMalesFileName, UNCOMPRESSED_FILE, "none");

					additionalFilteredByAllMalesFile.add(myAdditionalFilteredByAllMalesFile);
					// LOGGER.info("\t[MergeFiles.java] " + tmpAdditionalFilteredByAllXFile);

					rpanelAdditionalFilteredByAllMalesFile.add(additionalFilteredByAllMalesFile);

					String tmpAdditionalFilteredByAllFemalesFileName = testTypeName + "_" + rPanel
							+ FILTER_BY_ALL_FILENAME + this.endChr + "_females" + EXT_TXT_GZ;
					GenericFile myAdditionalFilteredByAllFemalesFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalFilteredByAllFemalesFileName, UNCOMPRESSED_FILE, "none");

					additionalFilteredByAllFemalesFile.add(myAdditionalFilteredByAllFemalesFile);
					// LOGGER.info("\t[MergeFiles.java] " + tmpAdditionalFilteredByAllXFile);

					rpanelAdditionalFilteredByAllFemalesFile.add(additionalFilteredByAllFemalesFile);
				}
			} // End of for refPanels

			// Now we have to build the list of reduced files for the type of Test. We store
			// this list
			this.testTypeReducedOutDir.add(rpanelListOutDir);
			this.testTypeReducedFile.add(rpanelReducedFile);

			this.testTypeReducedFilteredFile.add(rpanelReducedFilteredFile);
			this.testTypeReducedFilteredMalesFile.add(rpanelReducedFilteredMalesFile);
			this.testTypeReducedFilteredFemalesFile.add(rpanelReducedFilteredFemalesFile);

			this.testTypeReducedCondensedFile.add(rpanelReducedCondensedFile);
			this.testTypeReducedCondensedMalesFile.add(rpanelReducedCondensedMalesFile);
			this.testTypeReducedCondensedFemalesFile.add(rpanelReducedCondensedFemalesFile);

			this.testTypeCombinedReducedFilteredFile.add(rpanelCombinedReducedFilteredFile);
			this.testTypeCombinedReducedFilteredMalesFile.add(rpanelCombinedReducedFilteredMalesFile);
			this.testTypeCombinedReducedFilteredFemalesFile.add(rpanelCombinedReducedFilteredFemalesFile);
			this.testTypeCombinedReducedCondensedFile.add(rpanelCombinedReducedCondensedFile);

			this.testTypeFilteredByAllFile.add(rpanelFilteredByAllFile);
			this.testTypeFilteredByAllFemalesFile.add(rpanelFilteredByAllFemalesFile);
			this.testTypeFilteredByAllMalesFile.add(rpanelFilteredByAllMalesFile);

			this.testTypeCondensedFile.add(rpanelCondensedFile);
			this.testTypeCondensedFemalesFile.add(rpanelCondensedFemalesFile);
			this.testTypeCondensedMalesFile.add(rpanelCondensedMalesFile);

			this.testTypeCombinedFilteredByAllFile.add(rpanelCombinedFilteredByAllFile);
			this.testTypeCombinedFilteredByAllMalesFile.add(rpanelCombinedFilteredByAllMalesFile);
			this.testTypeCombinedFilteredByAllFemalesFile.add(rpanelCombinedFilteredByAllFemalesFile);

			this.testTypeCombinedCondensedFile.add(rpanelCombinedCondensedFile);
			this.testTypeCombinedCondensedMalesFile.add(rpanelCombinedCondensedMalesFile);
			this.testTypeCombinedCondensedFemalesFile.add(rpanelCombinedCondensedFemalesFile);

			this.testTypeAdditionalCondensedFile.add(rpanelAdditionalCondensedFile);
			this.testTypeAdditionalFilteredByAllFile.add(rpanelAdditionalFilteredByAllFile);

			if (this.endChr == 23) {
				this.testTypeAdditionalFilteredByAllMalesFile.add(rpanelAdditionalFilteredByAllMalesFile);
				this.testTypeAdditionalFilteredByAllFemalesFile.add(rpanelAdditionalFilteredByAllFemalesFile);
			}

			this.testTypeAdditionalCondensedIndex.add(rpanelAdditionalCondensedIndex);
			this.testTypeAdditionalFilteredByAllIndex.add(rpanelAdditionalFilteredByAllIndex);

		} // End of for test types
	}

	/**
	 * Method to access mergedGenDir information
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
		return this.testTypeReducedOutDir.get(testTypeIndex).get(rPanelIndex).get(i);
	}

	/**
	 * Method to access reducedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFile(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access reducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFilteredFile(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access reducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFilteredMalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeReducedFilteredMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFilteredFemalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeReducedFilteredFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedCondensedFile(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
	}

	/**
	 * Method to access reducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedCondensedMalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeReducedCondensedMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedCondensedFemalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeReducedCondensedFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCombinedReducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedFilteredFile(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCombinedReducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedFilteredMalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeCombinedReducedFilteredMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCombinedReducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedFilteredFemalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeCombinedReducedFilteredFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCombinedReducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedCondensedFile(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCombinedReducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedCondensedMalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedCombinedReducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedCondensedFemalesFile(int testTypeIndex, int rPanelIndex, int index) {
		return this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(0).get(index)
				.getFullName();
	}

	/**
	 * Method to access reducedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getTheLastReducedFile(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		int lastIndex = this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size() - 1;
		return this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(lastIndex).getFullName();
	}

	/**
	 * Method to access the index of the reducedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public int getTheLastReducedFileIndex(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		int lastIndex = this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size() - 1;
		return lastIndex;
	}

	/**
	 * Method to access print theLastReducedFiles
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 */
	public void printTheLastReducedFiles(int testTypeIndex, int rPanelIndex, int chromo) {
		LOGGER.info("\t[MergeFiles] theLastReducedFiles ARE:");
		for (int hh = 0; hh <= testTypeIndex; hh++) {
			for (int kk = 0; kk <= rPanelIndex; kk++) {
				for (int j = 0; j <= chromo - this.startChr; j++) {
					LOGGER.info("\t[MergeFiles] " + testTypeReducedFile.get(rPanelIndex).get(kk).get(j));
				}
			}
		}
	}

	/**
	 * Method to access filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getFilteredByAllFile(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
	}

	/**
	 * Method to access filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getFilteredByAllMalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeFilteredByAllMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getFilteredByAllFemalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeFilteredByAllFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedFilteredByAllFile(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
	}

	/**
	 * Method to access filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedFilteredByAllMalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeCombinedFilteredByAllMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedFilteredByAllFemalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeCombinedFilteredByAllFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCondensedFile(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
	}

	/**
	 * Method to access condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCondensedMalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeCondensedMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCondensedFemalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeCondensedFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedCondensedFile(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeCombinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
	}

	/**
	 * Method to access condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedCondensedMalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeCombinedCondensedMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedCondensedFemalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeCombinedCondensedFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access additionalCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalCondensedFile(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
	}

	/**
	 * Method to access additionalFilteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllFile(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
	}

	/**
	 * Method to access additionalFilteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllMalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeAdditionalFilteredByAllMalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access additionalFilteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllFemalesFile(int testTypeIndex, int rPanelIndex) {
		return this.testTypeAdditionalFilteredByAllFemalesFile.get(testTypeIndex).get(rPanelIndex).get(0).getFullName();
	}

	/**
	 * Method to access additionalFilteredByAllXFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllXMalesFile(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllMalesFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllMalesFile.get(testTypeIndex).get(rPanelIndex).get(index)
				.getFullName();
	}

	/**
	 * Method to access additionalFilteredByAllXFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllXFemalesFile(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllFemalesFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllFemalesFile.get(testTypeIndex).get(rPanelIndex).get(index)
				.getFullName();
	}

	/**
	 * Method to access the last CondensedFile of each testType and rPanel
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @return
	 */
	public String getFinalCondensedFile(int testTypeIndex, int rPanelIndex) {
		int lastIndex = this.testTypeAdditionalCondensedIndex.get(testTypeIndex).get(rPanelIndex);
		return this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(lastIndex - 1)
				.getFullName();
	}

	/**
	 * Method to access the last FilteredByAllFile of each testType and rPanel
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @return
	 */
	public String getFinalFilteredByAllFile(int testTypeIndex, int rPanelIndex) {
		int lastIndex = this.testTypeAdditionalFilteredByAllIndex.get(testTypeIndex).get(rPanelIndex);
		return this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(lastIndex - 1)
				.getFullName();
	}

	private void checkChromoIndex(int chromo) {
		if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
			LOGGER.fatal("[MergeFiles] Error, chromosome " + chromo + " does not exist");
			System.exit(1);
		}
	}

	private void checkIndex(int index, int maxIndex, int chromo) {
		if (index >= maxIndex) {
			LOGGER.fatal("[MergeFiles] Error, the index " + index + " is greater than the maximum " + maxIndex
					+ " in chromosome " + chromo);
			System.exit(1);
		}
	}

}
