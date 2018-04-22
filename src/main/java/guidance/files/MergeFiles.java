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
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeReducedCondensedFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeFilteredByAllFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCondensedFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalCondensedFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalFilteredByAllFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeAdditionalFilteredByAllXFile = new ArrayList<>();

	private ArrayList<ArrayList<Integer>> testTypeAdditionalCondensedIndex = new ArrayList<>();
	private ArrayList<ArrayList<Integer>> testTypeAdditionalFilteredByAllIndex = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedFilteredFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedCondensedFile = new ArrayList<>();

	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedFilteredByAllFile = new ArrayList<>();
	private ArrayList<ArrayList<ArrayList<GenericFile>>> testTypeCombinedCondensedFile = new ArrayList<>();

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
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelReducedCondensedFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelFilteredByAllFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCondensedFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelAdditionalCondensedFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelAdditionalFilteredByAllFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelAdditionalFilteredByAllXFile = new ArrayList<>();

			ArrayList<Integer> rpanelAdditionalCondensedIndex = new ArrayList<>();
			ArrayList<Integer> rpanelAdditionalFilteredByAllIndex = new ArrayList<>();

			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedFilteredFile = new ArrayList<>();
			ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedCondensedFile = new ArrayList<>();

			ArrayList<ArrayList<GenericFile>> rpanelCombinedFilteredByAllXFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCombinedFilteredByAllFile = new ArrayList<>();
			ArrayList<ArrayList<GenericFile>> rpanelCombinedCondensedFile = new ArrayList<>();

			for (int j = 0; j < refPanels.size(); j++) {
				rPanel = refPanels.get(j);
				String rpanelOutDir = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel;
				String rpanelOutDirSummary = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel
						+ File.separator + "summary";

				ArrayList<String> chromoListOutDir = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedFilteredFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListReducedCondensedFile = new ArrayList<>();

				ArrayList<GenericFile> chromoFilteredByAllFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCondensedFile = new ArrayList<>();

				ArrayList<GenericFile> additionalCondensedFile = new ArrayList<>();
				ArrayList<GenericFile> additionalFilteredByAllFile = new ArrayList<>();
				ArrayList<GenericFile> additionalFilteredByAllXFile = new ArrayList<>();

				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedFilteredFile = new ArrayList<>();
				ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedCondensedFile = new ArrayList<>();

				ArrayList<GenericFile> chromoCombinedFilteredByAllXFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCombinedFilteredByAllFile = new ArrayList<>();
				ArrayList<GenericFile> chromoCombinedCondensedFile = new ArrayList<>();

				boolean addExtraCount = false;
				for (int i = this.startChr; i <= this.endChr; i++) {
					int chromo = i;

					if (i != this.startChr) {
						LOGGER.info("mon :: TRUE :: i = " + i);
						addExtraCount = true;
					}

					String tmpChrDir = rpanelOutDir + File.separator + "Chr_" + chromo;
					chromoListOutDir.add(tmpChrDir);

					int maxSize = ChromoInfo.getMaxSize(chromo);
					int total_chunks = maxSize / chunkSize;
					int module = maxSize % chunkSize;
					if (module != 0) {
						total_chunks++;
					}

					ArrayList<GenericFile> listReducedFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedFilteredFile = new ArrayList<>();
					ArrayList<GenericFile> listReducedCondensedFile = new ArrayList<>();

					ArrayList<GenericFile> listCombinedReducedFilteredFile = new ArrayList<>();
					ArrayList<GenericFile> listCombinedReducedCondensedFile = new ArrayList<>();

					// Now we have to create the rest of file names that will be used to reduce
					// files
					int counter = 0;
					for (int deep = 0; deep < total_chunks - 1; deep++) {
						String tmpReducedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ "_reduce_file_" + counter + EXT_TXT_GZ;
						String tmpReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ REDUCE_FILTERED_FILENAME + counter + EXT_TXT_GZ;
						String tmpReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
								+ REDUCE_CONDENSED_FILENAME + counter + EXT_TXT_GZ;

						GenericFile myReducedFile = new GenericFile(tmpChrDir, tmpReducedFileName, UNCOMPRESSED_FILE,
								"none");
						GenericFile myReducedFilteredFile = new GenericFile(tmpChrDir, tmpReducedFilteredFileName,
								UNCOMPRESSED_FILE, "none");
						GenericFile myReducedCondensedFile = new GenericFile(tmpChrDir, tmpReducedCondensedFileName,
								UNCOMPRESSED_FILE, "none");

						listReducedFile.add(myReducedFile);
						listReducedFilteredFile.add(myReducedFilteredFile);
						listReducedCondensedFile.add(myReducedCondensedFile);

						if (j == 0) {
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
						counter++;
					}

					LOGGER.info("mon :: ABANS F :: chr " + chromo + " ::  size = "
							+ listCombinedReducedFilteredFile.size());

					LOGGER.info("mon :: ABANS C :: chr " + chromo + " ::  size = "
							+ listCombinedReducedCondensedFile.size());

					if (addExtraCount) {
						LOGGER.info("mon :: DINS EXTRA :: chr " + chromo + " ::  counter " + counter);

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

						addExtraCount = false;
					}

					LOGGER.info("mon :: DESPRES F :: chr " + chromo + " ::  size = "
							+ listCombinedReducedFilteredFile.size());

					LOGGER.info("mon :: DESPRES C :: chr " + chromo + " ::  size = "
							+ listCombinedReducedCondensedFile.size());

					LOGGER.info("\n\n");

					// Now we have built the list of reduce file for this chromosome. We store this
					// list
					chromoListReducedFile.add(listReducedFile);
					chromoListReducedFilteredFile.add(listReducedFilteredFile);
					chromoListReducedCondensedFile.add(listReducedCondensedFile);

					chromoListCombinedReducedFilteredFile.add(listCombinedReducedFilteredFile);
					chromoListCombinedReducedCondensedFile.add(listCombinedReducedCondensedFile);

					String tmpCombinedFilteredByAllFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
							+ "_combined_filtered_by_maf_info_hwe.txt.gz";
					
					GenericFile myCombinedFilteredByAllFile = new GenericFile(testTypeOutDir2,
							tmpCombinedFilteredByAllFileName, UNCOMPRESSED_FILE, "none");
					
					chromoCombinedFilteredByAllFile.add(myCombinedFilteredByAllFile);
					

					String tmpCombinedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
							+ "_combined_condensed.txt.gz";
					GenericFile myCombinedCondensedFile = new GenericFile(testTypeOutDir2, tmpCombinedCondensedFileName,
							UNCOMPRESSED_FILE, "none");
					chromoCombinedCondensedFile.add(myCombinedCondensedFile);

					String tmpFilteredByAllFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
							+ "_filtered_by_maf_info_hwe.txt.gz";
					GenericFile myFilteredByAllFile = new GenericFile(tmpChrDir, tmpFilteredByAllFileName,
							UNCOMPRESSED_FILE, "none");
					chromoFilteredByAllFile.add(myFilteredByAllFile);

					String tmpCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
							+ "_condensed.txt.gz";
					GenericFile myCondensedFile = new GenericFile(tmpChrDir, tmpCondensedFileName, UNCOMPRESSED_FILE,
							"none");
					chromoCondensedFile.add(myCondensedFile);

				} // End for Chromo
				rpanelListOutDir.add(chromoListOutDir);
				rpanelReducedFile.add(chromoListReducedFile);
				rpanelReducedFilteredFile.add(chromoListReducedFilteredFile);
				rpanelReducedCondensedFile.add(chromoListReducedCondensedFile);
				rpanelFilteredByAllFile.add(chromoFilteredByAllFile);
				rpanelCondensedFile.add(chromoCondensedFile);

				// Only for all the panels combined
				if (j == 0) {
					rpanelCombinedReducedFilteredFile.add(chromoListCombinedReducedFilteredFile);
					rpanelCombinedReducedCondensedFile.add(chromoListCombinedReducedCondensedFile);
					rpanelCombinedFilteredByAllXFile.add(chromoCombinedFilteredByAllXFile);
					rpanelCombinedFilteredByAllFile.add(chromoCombinedFilteredByAllFile);
					rpanelCombinedCondensedFile.add(chromoCombinedCondensedFile);
				}

				// Here we have to create an additional list of condensed files that will be
				// used when we execute
				// jointCondensedFiles Task, for all chromosomes.
				// The number of additional files is the number of chromosomes minus 1.
				int addCondensed = 0;
				for (int deep = this.startChr; deep < this.endChr; deep++) {
					String tmpAdditionalCondensedFileName = null;
					if (deep == (this.endChr - 1)) {
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
				// 23. (Chr 23 format fo
				// filteredByAllFile is different to the rest of chromosomes (thanks to
				// snptest).

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
					String tmpAdditionalFilteredByAllXFileName = testTypeName + "_" + rPanel + FILTER_BY_ALL_FILENAME
							+ this.endChr + EXT_TXT_GZ;
					GenericFile myAdditionalFilteredByAllXFile = new GenericFile(rpanelOutDirSummary,
							tmpAdditionalFilteredByAllXFileName, UNCOMPRESSED_FILE, "none");

					additionalFilteredByAllXFile.add(myAdditionalFilteredByAllXFile);
					// LOGGER.info("\t[MergeFiles.java] " + tmpAdditionalFilteredByAllXFile);

					rpanelAdditionalFilteredByAllXFile.add(additionalFilteredByAllXFile);
				}
			} // End of for refPanels

			// Now we have to build the list of reduced files for the type of Test. We store
			// this list
			this.testTypeReducedOutDir.add(rpanelListOutDir);
			this.testTypeReducedFile.add(rpanelReducedFile);
			this.testTypeReducedFilteredFile.add(rpanelReducedFilteredFile);
			this.testTypeReducedCondensedFile.add(rpanelReducedCondensedFile);

			this.testTypeCombinedReducedFilteredFile.add(rpanelCombinedReducedFilteredFile);
			this.testTypeCombinedReducedCondensedFile.add(rpanelCombinedReducedCondensedFile);

			this.testTypeFilteredByAllFile.add(rpanelFilteredByAllFile);
			this.testTypeCondensedFile.add(rpanelCondensedFile);

			this.testTypeCombinedFilteredByAllFile.add(rpanelCombinedFilteredByAllFile);
			this.testTypeCombinedCondensedFile.add(rpanelCombinedCondensedFile);

			this.testTypeAdditionalCondensedFile.add(rpanelAdditionalCondensedFile);
			this.testTypeAdditionalFilteredByAllFile.add(rpanelAdditionalFilteredByAllFile);

			if (this.endChr == 23) {
				this.testTypeAdditionalFilteredByAllXFile.add(rpanelAdditionalFilteredByAllXFile);
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
	 * Method to access reducedFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFileName(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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
	 * Method to set the finalStatus of the reducedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @param finalStatus
	 */
	public void setReducedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the reducedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
	}

	/**
	 * Method to access reducedFilteredFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFilteredFileName(int testTypeIndex, int rPanelIndex, int chromo, int index) {
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
	 * Method to set the finalStatus of the reducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @param finalStatus
	 */
	public void setReducedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		this.testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the reducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
	}

	/**
	 * Method to access reducedCondensedFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedCondensedFileName(int testTypeIndex, int rPanelIndex, int chromo, int index) {
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
	 * Method to set the finalStatus of the reducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @param finalStatus
	 */
	public void setReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		this.testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the reducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
	}

	/**
	 * Method to access reducedCombinedReducedFilteredFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedFilteredFileName(int testTypeIndex, int rPanelIndex, int chromo, int index) {
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
	 * Method to set the finalStatus of the reducedCombinedReducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @param finalStatus
	 */
	public void setCombinedReducedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		this.testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the reducedCombinedReducedFilteredFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.getFinalStatus();
	}

	/**
	 * Method to access reducedCombinedReducedCondensedFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedCondensedFileName(int testTypeIndex, int rPanelIndex, int chromo, int index) {
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
	 * Method to set the finalStatus of the reducedCombinedReducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @param finalStatus
	 */
	public void setCombinedReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the reducedCombinedReducedCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param index
	 * @return
	 */
	public String getCombinedReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo,
			int index) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).size();
		checkIndex(index, maxIndex, chromo);

		return this.testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index)
				.getFinalStatus();
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
	 * Method to access filterByAllFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getFilterByAllFileName(int testTypeIndex, int rPanelIndex, int chromo) {
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
	 * Method to set the finalStatus of the filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param finalStatus
	 */
	public void setFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		this.testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFinalStatus();
	}	

	/**
	 * Method to access filterByAllFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedFilterByAllFileName(int testTypeIndex, int rPanelIndex, int chromo) {
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
	 * Method to set the finalStatus of the filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param finalStatus
	 */
	public void setCombinedFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the filteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFinalStatus();
	}

	/**
	 * Method to access condensedFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCondensedFileName(int testTypeIndex, int rPanelIndex, int chromo) {
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
	 * Method to set the finalStatus of the condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param finalStatus
	 */
	public void setCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		this.testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFinalStatus();
	}

	/**
	 * Method to access condensedFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedCondensedFileName(int testTypeIndex, int rPanelIndex, int chromo) {
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
	 * Method to set the finalStatus of the condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @param finalStatus
	 */
	public void setCombinedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo,
			String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		this.testTypeCombinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the condensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param chromo
	 * @return
	 */
	public String getCombinedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int i = chromo - this.startChr;

		// Check the index
		int maxIndex = this.testTypeCombinedCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(i, maxIndex, chromo);

		return this.testTypeCombinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFinalStatus();
	}

	/**
	 * Method to access additionalCondensedFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalCondensedFileName(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
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
	 * Method to set the finalStatus of the additionalCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @param finalStatus
	 */
	public void setAdditionalCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int index,
			String finalStatus) {
		// Check the index
		int maxIndex = this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the additionalCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).getFinalStatus();
	}

	/**
	 * Method to access additionalFilteredByAllFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllFileName(int testTypeIndex, int rPanelIndex, int index) {
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
	public String getAdditionalFilteredByAllFile(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
	}

	/**
	 * Method to set the finalStatus of the additionalCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @param finalStatus
	 */
	public void setAdditionalFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int index,
			String finalStatus) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index)
				.setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the additionalFilteredByAllFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).getFinalStatus();
	}

	/**
	 * Method to access additionalFilteredByAllXFileName
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllXFileName(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
	}

	/**
	 * Method to access additionalFilteredByAllXFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllXFile(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
	}

	/**
	 * Method to set the finalStatus of the additionalCondensedFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @param finalStatus
	 */
	public void setAdditionalFilteredByAllXFileFinalStatus(int testTypeIndex, int rPanelIndex, int index,
			String finalStatus) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index)
				.setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the finalStatus of the additionalFilteredByAllXFile
	 * 
	 * @param testTypeIndex
	 * @param rPanelIndex
	 * @param index
	 * @return
	 */
	public String getAdditionalFilteredByAllXFileFinalStatus(int testTypeIndex, int rPanelIndex, int index) {
		// Check the index
		int maxIndex = this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).size();
		checkIndex(index, maxIndex, index);

		return this.testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index)
				.getFinalStatus();
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
