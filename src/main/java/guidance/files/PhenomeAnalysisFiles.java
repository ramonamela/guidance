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

import guidance.utils.ParseCmdLine;

import java.io.File;
import java.util.ArrayList;

public class PhenomeAnalysisFiles {

	private String phenomeAnalysisOutDir = null;

	private GenericFile topHitsAllPhenos = null;

	private ArrayList<GenericFile> phenomeMergedTopHits = new ArrayList<>();

	private GenericFile crossPhenoAll = null;
	private GenericFile crossPhenoRanges = null;
	private GenericFile crossPhenoAssocTop = null;

	/**
	 * New PhenomeAnalysisFiles instance
	 * 
	 * @param parsingArgs
	 * @param baseOutDir
	 * @param refPanels
	 */
	public PhenomeAnalysisFiles(ParseCmdLine parsingArgs, String baseOutDir) {
		String testTypeName = null;

		int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();

		String mixedCohort = parsingArgs.getCohort();

		this.phenomeAnalysisOutDir = baseOutDir + File.separator + "associations" + File.separator + "pheno_analysis"
				+ File.separator + mixedCohort;

		String basePrefixPhenomeName = "phenome_analysis_" + mixedCohort;
		
		topHitsAllPhenos = new GenericFile(basePrefixPhenomeName, "tophits_all_phenotypes.txt", "decompressed", "none");

		crossPhenoAll = new GenericFile(basePrefixPhenomeName, "cross_pheno_all.txt", "decompressed", "none");
		crossPhenoRanges = new GenericFile(basePrefixPhenomeName, "cross_pheno_ranges.txt", "decompressed", "none");
		crossPhenoAssocTop = new GenericFile(basePrefixPhenomeName, "cross_pheno_association_topvariants.txt", "decompressed", "none");

		for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
			testTypeName = parsingArgs.getTestTypeName(tt);
			phenomeMergedTopHits.add(new GenericFile(basePrefixPhenomeName, "tophits_merge_" + testTypeName + ".txt", "decompressed", "none"));
		}
	}

	/**
	 * Method to access correctedPvaluesFile information
	 * 
	 * @return
	 */
	public String getPhenomeAnalysisOutDir() {
		return this.phenomeAnalysisOutDir;
	}

	/**
	 * Method to access topHitsAllPhenos
	 * 
	 * @param index
	 * @return
	 */
	public String getTopHitsAllPhenos() {
		return this.topHitsAllPhenos.getFullName();
	}

	/**
	 * Method to access crossPhenoAll
	 * 
	 * @param index
	 * @return
	 */
	public String getCrossPhenoAll() {
		return this.crossPhenoAll.getFullName();
	}

	/**
	 * Method to access crossPhenoRanges
	 * 
	 * @param index
	 * @return
	 */
	public String getCrossPhenoRanges() {
		return this.crossPhenoRanges.getFullName();
	}

	/**
	 * Method to access crossPhenoAssocTop
	 * 
	 * @param index
	 * @return
	 */
	public String getCrossPhenoAssocTop() {
		return this.crossPhenoAssocTop.getFullName();
	}

	/**
	 * Method to access crossPhenoMergedTop
	 * 
	 * @param index
	 * @return
	 */
	public String getCrossPhenoMergedTop(int index) {
		return this.phenomeMergedTopHits.get(index).getFullName();
	}

}
