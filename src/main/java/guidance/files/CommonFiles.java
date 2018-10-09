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

import java.io.File;
import java.util.ArrayList;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import guidance.utils.ParseCmdLine;

public class CommonFiles {

	// Logger
	private static final Logger LOGGER = LogManager.getLogger("Console");

	private final int MAX_NUMBER_OF_CHROMOSOMES = 23;

	private int startChr;
	private int endChr;

	//private ArrayList<String> outputCasesDir = new ArrayList<String>();
	//private ArrayList<String> outputControlsDir = new ArrayList<String>();
	private ArrayList<String> outputMixedDir = new ArrayList<String>();

	private GenericFile mixedBedFile = null;
	private GenericFile mixedBimFile = null;
	private GenericFile mixedFamFile = null;

	private ArrayList<GenericFile> mixedByChrBedFile = new ArrayList<>();
	private ArrayList<GenericFile> mixedByChrBimFile = new ArrayList<>();
	private ArrayList<GenericFile> mixedByChrFamFile = new ArrayList<>();
	private ArrayList<GenericFile> mixedBedToBedLogFile = new ArrayList<>();

	private ArrayList<GenericFile> mixedSplitChr23MalesLogFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedSplitChr23FemalesLogFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedChr23BedMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedChr23BimMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedChr23FamMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedChr23BedFemalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedChr23BimFemalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedChr23FamFemalesFile = new ArrayList<GenericFile>();

	private String mixedBedDir = null;

	private ArrayList<GenericFile> mixedGenFile = new ArrayList<>();
	private ArrayList<GenericFile> mixedPairsFile = new ArrayList<>();
	private ArrayList<GenericFile> mixedSampleFile = new ArrayList<>();

	private ArrayList<String> mixedTypeSample = new ArrayList<>();
	private ArrayList<GenericFile> mixedGmapFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedPhasingHapsFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingHapsMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingHapsFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedPhasingSampleFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingNewSampleFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingSampleMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingSampleFemalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingNewSampleMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingNewSampleFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedPhasingLogFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingLogMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedPhasingLogFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedExcludedSnpsFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedExcludedSnpsMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedExcludedSnpsFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedFilteredHaplotypesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedFilteredHaplotypesSampleFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesSampleMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesSampleFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedFilteredHaplotypesLogFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesLogMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesLogFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFemalesFile = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFileBgzip = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfMalesFileBgzip = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFemalesFileBgzip = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFileBgzipIndexed = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfMalesFileBgzipIndexed = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFemalesFileBgzipIndexed = new ArrayList<GenericFile>();

	private ArrayList<GenericFile> mixedListOfSnpsFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedListOfSnpsMalesFile = new ArrayList<GenericFile>();
	private ArrayList<GenericFile> mixedListOfSnpsFemalesFile = new ArrayList<GenericFile>();

	/**
	 * A constructor for the class
	 * 
	 * @param parsingArgs
	 * @param myOutDir
	 */
	public CommonFiles(ParseCmdLine parsingArgs, String myOutDir) {
		String inputFormat = parsingArgs.getInputFormat();
		String cohort = parsingArgs.getCohort();

		// We create the names for mixed:
		String tmpOutDir = myOutDir + File.separator + cohort + File.separator + "common" + File.separator + "mixed";
		startChr = parsingArgs.getStart();
		endChr = parsingArgs.getEnd();
		
		LOGGER.info("Creating common files for chromosomes from  " + startChr + " to " + endChr);

		/* We create the input bed file names for the mixed */
		if (inputFormat.equals("BED")) {
			mixedBedDir = parsingArgs.getBedDir();
			String tmpFileName = parsingArgs.getBedFileName();
			mixedBedFile = new GenericFile(mixedBedDir, tmpFileName, "uncompressed", "none");
			// mixedBedFile.add(myMixedBedFile);

			tmpFileName = parsingArgs.getBimFileName();
			mixedBimFile = new GenericFile(mixedBedDir, tmpFileName, "uncompressed", "none");
			// mixedBimFile.add(myMixedBimFile);

			tmpFileName = parsingArgs.getFamFileName();
			mixedFamFile = new GenericFile(mixedBedDir, tmpFileName, "uncompressed", "none");
			// mixedFamFile.add(myMixedFamFile);
		}

		for (int i = startChr; i <= endChr; i++) {
			int chromo = i;

			// First: We create the output directory for this chromosome
			String theOutputDir = tmpOutDir + File.separator + "Chr_" + chromo;
			outputMixedDir.add(theOutputDir);

			String aTmpDir = null;
			String aTmpFileName = null;
			String aTmpSampleFileName = null;
			// String aTmpSampleFile = null;

			/* First, we ask for the inputFormat (GEN or BED) */
			if (inputFormat.equals("BED")) {
				/* We create the bed files name */
				String tmpBedFile = "mixed_" + cohort + "_chr_" + chromo + ".bed";
				GenericFile myMixedByChrBedFile = new GenericFile(theOutputDir, tmpBedFile, "decompressed", "none");
				mixedByChrBedFile.add(myMixedByChrBedFile);

				String tmpFile = "mixed_" + cohort + "_chr_" + chromo + ".bim";
				GenericFile myMixedByChrBimFile = new GenericFile(theOutputDir, tmpFile, "decompressed", "none");
				mixedByChrBimFile.add(myMixedByChrBimFile);

				tmpFile = "mixed_" + cohort + "_chr_" + chromo + ".fam";
				GenericFile myMixedByChrFamFile = new GenericFile(theOutputDir, tmpFile, "decompressed", "none");
				mixedByChrFamFile.add(myMixedByChrFamFile);

				tmpFile = "mixed_" + cohort + "_chr_" + chromo + ".log";
				GenericFile myMixedBedToBedLogFile = new GenericFile(theOutputDir, tmpFile, "decompressed", "none");
				mixedBedToBedLogFile.add(myMixedBedToBedLogFile);

				if (chromo == 23) {
					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_males.bed";
					GenericFile myMixedChr23BedMalesFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedChr23BedMalesFile.add(myMixedChr23BedMalesFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_males.bim";
					GenericFile myMixedChr23BimMalesFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedChr23BimMalesFile.add(myMixedChr23BimMalesFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_males.fam";
					GenericFile myMixedChr23FamMalesFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedChr23FamMalesFile.add(myMixedChr23FamMalesFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_females.bed";
					GenericFile myMixedChr23BedFemalesFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedChr23BedFemalesFile.add(myMixedChr23BedFemalesFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_females.bim";
					GenericFile myMixedChr23BimFemalesFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedChr23BimFemalesFile.add(myMixedChr23BimFemalesFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_females.fam";
					GenericFile myMixedChr23FamFemalesFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedChr23FamFemalesFile.add(myMixedChr23FamFemalesFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_males.log";
					GenericFile myMixedSplitChr23MalesLogFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedSplitChr23MalesLogFile.add(myMixedSplitChr23MalesLogFile);

					tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_females.log";
					GenericFile myMixedSplitChr23FemalesLogFile = new GenericFile(theOutputDir, tmpFile, "decompressed",
							"none");
					mixedSplitChr23FemalesLogFile.add(myMixedSplitChr23FemalesLogFile);
				}

			} else if (inputFormat.equals("GEN")) {
				/* We create the input gen file name */
				// mixedGenDir.add(parsingArgs.getChrDir("mixed"));
				aTmpDir = parsingArgs.getChrDir();
				aTmpFileName = parsingArgs.getGenFileName(chromo);
				GenericFile myMixedGenFile = new GenericFile(aTmpDir, aTmpFileName, "decompressed", "none");
				mixedGenFile.add(myMixedGenFile);
			} else {
				System.err.println(
						"[CommonFiles] Error, this type of input format: " + inputFormat + " does not exist!.");
				System.exit(1);
			}
			aTmpSampleFileName = parsingArgs.getSampleFileName();
			GenericFile myMixedSampleFile = new GenericFile(parsingArgs.getSampleDir(), aTmpSampleFileName,
					"decompressed", "none");
			mixedSampleFile.add(myMixedSampleFile);

			// We create the output pairs file name for mixed
			GenericFile myMixedPairsFile = new GenericFile(theOutputDir, "mixed_chr_" + chromo + ".pairs",
					"decompressed", "none");
			
			mixedPairsFile.add(myMixedPairsFile);

			GenericFile myMixedPhasingHapsFile = new GenericFile(theOutputDir,
					"mixed_phasing_chr_" + chromo + ".haps.gz", "compressed", "none");
			mixedPhasingHapsFile.add(myMixedPhasingHapsFile);

			GenericFile myMixedPhasingSampleFile = new GenericFile(theOutputDir,
					"mixed_phasing_chr_" + chromo + ".sample", "decompressed", "none");
			mixedPhasingSampleFile.add(myMixedPhasingSampleFile);
			
			GenericFile myMixedPhasingNewSampleFile = new GenericFile(theOutputDir,
					"new_mixed_phasing_chr_" + chromo + ".sample", "decompressed", "none");
			mixedPhasingNewSampleFile.add(myMixedPhasingNewSampleFile);

			GenericFile myMixedPhasingLogFile = new GenericFile(theOutputDir, "mixed_phasing_chr_" + chromo + ".log",
					"decompressed", "none");
			mixedPhasingLogFile.add(myMixedPhasingLogFile);

			GenericFile myMixedExcludedSnpsFile = new GenericFile(theOutputDir,
					"mixed_excluded_snps_chr_" + chromo + ".txt", "decompressed", "none");
			mixedExcludedSnpsFile.add(myMixedExcludedSnpsFile);

			GenericFile myMixedFilteredHapsFile = new GenericFile(theOutputDir,
					"mixed_phasing_filtered_chr_" + chromo + ".haps.gz", "compressed", "none");
			mixedFilteredHaplotypesFile.add(myMixedFilteredHapsFile);

			GenericFile myMixedFilteredHapsSampleFile = new GenericFile(theOutputDir,
					"mixed_phasing_filtered_chr_" + chromo + "_" + aTmpSampleFileName, "decompressed", "none");
			mixedFilteredHaplotypesSampleFile.add(myMixedFilteredHapsSampleFile);

			GenericFile myMixedFilteredHapsLogFile = new GenericFile(theOutputDir,
					"mixed_phasing_filtered_chr_" + chromo + ".log", "decompressed", "none");
			mixedFilteredHaplotypesLogFile.add(myMixedFilteredHapsLogFile);

			GenericFile myMixedFilteredHapsVcfFile = new GenericFile(theOutputDir,
					"mixed_phasing_filtered_chr_" + chromo + ".vcf", "decompressed", "none");
			mixedFilteredHaplotypesVcfFile.add(myMixedFilteredHapsVcfFile);

			GenericFile myMixedFilteredHapsVcfFileBgzip = new GenericFile(theOutputDir,
					"mixed_phasing_filtered_chr_" + chromo + ".vcf.gz", "decompressed", "none");
			mixedFilteredHaplotypesVcfFileBgzip.add(myMixedFilteredHapsVcfFileBgzip);

			GenericFile myMixedFilteredHapsVcfFileBgzipIndexed = new GenericFile(theOutputDir,
					"mixed_phasing_filtered_chr_" + chromo + ".vcf.gz.tbi", "decompressed", "none");
			mixedFilteredHaplotypesVcfFileBgzipIndexed.add(myMixedFilteredHapsVcfFileBgzipIndexed);

			GenericFile myMixedListOfSnpsFile = new GenericFile(theOutputDir,
					"mixed_filtered_snps_list_chr_" + chromo + ".txt", "decompressed", "none");
			mixedListOfSnpsFile.add(myMixedListOfSnpsFile);

			GenericFile myMixedGmapFile = new GenericFile(theOutputDir, "mixed_genetic_map_chr_" + chromo + ".txt",
					"decompressed", "none");
			mixedGmapFile.add(myMixedGmapFile);

			if (chromo == 23) {

				GenericFile myMixedPhasingHapsMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_chr_" + chromo + "_males.haps.gz", "compressed", "none");
				mixedPhasingHapsMalesFile.add(myMixedPhasingHapsMalesFile);

				GenericFile myMixedPhasingSampleMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_chr_" + chromo + "_males.sample", "decompressed", "none");
				mixedPhasingSampleMalesFile.add(myMixedPhasingSampleMalesFile);
				
				GenericFile myMixedPhasingNewSampleMalesFile = new GenericFile(theOutputDir,
						"new_mixed_phasing_chr_" + chromo + "_males.sample", "decompressed", "none");
				mixedPhasingNewSampleMalesFile.add(myMixedPhasingNewSampleMalesFile);

				GenericFile myMixedPhasingLogMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_chr_" + chromo + "_males.log", "decompressed", "none");
				mixedPhasingLogMalesFile.add(myMixedPhasingLogMalesFile);

				GenericFile myMixedExcludedSnpsMalesFile = new GenericFile(theOutputDir,
						"mixed_excluded_snps_chr_" + chromo + "_males.txt", "decompressed", "none");
				mixedExcludedSnpsMalesFile.add(myMixedExcludedSnpsMalesFile);

				GenericFile myMixedFilteredHapsMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_males.haps.gz", "compressed", "none");
				mixedFilteredHaplotypesMalesFile.add(myMixedFilteredHapsMalesFile);

				GenericFile myMixedFilteredHapsSampleMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_males_" + aTmpSampleFileName, "decompressed",
						"none");
				mixedFilteredHaplotypesSampleMalesFile.add(myMixedFilteredHapsSampleMalesFile);

				GenericFile myMixedFilteredHapsLogMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_males.log", "decompressed", "none");
				mixedFilteredHaplotypesLogMalesFile.add(myMixedFilteredHapsLogMalesFile);

				GenericFile myMixedFilteredHapsVcfMalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_males.vcf", "decompressed", "none");
				mixedFilteredHaplotypesVcfMalesFile.add(myMixedFilteredHapsVcfMalesFile);

				GenericFile myMixedFilteredHapsVcfMalesFileBgzip = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_males.vcf.gz", "decompressed", "none");
				mixedFilteredHaplotypesVcfMalesFileBgzip.add(myMixedFilteredHapsVcfMalesFileBgzip);

				GenericFile myMixedFilteredHapsVcfMalesFileBgzipIndexed = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_males.vcf.gz.tbi", "decompressed", "none");
				mixedFilteredHaplotypesVcfMalesFileBgzipIndexed.add(myMixedFilteredHapsVcfMalesFileBgzipIndexed);

				GenericFile myMixedListOfSnpsMalesFile = new GenericFile(theOutputDir,
						"mixed_filtered_snps_list_chr_" + chromo + "_males.txt", "decompressed", "none");
				mixedListOfSnpsMalesFile.add(myMixedListOfSnpsMalesFile);

				GenericFile myMixedPhasingHapsFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_chr_" + chromo + "_females.haps.gz", "compressed", "none");
				mixedPhasingHapsFemalesFile.add(myMixedPhasingHapsFemalesFile);

				GenericFile myMixedPhasingSampleFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_chr_" + chromo + "_females.sample", "decompressed", "none");
				mixedPhasingSampleFemalesFile.add(myMixedPhasingSampleFemalesFile);
				
				GenericFile myMixedPhasingNewSampleFemalesFile = new GenericFile(theOutputDir,
						"new_mixed_phasing_chr_" + chromo + "_females.sample", "decompressed", "none");
				mixedPhasingNewSampleFemalesFile.add(myMixedPhasingNewSampleFemalesFile);

				GenericFile myMixedPhasingLogFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_chr_" + chromo + "_females.log", "decompressed", "none");
				mixedPhasingLogFemalesFile.add(myMixedPhasingLogFemalesFile);

				GenericFile myMixedExcludedSnpsFemalesFile = new GenericFile(theOutputDir,
						"mixed_excluded_snps_chr_" + chromo + "_females.txt", "decompressed", "none");
				mixedExcludedSnpsFemalesFile.add(myMixedExcludedSnpsFemalesFile);

				GenericFile myMixedFilteredHapsFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_females.haps.gz", "compressed", "none");
				mixedFilteredHaplotypesFemalesFile.add(myMixedFilteredHapsFemalesFile);

				GenericFile myMixedFilteredHapsSampleFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_females_" + aTmpSampleFileName, "decompressed",
						"none");
				mixedFilteredHaplotypesSampleFemalesFile.add(myMixedFilteredHapsSampleFemalesFile);

				GenericFile myMixedFilteredHapsLogFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_females.log", "decompressed", "none");
				mixedFilteredHaplotypesLogFemalesFile.add(myMixedFilteredHapsLogFemalesFile);

				GenericFile myMixedFilteredHapsVcfFemalesFile = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_females.vcf", "decompressed", "none");
				mixedFilteredHaplotypesVcfFemalesFile.add(myMixedFilteredHapsVcfFemalesFile);

				GenericFile myMixedFilteredHapsVcfFemalesFileBgzip = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_females.vcf.gz", "decompressed", "none");
				mixedFilteredHaplotypesVcfFemalesFileBgzip.add(myMixedFilteredHapsVcfFemalesFileBgzip);

				GenericFile myMixedFilteredHapsVcfFemalesFileBgzipIndexed = new GenericFile(theOutputDir,
						"mixed_phasing_filtered_chr_" + chromo + "_females.vcf.gz.tbi", "decompressed", "none");
				mixedFilteredHaplotypesVcfFemalesFileBgzipIndexed.add(myMixedFilteredHapsVcfFemalesFileBgzipIndexed);

				GenericFile myMixedListOfSnpsFemalesFile = new GenericFile(theOutputDir,
						"mixed_filtered_snps_list_chr_" + chromo + "_females.txt", "decompressed", "none");
				mixedListOfSnpsFemalesFile.add(myMixedListOfSnpsFemalesFile);

			}

		}

	}

	/**
	 * Method to access bedFile information
	 * 
	 * @return
	 */
	public String getBedFile() {
		return mixedBedFile.getFullName();
	}

	/**
	 * Method to access byChrBedFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getByChrBedFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedByChrBedFile.get(index).getFullName();
	}

	/**
	 * Method to access bimFile information
	 * 
	 * @return
	 */
	public String getBimFile() {
		return mixedBimFile.getFullName();
	}

	/**
	 * Method to access byChrBimFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getByChrBimFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedByChrBimFile.get(index).getFullName();
	}

	/**
	 * Method to access famFile information
	 * 
	 * @return
	 */
	public String getFamFile() {
		return mixedFamFile.getFullName();
	}

	/**
	 * Method to access byChrFamFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getByChrFamFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedByChrFamFile.get(index).getFullName();
	}

	/**
	 * Method to access BedToBedLogFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getBedToBedLogFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedBedToBedLogFile.get(index).getFullName();
	}

	/**
	 * Method to access chr 23 Bed males file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getChr23BedMalesFile() {
		return mixedChr23BedMalesFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 Bim males file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getChr23BimMalesFile() {
		return mixedChr23BimMalesFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 Fam males file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getChr23FamMalesFile() {
		return mixedChr23FamMalesFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 Bed females file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getChr23BedFemalesFile() {
		return mixedChr23BedFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 Bim females file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getChr23BimFemalesFile() {
		return mixedChr23BimFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 Fam females file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getChr23FamFemalesFile() {
		return mixedChr23FamFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 split males log file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getSplitChr23MalesLogFile() {
		return mixedSplitChr23MalesLogFile.get(0).getFullName();
	}

	/**
	 * Method to access chr 23 split females log file information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getSplitChr23FemalesLogFile() {
		return mixedSplitChr23FemalesLogFile.get(0).getFullName();
	}

	/**
	 * Method to access genDir information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getGenDir(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedGenFile.get(index).getDir();
	}

	/**
	 * Method to access genFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getGenFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedGenFile.get(index).getFullName();
	}

	/**
	 * Method to access the final status of a genFile
	 * 
	 * @param chromo
	 * @return
	 */
	public String getGenFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedGenFile.get(index).getFinalStatus();
	}

	/**
	 * Method to access pairsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPairsFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		System.out.println("Serving " + mixedPairsFile.get(index).getFullName());
		return mixedPairsFile.get(index).getFullName();
	}

	/**
	 * Method to set the finalStatus of pairsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setPairsFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedPairsFile.get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access pairsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPairsFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPairsFile.get(index).getFinalStatus();
	}

	/**
	 * Method to access sampleDir information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getSampleDir(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedSampleFile.get(index).getDir();
	}

	/**
	 * Method to access sampleFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getSampleFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedSampleFile.get(index).getFullName();
	}

	/**
	 * Method to access the final status of a sampleFile
	 * 
	 * @param chromo
	 * @return
	 */
	public String getSampleFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedSampleFile.get(index).getFinalStatus();
	}

	/**
	 * Method to access typeSample information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getTypeSample(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedTypeSample.get(index);
	}

	/**
	 * Method to access shapeitHapsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingHapsFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingHapsFile.get(index).getFullName();
	}

	/**
	 * Method to set finalStatus of shapeitHapsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setPhasingHapsFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedPhasingHapsFile.get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set finalStatus of shapeit Sample file
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setPhasingSampleFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedPhasingSampleFile.get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to access the final status information of shapeitHapsFile
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingHapsFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingHapsFile.get(index).getFinalStatus();
	}

	/**
	 * Method to access shapeitSampleFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingSampleFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingSampleFile.get(index).getFullName();
	}
	
	/**
	 * Method to access shapeitSampleFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingNewSampleFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingNewSampleFile.get(index).getFullName();
	}

	/**
	 * Method to access the final status information of a shapeitSampleFile
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingSampleFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingSampleFile.get(index).getFinalStatus();
	}

	/**
	 * Method to access shapeitLogFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingLogFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingLogFile.get(index).getFullName();
	}

	/**
	 * Method to access the final status information of a shapeitLogFile
	 * 
	 * @param chromo
	 * @return
	 */
	public String getPhasingLogFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedPhasingLogFile.get(index).getFinalStatus();
	}

	// Method to access getPhasingHapsMalesFile
	public String getPhasingHapsMalesFile() {
		return mixedPhasingHapsMalesFile.get(0).getFullName();
	}

	// Method to access getPhasingSampleMalesFile
	public String getPhasingSampleMalesFile() {
		return mixedPhasingSampleMalesFile.get(0).getFullName();
	}
	
	// Method to access getPhasingSampleMalesFile
	public String getPhasingNewSampleMalesFile() {
		return mixedPhasingNewSampleMalesFile.get(0).getFullName();
	}

	// Method to access getPhasingLogMalesFile
	public String getPhasingLogMalesFile() {
		return mixedPhasingLogMalesFile.get(0).getFullName();
	}

	// Method to access getPhasingHapsFemalesFile
	public String getPhasingHapsFemalesFile() {
		return mixedPhasingHapsFemalesFile.get(0).getFullName();
	}

	// Method to access getPhasingSampleFemalesFile
	public String getPhasingSampleFemalesFile() {
		return mixedPhasingSampleFemalesFile.get(0).getFullName();
	}
	
	// Method to access getPhasingSampleFemalesFile
	public String getPhasingNewSampleFemalesFile() {
		return mixedPhasingNewSampleFemalesFile.get(0).getFullName();
	}

	// Method to access getPhasingLogFemalesFile
	public String getPhasingLogFemalesFile() {
		return mixedPhasingLogFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access excludedSnpsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getExcludedSnpsFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedExcludedSnpsFile.get(index).getFullName();
	}

	/**
	 * Method to access the final status information of a excludedSnpsFile
	 * 
	 * @param chromo
	 * @return
	 */
	public String getExcludedSnpsFileFinalStatus(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedExcludedSnpsFile.get(index).getFinalStatus();
	}

	// Method to access getExcludedSnpsFileMales
	public String getExcludedSnpsMalesFile() {
		return mixedExcludedSnpsMalesFile.get(0).getFullName();
	}

	// Method to access getExcludedSnpsFemalesFile
	public String getExcludedSnpsFemalesFile() {
		return mixedExcludedSnpsFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access shapeitHapsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getFilteredHaplotypesFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedFilteredHaplotypesFile.get(index).getFullName();
	}

	/**
	 * Method to set finalStatus of shapeitHapsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setFilteredHaplotypesFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedFilteredHaplotypesFile.get(index).setFinalStatus(finalStatus);
	}

	// Method to access filteredHaplotypesMalesFile
	public String getFilteredHaplotypesMalesFile() {
		return mixedFilteredHaplotypesMalesFile.get(0).getFullName();
	}

	// Method to access filteredHaplotypesFemalesFile
	public String getFilteredHaplotypesFemalesFile() {
		return mixedFilteredHaplotypesFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access shapeitHapsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getFilteredHaplotypesSampleFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedFilteredHaplotypesSampleFile.get(index).getFullName();
	}

	/**
	 * Method to set finalStatus of shapeitHapsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setFilteredHaplotypesSampleFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedFilteredHaplotypesSampleFile.get(index).setFinalStatus(finalStatus);
	}

	// Method to access filteredHaplotypesSampleMalesFile
	public String getFilteredHaplotypesSampleMalesFile() {
		return mixedFilteredHaplotypesSampleMalesFile.get(0).getFullName();
	}

	// Method to access filteredHaplotypesSampleFemalesFile
	public String getFilteredHaplotypesSampleFemalesFile() {
		return mixedFilteredHaplotypesSampleFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access filteredHaplotypesLogFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getFilteredHaplotypesLogFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedFilteredHaplotypesLogFile.get(index).getFullName();
	}

	/**
	 * Method to set finalStatus of shapeitHapsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setFilteredHaplotypesLogFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedFilteredHaplotypesLogFile.get(index).setFinalStatus(finalStatus);
	}

	// Method to access filteredHaplotypesLogMalesFile
	public String getFilteredHaplotypesLogMalesFile() {
		return mixedFilteredHaplotypesLogMalesFile.get(0).getFullName();
	}

	// Method to access filteredHaplotypesLogFemalesFile
	public String getFilteredHaplotypesLogFemalesFile() {
		return mixedFilteredHaplotypesLogFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to access filteredHaplotypesVcfFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getFilteredHaplotypesVcfFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedFilteredHaplotypesVcfFile.get(index).getFullName();
	}

	// Method to access filteredHaplotypesVcfMalesFile
	public String getFilteredHaplotypesVcfMalesFile() {
		return mixedFilteredHaplotypesVcfMalesFile.get(0).getFullName();
	}

	// Method to access filteredHaplotypesVcfFemalesFile
	public String getFilteredHaplotypesVcfFemalesFile() {
		return mixedFilteredHaplotypesVcfFemalesFile.get(0).getFullName();
	}

	/**
	 * Method to set finalStatus of setFilteredHaplotypesVcfFileFinalStatus
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setFilteredHaplotypesVcfFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedFilteredHaplotypesVcfFile.get(index).setFinalStatus(finalStatus);
	}

	// Method to access filteredHaplotypesVcfFile information
	public String getFilteredHaplotypesVcfFileBgzip(int chromo) {
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedFilteredHaplotypesVcfFileBgzip.get(index).getFullName();
	}

	// Method to set finalStatus of setFilteredHaplotypesVcfFileFinalStatus
	public void setFilteredHaplotypesVcfFileBgzipFinalStatus(int chromo, String finalStatus) {
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedFilteredHaplotypesVcfFileBgzip.get(index).setFinalStatus(finalStatus);
	}

	// Method to access filteredHaplotypesVcfMalesBgzipFile
	public String getFilteredHaplotypesVcfMalesFileBgzip() {
		return mixedFilteredHaplotypesVcfMalesFileBgzip.get(0).getFullName();
	}

	// Method to access filteredHaplotypesVcfFemalesFileBgzip
	public String getFilteredHaplotypesVcfFemalesFileBgzip() {
		return mixedFilteredHaplotypesVcfFemalesFileBgzip.get(0).getFullName();
	}

	// Method to access filteredHaplotypesVcfFile information
	public String getFilteredHaplotypesVcfFileBgzipIndexed(int chromo) {
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedFilteredHaplotypesVcfFileBgzipIndexed.get(index).getFullName();
	}

	// Method to set finalStatus of setFilteredHaplotypesVcfFileFinalStatus
	public void setFilteredHaplotypesVcfFileBgzipIndexedFinalStatus(int chromo, String finalStatus) {
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedFilteredHaplotypesVcfFileBgzipIndexed.get(index).setFinalStatus(finalStatus);
	}

	// Method to access filteredHaplotypesVcfMalesBgzipFileIndexed
	public String getFilteredHaplotypesVcfMalesFileBgzipIndexed() {
		return mixedFilteredHaplotypesVcfMalesFileBgzipIndexed.get(0).getFullName();
	}

	// Method to access filteredHaplotypesVcfFemalesFileBgzipIndexed
	public String getFilteredHaplotypesVcfFemalesFileBgzipIndexed() {
		return mixedFilteredHaplotypesVcfFemalesFileBgzipIndexed.get(0).getFullName();
	}

	/**
	 * Method to access listOfSnpsFile information
	 * 
	 * @param chromo
	 * @return
	 */
	public String getListOfSnpsFile(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedListOfSnpsFile.get(index).getFullName();
	}

	/**
	 * Method to set finalStatus of listOfSnpsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setListOfSnpsFileFinalStatus(int chromo, String finalStatus) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedListOfSnpsFile.get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to set finalStatus of listOfSnpsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public String getGmapFile(int chromo) {
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		return mixedGmapFile.get(index).getFullName();
	}

	/**
	 * Method to set finalStatus of listOfSnpsFile
	 * 
	 * @param chromo
	 * @param finalStatus
	 */
	public void setGmapFileFinalStatus(int chromo, String finalStatus) {
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		mixedGmapFile.get(index).setFinalStatus(finalStatus);
	}

	/**
	 * Method to print files information
	 * 
	 * @param chromo
	 */
	public void printCommonFiles(int chromo) {
		// Check that chromo index is within the bounds
		checkChromoIndex(chromo);

		int index = chromo - startChr;
		LOGGER.info("-------------------------------------------------");
		LOGGER.info("Mixed files information for the chromosome " + chromo + "(" + index + ")");
		LOGGER.info("mixedGenFile         : " + mixedGenFile.get(index).getFullName());

		LOGGER.info("mixedPairsFile       : " + mixedPairsFile.get(index).getFullName());
		LOGGER.info("mixedSampleFile      : " + mixedSampleFile.get(index));
		LOGGER.info("mixedTypeSample      : " + mixedTypeSample.get(index));

		LOGGER.info("mixedPhasingHapsFile   : " + mixedPhasingHapsFile.get(index).getFullName());
		LOGGER.info("mixedPhasingSampleFile : " + mixedPhasingSampleFile.get(index).getFullName());
		LOGGER.info("mixedPhasingLogFile    : " + mixedPhasingLogFile.get(index).getFullName());
		LOGGER.info("mixedExcludedSnpsFile :  " + mixedExcludedSnpsFile.get(index).getFullName());

		LOGGER.info("-------------------------------------------------");
	}

	private void checkChromoIndex(int chromo) {
		if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
			LOGGER.fatal("[CommonFiles] Error, chromosome " + chromo + "does not exist");
			System.exit(1);
		}
	}
}
