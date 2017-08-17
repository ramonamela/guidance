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

    private ArrayList<String> outputMixedDir = new ArrayList<>();

    private GenericFile mixedBedFile = null;
    private GenericFile mixedBimFile = null;
    private GenericFile mixedFamFile = null;

    private ArrayList<GenericFile> mixedByChrBedFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedByChrBimFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedByChrFamFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedBedToBedLogFile = new ArrayList<>();

    private String mixedBedDir = null;

    private ArrayList<GenericFile> mixedGenFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedPairsFile = new ArrayList<>();

    private ArrayList<GenericFile> mixedSampleFile = new ArrayList<>();
    private ArrayList<String> mixedTypeSample = new ArrayList<>();

    private ArrayList<GenericFile> mixedShapeitHapsFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedShapeitSampleFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedShapeitLogFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedExcludedSnpsFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedFilteredHaplotypesFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedFilteredHaplotypesSampleFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedFilteredHaplotypesLogFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedFilteredHaplotypesVcfFile = new ArrayList<>();
    private ArrayList<GenericFile> mixedListOfSnpsFile = new ArrayList<>();


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
        String tmpOutDir = myOutDir + "/" + cohort + "/common/mixed";
        startChr = parsingArgs.getStart();
        endChr = parsingArgs.getEnd();

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
            String theOutputDir = tmpOutDir + "/Chr_" + chromo;
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

                tmpFile = "mixed_" + cohort + "_chr_" + chromo + "_bed2bed.log";
                GenericFile myMixedBedToBedLogFile = new GenericFile(theOutputDir, tmpFile, "decompressed", "none");
                mixedBedToBedLogFile.add(myMixedBedToBedLogFile);

            } else if (inputFormat.equals("GEN")) {
                /* We create the input gen file name */
                // mixedGenDir.add(parsingArgs.getChrDir("mixed"));
                aTmpDir = parsingArgs.getChrDir();
                aTmpFileName = parsingArgs.getGenFileName(chromo);
                GenericFile myMixedGenFile = new GenericFile(aTmpDir, aTmpFileName, "decompressed", "none");
                mixedGenFile.add(myMixedGenFile);
            } else {
                System.err.println("[CommonFiles] Error, this type of input format: " + inputFormat + " does not exist!.");
                System.exit(1);
            }

            aTmpSampleFileName = parsingArgs.getSampleFileName();
            GenericFile myMixedSampleFile = new GenericFile(parsingArgs.getSampleDir(), aTmpSampleFileName, "decompressed", "none");
            mixedSampleFile.add(myMixedSampleFile);

            // We create the output pairs file name for mixed
            GenericFile myMixedPairsFile = new GenericFile(theOutputDir, "mixed_chr_" + chromo + ".pairs", "decompressed", "none");
            mixedPairsFile.add(myMixedPairsFile);

            GenericFile myMixedShapeitHapsFile = new GenericFile(theOutputDir, "mixed_shapeit_chr_" + chromo + ".haps.gz", "compressed",
                    "none");
            mixedShapeitHapsFile.add(myMixedShapeitHapsFile);

            GenericFile myMixedShapeitSampleFile = new GenericFile(theOutputDir, "mixed_shapeit_chr_" + chromo + "_" + aTmpSampleFileName,
                    "decompressed", "none");
            mixedShapeitSampleFile.add(myMixedShapeitSampleFile);

            GenericFile myMixedShapeitLogFile = new GenericFile(theOutputDir, "mixed_shapeit_chr_" + chromo + ".log", "decompressed",
                    "none");
            mixedShapeitLogFile.add(myMixedShapeitLogFile);

            GenericFile myMixedExcludedSnpsFile = new GenericFile(theOutputDir, "mixed_excluded_snps_chr_" + chromo + ".txt",
                    "decompressed", "none");
            mixedExcludedSnpsFile.add(myMixedExcludedSnpsFile);

            GenericFile myMixedFilteredHapsFile = new GenericFile(theOutputDir, "mixed_shapeit_filtered_chr_" + chromo + ".haps.gz",
                    "compressed", "none");
            mixedFilteredHaplotypesFile.add(myMixedFilteredHapsFile);

            GenericFile myMixedFilteredHapsSampleFile = new GenericFile(theOutputDir,
                    "mixed_shapeit_filtered_chr_" + chromo + "_" + aTmpSampleFileName, "decompressed", "none");
            mixedFilteredHaplotypesSampleFile.add(myMixedFilteredHapsSampleFile);

            GenericFile myMixedFilteredHapsLogFile = new GenericFile(theOutputDir, "mixed_shapeit_filtered_chr_" + chromo + ".log",
                    "decompressed", "none");
            mixedFilteredHaplotypesLogFile.add(myMixedFilteredHapsLogFile);

            GenericFile myMixedFilteredHapsVcfFile = new GenericFile(theOutputDir, "mixed_shapeit_filtered_chr_" + chromo + ".vcf.gz",
                    "decompressed", "none");
            mixedFilteredHaplotypesVcfFile.add(myMixedFilteredHapsVcfFile);

            GenericFile myMixedListOfSnpsFile = new GenericFile(theOutputDir, "mixed_filtered_snps_list_chr_" + chromo + ".txt",
                    "decompressed", "none");
            mixedListOfSnpsFile.add(myMixedListOfSnpsFile);
        }
        /*
         * for(int cromo=startChr; cromo<=endChr;cromo++) { printCommonFiles(cromo); }
         */
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
     * Method to access genFileName information
     * 
     * @param chromo
     * @return
     */
    public String getGenFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedGenFile.get(index).getName();
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
     * Method to access pairsFileName information
     * 
     * @param chromo
     * @return
     */
    public String getPairsFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedPairsFile.get(index).getName();
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
     * Method to access sampleFileName information
     * 
     * @param chromo
     * @return
     */
    public String getSampleFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedSampleFile.get(index).getName();
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
     * Method to access shapeitHapsFileName information
     * 
     * @param chromo
     * @return
     */
    public String getShapeitHapsFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitHapsFile.get(index).getName();
    }

    /**
     * Method to access shapeitHapsFile information
     * 
     * @param chromo
     * @return
     */
    public String getShapeitHapsFile(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitHapsFile.get(index).getFullName();
    }

    /**
     * Method to set finalStatus of shapeitHapsFile
     * 
     * @param chromo
     * @param finalStatus
     */
    public void setShapeitHapsFileFinalStatus(int chromo, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        mixedShapeitHapsFile.get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to set finalStatus of shapeit Sample file
     * 
     * @param chromo
     * @param finalStatus
     */
    public void setShapeitSampleFileFinalStatus(int chromo, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        mixedShapeitSampleFile.get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the final status information of shapeitHapsFile
     * 
     * @param chromo
     * @return
     */
    public String getShapeitHapsFileFinalStatus(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitHapsFile.get(index).getFinalStatus();
    }

    /**
     * Method to access shapeitSampleFileName information
     * 
     * @param chromo
     * @return
     */
    public String getShapeitSampleFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitSampleFile.get(index).getName();
    }

    /**
     * Method to access shapeitSampleFile information
     * 
     * @param chromo
     * @return
     */
    public String getShapeitSampleFile(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitSampleFile.get(index).getFullName();
    }

    /**
     * Method to access the final status information of a shapeitSampleFile
     * 
     * @param chromo
     * @return
     */
    public String getShapeitSampleFileFinalStatus(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitSampleFile.get(index).getFinalStatus();
    }

    /**
     * Method to access shapeitLogFileName information
     * 
     * @param chromo
     * @return
     */
    public String getShapeitLogFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitLogFile.get(index).getName();
    }

    /**
     * Method to access shapeitLogFile information
     * 
     * @param chromo
     * @return
     */
    public String getShapeitLogFile(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitLogFile.get(index).getFullName();
    }

    /**
     * Method to access the final status information of a shapeitLogFile
     * 
     * @param chromo
     * @return
     */
    public String getShapeitLogFileFinalStatus(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedShapeitLogFile.get(index).getFinalStatus();
    }

    /**
     * Method to access excludedSnpsFileName information
     * 
     * @param chromo
     * @return
     */
    public String getExcludedSnpsFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedExcludedSnpsFile.get(index).getName();
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

    /**
     * Method to access filteredHaplotypesFileName information
     * 
     * @param chromo
     * @return
     */
    public String getFilteredHaplotypesFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedFilteredHaplotypesFile.get(index).getName();
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

    /**
     * Method to access filteredHaplotypesSampleFileName information
     * 
     * @param chromo
     * @return
     */
    public String getFilteredHaplotypesSampleFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedFilteredHaplotypesSampleFile.get(index).getName();
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

    /**
     * Method to access filteredHaplotypesLogFileName information
     * 
     * @param chromo
     * @return
     */
    public String getFilteredHaplotypesLogFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedFilteredHaplotypesLogFile.get(index).getName();
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

    /**
     * Method to access filteredHaplotypesVcfFileName information
     * 
     * @param chromo
     * @return
     */
    public String getFilteredHaplotypesVcfFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedFilteredHaplotypesVcfFile.get(index).getName();
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

    /**
     * Method to access listOfSnpsFileName information
     * 
     * @param chromo
     * @return
     */
    public String getListOfSnpsFileName(int chromo) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        int index = chromo - startChr;
        return mixedListOfSnpsFile.get(index).getName();
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
        LOGGER.info("mixedSampleFileName  : " + mixedSampleFile.get(index).getName());
        LOGGER.info("mixedSampleFile      : " + mixedSampleFile.get(index));
        LOGGER.info("mixedTypeSample      : " + mixedTypeSample.get(index));

        LOGGER.info("mixedShapeitHapsFile   : " + mixedShapeitHapsFile.get(index).getFullName());
        LOGGER.info("mixedShapeitSampleFile : " + mixedShapeitSampleFile.get(index).getFullName());
        LOGGER.info("mixedShapeitLogFile    : " + mixedShapeitLogFile.get(index).getFullName());
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
