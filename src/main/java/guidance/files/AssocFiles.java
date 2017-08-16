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

import guidance.utils.ChromoInfo;
import guidance.utils.ParseCmdLine;

import java.io.File;
import java.util.ArrayList;


public class AssocFiles {

    private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;
    private static final int MIN_LIMIT = 1;
    private static final int MAX_LIMIT = 252_000_000;

    private ArrayList<ArrayList<ArrayList<String>>> outDir = new ArrayList<>();

    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestOutFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestLogFileName = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> snptestLogFile = new ArrayList<>();

    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryFilteredFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> summaryCondensedFile = new ArrayList<>();

    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedFilteredFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> combinedCondensedFile = new ArrayList<>();

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
    public AssocFiles(ParseCmdLine parsingArgs, ChromoInfo generalChromoInfo, String baseOutDir, List<String> refPanels) {
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
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestOutFile = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSnptestLogFile = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFile = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryFilteredFile = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListSummaryCondensedFile = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedFilteredFile = new ArrayList<>();
            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelListCombinedCondensedFile = new ArrayList<>();

            String rPanel = null;
            String testTypeOutDir2 = baseOutDir + File.separator + "associations" + File.separator + testTypeName + File.separator
                    + mixedCohort + "_combined_panels";
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
                ArrayList<ArrayList<GenericFile>> chromoListCombinedCondensedFile = new ArrayList<>();

                for (int i = this.startChr; i <= this.endChr; i++) {
                    int chromo = i;
                    int maxSize = generalChromoInfo.getMaxSize(chromo);
                    int total_chunks = maxSize / chunkSize;
                    int module = maxSize % chunkSize;
                    if (module != 0) {
                        total_chunks++;
                    }
                    int lim1 = 1;
                    int lim2 = lim1 + chunkSize - 1;

                    ArrayList<GenericFile> chunkListCombinedFilteredFile = new ArrayList<>();
                    ArrayList<GenericFile> chunkListCombinedCondensedFile = new ArrayList<>();

                    for (int k = 0; k < total_chunks; k++) {
                        String tmpCombinedFilteredFileName = prefixFilteredName + "_chr_" + chromo + "_" + lim1 + "_" + lim2
                                + "_combined.txt.gz";
                        GenericFile myChunkListCombinedFilteredFile = new GenericFile(testTypeOutDir2, tmpCombinedFilteredFileName,
                                "uncompressed", "none");
                        chunkListCombinedFilteredFile.add(myChunkListCombinedFilteredFile);

                        String tmpCombinedCondensedFileName = prefixCondensedName + "_chr_" + chromo + "_" + lim1 + "_" + lim2
                                + "_combined.txt.gz";
                        GenericFile myChunkListCombinedCondensedFile = new GenericFile(testTypeOutDir2, tmpCombinedCondensedFileName,
                                "uncompressed", "none");
                        chunkListCombinedCondensedFile.add(myChunkListCombinedCondensedFile);

                        lim1 = lim1 + chunkSize;
                        lim2 = lim2 + chunkSize;
                    }

                    chromoListCombinedFilteredFile.add(chunkListCombinedFilteredFile);
                    chromoListCombinedCondensedFile.add(chunkListCombinedCondensedFile);

                } // End for Chromo

                rpanelListCombinedFilteredFile.add(chromoListCombinedFilteredFile);
                rpanelListCombinedCondensedFile.add(chromoListCombinedCondensedFile);

            } // End of for refPanels

            this.combinedFilteredFile.add(rpanelListCombinedFilteredFile);
            this.combinedCondensedFile.add(rpanelListCombinedCondensedFile);

            for (int j = 0; j < refPanels.size(); j++) {
                // String rpanel = refPanels.get(j);
                rPanel = refPanels.get(j);
                String rpanelOutDir = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel;

                ArrayList<String> chromoListOutDir = new ArrayList<>();
                ArrayList<ArrayList<GenericFile>> chromoListSnptestOutFile = new ArrayList<>();
                ArrayList<ArrayList<GenericFile>> chromoListSnptestLogFile = new ArrayList<>();
                ArrayList<ArrayList<GenericFile>> chromoListSummaryFile = new ArrayList<>();
                ArrayList<ArrayList<GenericFile>> chromoListSummaryFilteredFile = new ArrayList<>();
                ArrayList<ArrayList<GenericFile>> chromoListSummaryCondensedFile = new ArrayList<>();

                for (int i = this.startChr; i <= this.endChr; i++) {
                    int chromo = i;

                    String tmpChrDir = rpanelOutDir + File.separator + "Chr_" + chromo;
                    chromoListOutDir.add(tmpChrDir);

                    int maxSize = generalChromoInfo.getMaxSize(chromo);
                    int total_chunks = maxSize / chunkSize;
                    int module = maxSize % chunkSize;
                    if (module != 0)
                        total_chunks++;
                    int lim1 = 1;
                    int lim2 = lim1 + chunkSize - 1;

                    ArrayList<GenericFile> chunkListSnptestOutFile = new ArrayList<>();
                    ArrayList<GenericFile> chunkListSnptestLogFile = new ArrayList<>();
                    ArrayList<GenericFile> chunkListSummaryFile = new ArrayList<>();

                    ArrayList<GenericFile> chunkListSummaryFilteredFile = new ArrayList<>();
                    ArrayList<GenericFile> chunkListSummaryCondensedFile = new ArrayList<>();
                    for (int k = 0; k < total_chunks; k++) {
                        // Now we have to create the impute files for this iteration
                        String tmpSnptestOutFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_" + lim1 + "_" + lim2
                                + "_snptest.out.gz";
                        GenericFile myChunkListSnptestOutFile = new GenericFile(tmpChrDir, tmpSnptestOutFileName, "uncompressed", "none");
                        chunkListSnptestOutFile.add(myChunkListSnptestOutFile);

                        String tmpSnptestLogFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_" + lim1 + "_" + lim2
                                + "_snptest.log";
                        GenericFile myChunkListSnptestLogFile = new GenericFile(tmpChrDir, tmpSnptestLogFileName, "uncompressed", "none");
                        chunkListSnptestLogFile.add(myChunkListSnptestLogFile);

                        String tmpSummaryFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_" + lim1 + "_" + lim2
                                + "_summary.txt.gz";
                        GenericFile myChunkListSummaryFile = new GenericFile(tmpChrDir, tmpSummaryFileName, "uncompressed", "none");
                        chunkListSummaryFile.add(myChunkListSummaryFile);

                        String tmpSummaryFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_" + lim1 + "_" + lim2
                                + "_summary_filtered.txt.gz";
                        GenericFile myChunkListSummaryFilteredFile = new GenericFile(tmpChrDir, tmpSummaryFilteredFileName, "uncompressed",
                                "none");
                        chunkListSummaryFilteredFile.add(myChunkListSummaryFilteredFile);

                        String tmpSummaryCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_" + lim1 + "_" + lim2
                                + "_summary_condensed.txt.gz";
                        GenericFile myChunkListSummaryCondensedFile = new GenericFile(tmpChrDir, tmpSummaryCondensedFileName,
                                "uncompressed", "none");
                        chunkListSummaryCondensedFile.add(myChunkListSummaryCondensedFile);

                        lim1 = lim1 + chunkSize;
                        lim2 = lim2 + chunkSize;
                    }

                    chromoListSnptestOutFile.add(chunkListSnptestOutFile);
                    chromoListSnptestLogFile.add(chunkListSnptestLogFile);
                    chromoListSummaryFile.add(chunkListSummaryFile);

                    chromoListSummaryFilteredFile.add(chunkListSummaryFilteredFile);
                    chromoListSummaryCondensedFile.add(chunkListSummaryCondensedFile);

                } // End for chromo

                rpanelListOutDir.add(chromoListOutDir);

                rpanelListSnptestOutFile.add(chromoListSnptestOutFile);
                rpanelListSnptestLogFile.add(chromoListSnptestLogFile);
                rpanelListSummaryFile.add(chromoListSummaryFile);

                rpanelListSummaryFilteredFile.add(chromoListSummaryFilteredFile);
                rpanelListSummaryCondensedFile.add(chromoListSummaryCondensedFile);

            } // End of for refPanels

            this.snptestOutFile.add(rpanelListSnptestOutFile);
            this.snptestLogFile.add(rpanelListSnptestLogFile);
            this.summaryFile.add(rpanelListSummaryFile);

            this.summaryFilteredFile.add(rpanelListSummaryFilteredFile);
            this.summaryCondensedFile.add(rpanelListSummaryCondensedFile);

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
        // int i= chromo - 1 ;
        return this.outDir.get(testTypeIndex).get(rPanelIndex).get(i);
    }

    /**
     * Method to access snptestOutFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSnptestOutFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
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
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of snptestOutFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setSnptestOutFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
            String finalStatus) {

        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        this.snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of snptestOutFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSnptestOutFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
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
    public String getSummaryFilteredFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
    }

    /**
     * Method to access summaryFilteredFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSummaryFilteredFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of summaryFilteredFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setSummaryFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
            String finalStatus) {

        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        this.summaryFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to get the finalStatus of summaryFilteredFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSummaryFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
    }

    /**
     * Method to access combinedFilteredFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getCombinedFilteredFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.combinedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
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
    public String getCombinedFilteredFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.combinedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of combinedFilteredFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setCombinedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
            String finalStatus) {

        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        this.combinedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to get the finalStatus of combinedFilteredFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getCombinedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.combinedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
    }

    /** Condensed ------------------------------------- */

    /**
     * Method to access summaryCondensedFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSummaryCondensedFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
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
    public String getSummaryCondensedFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of summaryCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setSummaryCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
            String finalStatus) {

        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        this.summaryCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to get the finalStatus of summaryCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSummaryCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
    }

    /**
     * Method to access combinedCondensedFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getCombinedCondensedFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.combinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
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
    public String getCombinedCondensedFile(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.combinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of combinedCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setCombinedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
            String finalStatus) {

        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        this.combinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to get the finalStatus of combinedCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getCombinedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.combinedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
    }

    /**
     * Method to access snptestLogFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSnptestLogFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.snptestLogFileName.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
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
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.snptestLogFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to access snptestLogFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSummaryFileName(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getName();
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
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of snptestLogFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setSummaryFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize,
            String finalStatus) {

        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        this.summaryFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of snptestLogFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getSummaryFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        // The offset is because the array start in position 0
        // int i = chromo -1;
        int i = chromo - this.startChr;
        int index = lim1 / chunkSize;
        return this.summaryFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
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
            System.out.println("-------------------------------------------------");
            System.out.println("Assoc files information for the chromosome " + chromo);
            System.out.println("outDir                  : " + outDir.get(testTypeIndex).get(rPanelIndex).get(i));
            System.out.println("snptestOutFile    : " + snptestOutFile.get(testTypeIndex).get(rPanelIndex).get(i).get(indexLow));
            System.out.println("snptestLogFile    : " + snptestLogFile.get(testTypeIndex).get(rPanelIndex).get(i).get(indexLow));
            System.out.println("-------------------------------------------------");
        }
    }

    private void checkChromoIndex(int chromo) {
        if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
            System.err.println("[AssocFiles] Error, chromosome " + chromo + "does not exist");
            System.exit(1);
        }
    }

    private void checkLimits(int chromo, int lim1, int lim2) {
        if ((lim1 < MIN_LIMIT) || (lim2 > MAX_LIMIT)) {
            System.err.println("[AssocFiles] Error, Chunk " + lim1 + "_" + lim2 + "does not exist for chromosome " + chromo);
            System.exit(1);
        }
    }

}
