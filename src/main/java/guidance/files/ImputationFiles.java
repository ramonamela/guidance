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
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedInfoFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedSummaryFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedWarningsFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedLogFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFileLogFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> filteredFileRsIdFile = new ArrayList<>();

    private ArrayList<ArrayList<String>> imputedMMOutDir = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMInfoFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMDraftFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMErateFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMRecFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMDoseFile = new ArrayList<>();
    private ArrayList<ArrayList<ArrayList<GenericFile>>> imputedMMLogFile = new ArrayList<>();


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
            System.err.println("[ImputationFiles] Error, this imputation tool (" + imputationTool + ") is not supported yet!.");
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
            ArrayList<String> chromoListImputedOutDir = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedInfoFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedSummaryFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedWarningsFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedLogFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListFilteredFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListFilteredLogFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdFile = new ArrayList<>();

            // int maxSize = chromoInformation.getMaxSize(i);
            for (int i = this.startChr; i <= this.endChr; i++) {
                int chromo = i;
                int lim1 = 1;
                int lim2 = lim1 + chunkSize - 1;
                String tmpChrDir = mixOutDir + File.separator + "Chr_" + chromo;
                chromoListImputedOutDir.add(tmpChrDir);

                ArrayList<GenericFile> chunkListImputedFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedInfoFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedSummaryFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedWarningsFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedLogFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListFilteredFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListFilteredLogFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListFilteredRsIdFile = new ArrayList<>();

                int numberOfChunks = ChromoInfo.getMaxSize(chromo) / chunkSize;
                int module = ChromoInfo.getMaxSize(chromo) % chunkSize;
                if (module != 0)
                    numberOfChunks++;

                for (int k = 0; k < numberOfChunks; k++) {
                    // Now we have to create the impute files for this iteration
                    // String imputedFileName = parsingArgs.getGenFileName(chromo) + "_" + rPanel + "_" + lim1 + "_"
                    // + lim2 +".impute";
                    String imputedFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2 + ".impute";
                    GenericFile myChunkListImputedFile = new GenericFile(tmpChrDir, imputedFileName + ".gz", "compressed", "none");
                    chunkListImputedFile.add(myChunkListImputedFile);

                    GenericFile myChunkListImputedInfoFile = new GenericFile(tmpChrDir, imputedFileName + "_info", "compressed", "none");
                    chunkListImputedInfoFile.add(myChunkListImputedInfoFile);

                    GenericFile myChunkListImputedSummaryFile = new GenericFile(tmpChrDir, imputedFileName + "_summary", "decompressed",
                            "none");
                    chunkListImputedSummaryFile.add(myChunkListImputedSummaryFile);

                    GenericFile myChunkListImputedWarningsFile = new GenericFile(tmpChrDir, imputedFileName + "_warnings", "decompressed",
                            "none");
                    chunkListImputedWarningsFile.add(myChunkListImputedWarningsFile);

                    GenericFile myChunkListImputedLogFile = new GenericFile(tmpChrDir, imputedFileName + ".log", "decompressed", "none");
                    chunkListImputedLogFile.add(myChunkListImputedLogFile);

                    String filteredFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2 + "_filtered.impute";
                    GenericFile myChunkListFilteredFile = new GenericFile(tmpChrDir, filteredFileName, "compressed", "none");
                    chunkListFilteredFile.add(myChunkListFilteredFile);

                    String filteredLogFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2 + "_filtered.impute.log";
                    GenericFile myChunkListFilteredLogFile = new GenericFile(tmpChrDir, filteredLogFileName, "decompressed", "none");
                    chunkListFilteredLogFile.add(myChunkListFilteredLogFile);

                    // String filteredRsIdFileName = parsingArgs.getGenFileName(chromo) + "_" + rPanel + "_" + lim1
                    // + "_" + lim2 +"_filtered_rsid.txt";
                    String filteredRsIdFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2 + "_filtered_rsid.txt";
                    GenericFile myChunkListFilteredRsIdFile = new GenericFile(tmpChrDir, filteredRsIdFileName, "compressed", "none");
                    chunkListFilteredRsIdFile.add(myChunkListFilteredRsIdFile);

                    lim1 = lim1 + chunkSize;
                    lim2 = lim2 + chunkSize;
                }
                chromoListImputedFile.add(chunkListImputedFile);
                chromoListImputedInfoFile.add(chunkListImputedInfoFile);
                chromoListImputedSummaryFile.add(chunkListImputedSummaryFile);
                chromoListImputedWarningsFile.add(chunkListImputedWarningsFile);
                chromoListImputedLogFile.add(chunkListImputedLogFile);
                chromoListFilteredFile.add(chunkListFilteredFile);
                chromoListFilteredLogFile.add(chunkListFilteredLogFile);
                chromoListFilteredRsIdFile.add(chunkListFilteredRsIdFile);
            }
            this.imputedOutDir.add(chromoListImputedOutDir);
            this.imputedFile.add(chromoListImputedFile);

            this.imputedInfoFile.add(chromoListImputedInfoFile);
            this.imputedSummaryFile.add(chromoListImputedSummaryFile);
            this.imputedWarningsFile.add(chromoListImputedWarningsFile);
            this.imputedLogFile.add(chromoListImputedLogFile);
            this.filteredFile.add(chromoListFilteredFile);
            this.filteredFileLogFile.add(chromoListFilteredLogFile);
            this.filteredFileRsIdFile.add(chromoListFilteredRsIdFile);
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
            ArrayList<String> chromoListImputedMMOutDir = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMInfoFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMDraftFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMErateFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMRecFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMDoseFile = new ArrayList<>();
            ArrayList<ArrayList<GenericFile>> chromoListImputedMMLogFile = new ArrayList<>();
            // ArrayList<ArrayList<GenericFile>> chromoListFilteredFile = new ArrayList<>();
            // ArrayList<ArrayList<GenericFile>> chromoListFilteredLogFile = new ArrayList<>();
            // ArrayList<ArrayList<GenericFile>> chromoListFilteredRsIdFile = new ArrayList<>();

            for (int i = this.startChr; i <= this.endChr; i++) {
                int chromo = i;
                int lim1 = 1;
                int lim2 = lim1 + chunkSize - 1;
                String tmpChrDir = mixOutDir + File.separator + "Chr_" + chromo;
                chromoListImputedMMOutDir.add(tmpChrDir);

                ArrayList<GenericFile> chunkListImputedMMFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedMMInfoFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedMMDraftFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedMMErateFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedMMRecFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedMMDoseFile = new ArrayList<>();
                ArrayList<GenericFile> chunkListImputedMMLogFile = new ArrayList<>();

                int numberOfChunks = ChromoInfo.getMaxSize(chromo) / chunkSize;
                int module = ChromoInfo.getMaxSize(chromo) % chunkSize;
                if (module != 0) {
                    numberOfChunks++;
                }

                for (int k = 0; k < numberOfChunks; k++) {
                    // Now we have to create the impute files for this iteration
                    // String imputedFileName = parsingArgs.getGenFileName(chromo) + "_" + rPanel + "_" + lim1 + "_"
                    // + lim2 +".impute";
                    String imputedMMFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2 + "_minimac";
                    GenericFile myChunkListImputedMMFile = new GenericFile(tmpChrDir, imputedMMFileName, "decompressed", "none");
                    chunkListImputedMMFile.add(myChunkListImputedMMFile);

                    GenericFile myChunkListImputedMMInfoFile = new GenericFile(tmpChrDir, imputedMMFileName + ".info.gz", "compressed",
                            "none");
                    chunkListImputedMMInfoFile.add(myChunkListImputedMMInfoFile);

                    GenericFile myChunkListImputedMMDraftFile = new GenericFile(tmpChrDir, imputedMMFileName + ".info.draft.gz",
                            "compressed", "none");
                    chunkListImputedMMDraftFile.add(myChunkListImputedMMDraftFile);

                    GenericFile myChunkListImputedMMErateFile = new GenericFile(tmpChrDir, imputedMMFileName + ".erate.gz", "compressed",
                            "none");
                    chunkListImputedMMErateFile.add(myChunkListImputedMMErateFile);

                    GenericFile myChunkListImputedMMRecFile = new GenericFile(tmpChrDir, imputedMMFileName + ".rec.gz", "compressed",
                            "none");
                    chunkListImputedMMRecFile.add(myChunkListImputedMMRecFile);

                    GenericFile myChunkListImputedMMDoseFile = new GenericFile(tmpChrDir, imputedMMFileName + ".dose.gz", "compressed",
                            "none");
                    chunkListImputedMMDoseFile.add(myChunkListImputedMMDoseFile);

                    GenericFile myChunkListImputedMMLogFile = new GenericFile(tmpChrDir, imputedMMFileName + ".log", "decompressed",
                            "none");
                    chunkListImputedMMLogFile.add(myChunkListImputedMMLogFile);

                    // String filteredFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
                    // +"_filtered.impute.gz";
                    // GenericFile myChunkListFilteredFile = new GenericFile(tmpChrDir, filteredFileName,
                    // "compressed", "none");
                    // chunkListFilteredFile.add(myChunkListFilteredFile);

                    // String filteredLogFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
                    // +"_filtered.impute.log";
                    // GenericFile myChunkListFilteredLogFile = new GenericFile(tmpChrDir, filteredLogFileName,
                    // "decompressed", "none");
                    // chunkListFilteredLogFile.add(myChunkListFilteredLogFile);

                    // //String filteredRsIdFileName = parsingArgs.getGenFileName(chromo) + "_" + rPanel + "_" +
                    // lim1 + "_" + lim2 +"_filtered_rsid.txt";
                    // String filteredRsIdFileName = "chr_" + chromo + "_mixed_" + rPanel + "_" + lim1 + "_" + lim2
                    // +"_filtered_rsid.txt";
                    // GenericFile myChunkListFilteredRsIdFile = new GenericFile(tmpChrDir, filteredRsIdFileName,
                    // "compressed", "none");
                    // chunkListFilteredRsIdFile.add(myChunkListFilteredRsIdFile);

                    lim1 = lim1 + chunkSize;
                    lim2 = lim2 + chunkSize;
                }
                chromoListImputedMMFile.add(chunkListImputedMMFile);
                chromoListImputedMMInfoFile.add(chunkListImputedMMInfoFile);
                chromoListImputedMMDraftFile.add(chunkListImputedMMDraftFile);
                chromoListImputedMMErateFile.add(chunkListImputedMMErateFile);
                chromoListImputedMMRecFile.add(chunkListImputedMMRecFile);
                chromoListImputedMMDoseFile.add(chunkListImputedMMDoseFile);
                chromoListImputedMMLogFile.add(chunkListImputedMMLogFile);
                // chromoListFilteredFile.add(chunkListFilteredFile);
                // chromoListFilteredLogFile.add(chunkListFilteredLogFile);
                // chromoListFilteredRsIdFile.add(chunkListFilteredRsIdFile);
            }
            this.imputedMMOutDir.add(chromoListImputedMMOutDir);
            this.imputedMMFile.add(chromoListImputedMMFile);

            this.imputedMMInfoFile.add(chromoListImputedMMInfoFile);
            this.imputedMMDraftFile.add(chromoListImputedMMDraftFile);
            this.imputedMMErateFile.add(chromoListImputedMMErateFile);
            this.imputedMMRecFile.add(chromoListImputedMMRecFile);
            this.imputedMMDoseFile.add(chromoListImputedMMDoseFile);
            this.imputedMMLogFile.add(chromoListImputedMMLogFile);
            // filteredFile.add(chromoListFilteredFile);
            // filteredFileLogFile.add(chromoListFilteredLogFile);
            // filteredFileRsIdFile.add(chromoListFilteredRsIdFile);
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
     * Method to access the imputedFileName
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getImputedFileName(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        // The offset is because the array start in position 0.
        return this.imputedFile.get(rPanelIndex).get(indexChr).get(indexChunk).getName();
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
     * Method to set the finalStatus of the imputedFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setImputedFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
    public void setImputedInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
     * Method to access the imputedLogFileName
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getImputedLogFileName(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        return this.imputedLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).getName();
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
     * Method to access the imputedMMFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getImputedMMFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        return this.imputedMMFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
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
    public void setImputedMMFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedMMFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
     * Method to set the finalStatus of the imputedFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setImputedMMInfoFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - startChr;
        int indexChunk = lim1 / chunkSize;
        imputedMMInfoFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the imputedMMDraftFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getImputedMMDraftFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        return this.imputedMMDraftFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
    }

    /**
     * Method to set the finalStatus of the imputeDraftFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setImputedMMDraftFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedMMDraftFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
     * Method to set the finalStatus of the imputedMMErateFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setImputedMMErateFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedMMErateFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
     * Method to set the finalStatus of the imputeMMRecFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setImputedMMRecFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedMMRecFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the imputedMMDoseFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getImputedMMDoseFile(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        return this.imputedMMDoseFile.get(rPanelIndex).get(indexChr).get(indexChunk).getFullName();
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
    public void setImputedMMDoseFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedMMDoseFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
     * Method to set the finalStatus of the imputedMMLogFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setImputedMMLogFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.imputedMMLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the filteredFileName
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getFilteredFileName(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        // The offset is because the array start in position 0.
        return this.filteredFile.get(rPanelIndex).get(indexChr).get(indexChunk).getName();
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
     * Method to set the finalStatus of filteredFile
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @param finalStatus
     */
    public void setFilteredFileFinalStatus(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        this.filteredFile.get(rPanelIndex).get(indexChr).get(indexChunk).setFinalStatus(finalStatus);
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
     * Method to access the filteredLogFileName
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getFilteredLogFileName(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        // The offset is because the array start in position 0.
        return this.filteredFileLogFile.get(rPanelIndex).get(indexChr).get(indexChunk).getName();
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
     * Method to access the filteredRsIdFileName
     * 
     * @param rPanelIndex
     * @param chromo
     * @param lim1
     * @param lim2
     * @param chunkSize
     * @return
     */
    public String getFilteredRsIdFileName(int rPanelIndex, int chromo, int lim1, int lim2, int chunkSize) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // Check limits are within the bounds
        checkLimits(chromo, lim1, lim2);

        int indexChr = chromo - this.startChr;
        int indexChunk = lim1 / chunkSize;
        // The offset is because the array start in position 0.
        return this.filteredFileRsIdFile.get(rPanelIndex).get(indexChr).get(indexChunk).getName();
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
        LOGGER.info("OutDir[" + rPanelIndex + "][Chr_" + chromo + "]=" + this.imputedOutDir.get(rPanelIndex).get(indexChr));

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
            LOGGER.fatal("[ImputationFiles] Error, Chunk " + lim1 + "_" + lim2 + "does not exist for chromosome " + chromo);
            System.exit(1);
        }
    }

}
