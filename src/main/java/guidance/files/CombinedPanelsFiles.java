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
import java.util.ArrayList;

import java.io.File;

import guidance.utils.ChromoInfo;
import guidance.utils.ParseCmdLine;


public class CombinedPanelsFiles {

    private ArrayList<String> testTypeListOutDir = new ArrayList<>();

    // Test-combined files per chromosome
    private ArrayList<ArrayList<GenericFile>> testTypeCombinedFilteredByAllChromoFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeCombinedFilteredByAllXChromoFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeCombinedCondensedChromoFile = new ArrayList<>();

    // Test-combined global files
    private ArrayList<GenericFile> testTypeCombinedFilteredByAllFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeCombinedFilteredByAllXFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeCombinedCondensedFile = new ArrayList<>();

    private ArrayList<GenericFile> testTypeTopHitsFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeCorrectedPvaluesFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeQqPlotPdfFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeManhattanPdfFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeQqPlotTiffFile = new ArrayList<>();
    private ArrayList<GenericFile> testTypeManhattanTiffFile = new ArrayList<>();


    /**
     * CombinedPanelsFiles new Instance
     * 
     * @param parsingArgs
     * @param baseOutDir
     * @param refPanels
     */
    public CombinedPanelsFiles(ParseCmdLine parsingArgs, String baseOutDir, List<String> refPanels) {
        int startChr = parsingArgs.getStart();
        int endChr = parsingArgs.getEnd();
        int endChrNormal = endChr;
        if (endChr == 23) {
            endChrNormal = 22;
        }

        String startChrS = Integer.toString(startChr);
        String endChrS = Integer.toString(endChr);
        String endChrNormalS = Integer.toString(endChrNormal);

        int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();

        // We create the first directory name: the cohort directory
        String mixedCohort = parsingArgs.getCohort();

        for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
            String testTypeName = parsingArgs.getTestTypeName(tt);
            String testTypeOutDir = baseOutDir + File.separator + "associations" + File.separator + testTypeName + File.separator
                    + mixedCohort + "_combined_panels";

            ArrayList<GenericFile> combinedFilteredByAllFile = new ArrayList<>();
            ArrayList<GenericFile> combinedFilteredByAllXFile = new ArrayList<>();
            ArrayList<GenericFile> combinedCondensedAllFile = new ArrayList<>();

            String rPanel = null;
            for (int j = 0; j < refPanels.size(); j++) {
                rPanel = refPanels.get(j);
                testTypeOutDir = testTypeOutDir + "_" + rPanel;
            }

            rPanel = refPanels.get(0);
            String prefixFilteredName = "filteredByAll_results_" + testTypeName + "_" + mixedCohort;
            String prefixCondensedName = "condensed_results_" + testTypeName + "_" + mixedCohort;
            String prefixTopHitsName = "tophits_" + testTypeName + "_" + mixedCohort;
            String prefixQqPlotName = "QQplot_" + testTypeName + "_" + mixedCohort;
            String prefixManhattanName = "manhattan_" + testTypeName + "_" + mixedCohort;
            String prefixCorrectedPvaluesName = "corrected_pvalues_" + testTypeName + "_" + mixedCohort;

            for (int chr = startChr; chr <= endChr; ++chr) {
                int minSize = ChromoInfo.getMinSize(chr);
                int maxSize = ChromoInfo.getMaxSize(chr);
                int chunkSize = parsingArgs.getChunkSize();

                for (int j = minSize; j < maxSize; j = j + chunkSize) {

                    String tmpCombinedFilteredByAllFileName = null;
                    if (startChr == endChrNormal) {
                        tmpCombinedFilteredByAllFileName = prefixFilteredName + "_chr_" + startChrS + "_chunk" + j + ".txt.gz";
                    } else {
                        tmpCombinedFilteredByAllFileName = prefixFilteredName + "_chr_" + startChrS + "_to_" + endChrNormalS + "_chunk" + j
                                + ".txt.gz";
                    }
                    GenericFile myCombinedFilteredByAllFile = new GenericFile(testTypeOutDir, tmpCombinedFilteredByAllFileName,
                            "compressed", "none");
                    combinedFilteredByAllFile.add(myCombinedFilteredByAllFile);

                    // If we are going to process chr 23, then prepare file names for it.
                    if (chr == 23) {
                        String tmpCombinedFilteredByAllXFileName = prefixFilteredName + "_chr_" + chr + "_chunk" + j + ".txt.gz";
                        GenericFile myCombinedFilteredByAllXFile = new GenericFile(testTypeOutDir, tmpCombinedFilteredByAllXFileName,
                                "compressed", "none");
                        combinedFilteredByAllXFile.add(myCombinedFilteredByAllXFile);
                    }

                    String tmpCombinedCondensedFileName = null;
                    if (startChr == endChrNormal) {
                        tmpCombinedCondensedFileName = prefixCondensedName + "_chr_" + startChrS + "_chunk" + j + ".txt.gz";
                    } else {
                        tmpCombinedCondensedFileName = prefixCondensedName + "_chr_" + startChrS + "_to_" + endChrNormalS + "_chunk" + j
                                + ".txt.gz";
                    }
                    GenericFile myCombinedCondensedFile = new GenericFile(testTypeOutDir, tmpCombinedCondensedFileName, "compressed",
                            "none");
                    combinedCondensedAllFile.add(myCombinedCondensedFile);
                } // End of for chunk
            } // End of for chromo

            // Now we have to build the list of reduced files for the type of Test. We store this list
            this.testTypeListOutDir.add(testTypeOutDir);
            this.testTypeCombinedFilteredByAllChromoFile.add(combinedFilteredByAllFile);
            if (endChr == 23) {
                this.testTypeCombinedFilteredByAllXChromoFile.add(combinedFilteredByAllXFile);
            }
            this.testTypeCombinedCondensedChromoFile.add(combinedCondensedAllFile);

            String tmpCombinedFilteredByAllFileName = null;
            if (startChr == endChr) {
                tmpCombinedFilteredByAllFileName = prefixFilteredName + "_chr_" + startChrS + ".txt.gz";
            } else {
                tmpCombinedFilteredByAllFileName = prefixFilteredName + "_chr_" + startChrS + "_to_" + endChrS + ".txt.gz";
            }
            GenericFile myFilteredByAllFile = new GenericFile(testTypeOutDir, tmpCombinedFilteredByAllFileName, "compressed", "none");
            this.testTypeCombinedFilteredByAllFile.add(myFilteredByAllFile);

            if (endChr == 23) {
                String tmpCombinedFilteredByAllXFileName = null;
                if (startChr == endChr) {
                    tmpCombinedFilteredByAllXFileName = prefixFilteredName + "_chr_" + startChrS + ".txt.gz";
                } else {
                    tmpCombinedFilteredByAllXFileName = prefixFilteredName + "_chr_" + startChrS + "_to_" + endChrS + ".txt.gz";
                }
                GenericFile myFilteredByAllXFile = new GenericFile(testTypeOutDir, tmpCombinedFilteredByAllXFileName, "compressed", "none");
                this.testTypeCombinedFilteredByAllXFile.add(myFilteredByAllXFile);
            }

            String tmpCombinedCondensedFileName = null;
            if (startChr == endChr) {
                tmpCombinedCondensedFileName = prefixCondensedName + "_chr_" + startChrS + ".txt.gz";
            } else {
                tmpCombinedCondensedFileName = prefixCondensedName + "_chr_" + startChrS + "_to_" + endChrS + ".txt.gz";
            }
            GenericFile myCombinedCondensedFile = new GenericFile(testTypeOutDir, tmpCombinedCondensedFileName, "compressed", "none");
            this.testTypeCombinedCondensedFile.add(myCombinedCondensedFile);

            String tmpTopHitsFileName = null;
            if (startChr == endChr) {
                tmpTopHitsFileName = prefixTopHitsName + "_chr_" + startChrS + ".txt.gz";
            } else {
                tmpTopHitsFileName = prefixTopHitsName + "_chr_" + startChrS + "_to_" + endChrS + ".txt.gz";
            }
            GenericFile myTopHitsFile = new GenericFile(testTypeOutDir, tmpTopHitsFileName, "compressed", "none");
            this.testTypeTopHitsFile.add(myTopHitsFile);

            String tmpQqPlotPdfFileName = null;
            if (startChr == endChr) {
                tmpQqPlotPdfFileName = prefixQqPlotName + "_chr_" + startChrS + ".pdf";
            } else {
                tmpQqPlotPdfFileName = prefixQqPlotName + "_chr_" + startChrS + "_to_" + endChrS + ".pdf";
            }
            GenericFile myQqPlotPdfFile = new GenericFile(testTypeOutDir, tmpQqPlotPdfFileName, "compressed", "none");
            this.testTypeQqPlotPdfFile.add(myQqPlotPdfFile);

            String tmpQqPlotTiffFileName = null;
            if (startChr == endChr) {
                tmpQqPlotTiffFileName = prefixQqPlotName + "_chr_" + startChrS + ".tiff";
            } else {
                tmpQqPlotTiffFileName = prefixQqPlotName + "_chr_" + startChrS + "_to_" + endChrS + ".tiff";
            }
            GenericFile myQqPlotTiffFile = new GenericFile(testTypeOutDir, tmpQqPlotTiffFileName, "compressed", "none");
            this.testTypeQqPlotTiffFile.add(myQqPlotTiffFile);

            String tmpManhattanPdfFileName = null;
            if (startChr == endChr) {
                tmpManhattanPdfFileName = prefixManhattanName + "_chr_" + startChrS + ".pdf";
            } else {
                tmpManhattanPdfFileName = prefixManhattanName + "_chr_" + startChrS + "_to_" + endChrS + ".pdf";
            }
            GenericFile myManhattanPdfFile = new GenericFile(testTypeOutDir, tmpManhattanPdfFileName, "compressed", "none");
            this.testTypeManhattanPdfFile.add(myManhattanPdfFile);

            String tmpManhattanTiffFileName = null;
            if (startChr == endChr) {
                tmpManhattanTiffFileName = prefixManhattanName + "_chr_" + startChrS + ".tiff";
            } else {
                tmpManhattanTiffFileName = prefixManhattanName + "_chr_" + startChrS + "_to_" + endChrS + ".tiff";
            }
            GenericFile myManhattanTiffFile = new GenericFile(testTypeOutDir, tmpManhattanTiffFileName, "compressed", "none");
            this.testTypeManhattanTiffFile.add(myManhattanTiffFile);

            String tmpCorrectedPvaluesFileName = null;
            if (startChr == endChr) {
                tmpCorrectedPvaluesFileName = prefixCorrectedPvaluesName + "_chr_" + startChrS + ".txt";
            } else {
                tmpCorrectedPvaluesFileName = prefixCorrectedPvaluesName + "_chr_" + startChrS + "_to_" + endChrS + ".txt";
            }

            GenericFile myCorrectedPvaluesFile = new GenericFile(testTypeOutDir, tmpCorrectedPvaluesFileName, "compressed", "none");
            this.testTypeCorrectedPvaluesFile.add(myCorrectedPvaluesFile);

        } // End of for test types
    }

    /**
     * Method to access correctedPvaluesFile information
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCombinedOutDir(int testTypeIndex) {
        return this.testTypeListOutDir.get(testTypeIndex);
    }

    /**
     * Method to access testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllFile(int testTypeIndex) {
        return this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeCorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @param index
     * @param finalStatus
     */
    public void setCombinedFilteredByAllFileFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllFileFinalStatus(int testTypeIndex) {
        return this.testTypeCombinedFilteredByAllFile.get(testTypeIndex).getFinalStatus();
    }

    /**
     * Method to access testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllXFile(int testTypeIndex) {
        return this.testTypeCombinedFilteredByAllXFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeCorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @param index
     * @param finalStatus
     */
    public void setCombinedFilteredByAllXFileFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeCombinedFilteredByAllXFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllXFileFinalStatus(int testTypeIndex) {
        return this.testTypeCombinedFilteredByAllXFile.get(testTypeIndex).getFinalStatus();
    }

    /**
     * Method to access testTypeCombinedCondensedFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCombinedCondensedFile(int testTypeIndex) {
        return this.testTypeCombinedCondensedFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeCombinedCondensedFile
     * 
     * @param testTypeIndex
     * @param finalStatus
     */
    public void setCombinedCondensedFileFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeCombinedCondensedFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeCombinedCondensedFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCombinedCondensedFileFinalStatus(int testTypeIndex) {
        return this.testTypeCombinedCondensedFile.get(testTypeIndex).getFinalStatus();
    }

    /**
     * Method to access testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllChromoFile(int testTypeIndex, int index) {
        return this.testTypeCombinedFilteredByAllChromoFile.get(testTypeIndex).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeCorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @param index
     * @param finalStatus
     */
    public void setCombinedFilteredByAllFileFinalStatus(int testTypeIndex, int index, String finalStatus) {
        this.testTypeCombinedFilteredByAllChromoFile.get(testTypeIndex).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllFileFinalStatus(int testTypeIndex, int index) {
        return this.testTypeCombinedFilteredByAllChromoFile.get(testTypeIndex).get(index).getFinalStatus();
    }

    /**
     * Method to access testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllXChromoFile(int testTypeIndex, int index) {
        return this.testTypeCombinedFilteredByAllXChromoFile.get(testTypeIndex).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeCorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @param index
     * @param finalStatus
     */
    public void setCombinedFilteredByAllXChromoFileFinalStatus(int testTypeIndex, int index, String finalStatus) {
        this.testTypeCombinedFilteredByAllXChromoFile.get(testTypeIndex).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeCombinedFilteredByAllFile
     * 
     * @param testTypeIndex
     * @param index
     * @return
     */
    public String getCombinedFilteredByAllXChromoFileFinalStatus(int testTypeIndex, int index) {
        return this.testTypeCombinedFilteredByAllXChromoFile.get(testTypeIndex).get(index).getFinalStatus();
    }

    /**
     * Method to access testTypeCombinedCondensedFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCombinedCondensedChromoFile(int testTypeIndex, int index) {
        return this.testTypeCombinedCondensedChromoFile.get(testTypeIndex).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeCombinedCondensedFile
     * 
     * @param testTypeIndex
     * @param finalStatus
     */
    public void setCombinedCondensedChromoFileFinalStatus(int testTypeIndex, String finalStatus, int index) {
        this.testTypeCombinedCondensedChromoFile.get(testTypeIndex).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeCombinedCondensedFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCombinedCondensedChromoFileFinalStatus(int testTypeIndex, int index) {
        return this.testTypeCombinedCondensedChromoFile.get(testTypeIndex).get(index).getFinalStatus();
    }

    /**
     * Method to access testTypeTopHitsFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getTopHitsFile(int testTypeIndex) {
        return this.testTypeTopHitsFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the testTypeTopHitsFile
     * 
     * @param testTypeIndex
     * @param finalStatus
     */
    public void setTopHitsFileFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeTopHitsFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the testTypeTopHitsFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getTopHitsFileFinalStatus(int testTypeIndex) {
        return this.testTypeTopHitsFile.get(testTypeIndex).getFinalStatus();
    }

    /**
     * Method to access qqPlotPdfFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getQqPlotPdfFile(int testTypeIndex) {
        return this.testTypeQqPlotPdfFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the qqPlotPdfFile
     * 
     * @param testTypeIndex
     * @param finalStatus
     */
    public void setQqPlotPdfFileFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeQqPlotPdfFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the qqPlotPdfFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getQqPlotPdfFileFinalStatus(int testTypeIndex) {
        return this.testTypeQqPlotPdfFile.get(testTypeIndex).getFinalStatus();
    }

    /**
     * Method to access qqPlotTiffFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getQqPlotTiffFile(int testTypeIndex) {
        return this.testTypeQqPlotTiffFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the qqPlotPdfFile
     * 
     * @param testTypeIndex
     * @param finalStatus
     */
    public void setQqPlotTiffFileFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeQqPlotTiffFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the qqPlotTiffFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getQqPlotTiffFileFinalStatus(int testTypeIndex) {
        return this.testTypeQqPlotTiffFile.get(testTypeIndex).getFinalStatus();
    }

    /**
     * Method to access manhattanPdfFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getManhattanPdfFile(int testTypeIndex) {
        return this.testTypeManhattanPdfFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to access manhattanTiffFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getManhattanTiffFile(int testTypeIndex) {
        return this.testTypeManhattanTiffFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to access CorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCorrectedPvaluesFile(int testTypeIndex) {
        return this.testTypeCorrectedPvaluesFile.get(testTypeIndex).getFullName();
    }

    /**
     * Method to set the finalStatus of the CorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @param finalStatus
     */
    public void setCorrectedPvaluesFinalStatus(int testTypeIndex, String finalStatus) {
        this.testTypeCorrectedPvaluesFile.get(testTypeIndex).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the CorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @return
     */
    public String getCorrectedPvaluesFileFinalStatus(int testTypeIndex) {
        return this.testTypeCorrectedPvaluesFile.get(testTypeIndex).getFinalStatus();
    }

}
