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


public class ResultsFiles {

    private ArrayList<ArrayList<String>> testTypeListOutDir = new ArrayList<>();

    private ArrayList<ArrayList<GenericFile>> testTypeTopHitsFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeCrossRangesFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeCorrectedPvaluesFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeQqPlotPdfFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeManhattanPdfFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeQqPlotTiffFile = new ArrayList<>();
    private ArrayList<ArrayList<GenericFile>> testTypeManhattanTiffFile = new ArrayList<>();


    /**
     * New ResultFiles instance
     * 
     * @param parsingArgs
     * @param baseOutDir
     * @param refPanels
     */
    public ResultsFiles(ParseCmdLine parsingArgs, String baseOutDir, List<String> refPanels) {
        int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();

        // We create the first directory name: the cohort directory
        String mixedCohort = parsingArgs.getCohort();

        for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
            String testTypeName = parsingArgs.getTestTypeName(tt);
            String testTypeOutDir = baseOutDir + File.separator + "associations" + File.separator + testTypeName;

            ArrayList<String> rpanelListOutDir = new ArrayList<>();

            ArrayList<GenericFile> rpanelTopHitsFile = new ArrayList<>();
            ArrayList<GenericFile> rpanelCrossRangesFile = new ArrayList<>();
            ArrayList<GenericFile> rpanelCorrectedPvaluesFile = new ArrayList<>();

            ArrayList<GenericFile> rpanelQqPlotPdfFile = new ArrayList<>();
            ArrayList<GenericFile> rpanelManhattanPdfFile = new ArrayList<>();

            ArrayList<GenericFile> rpanelQqPlotTiffFile = new ArrayList<>();
            ArrayList<GenericFile> rpanelManhattanTiffFile = new ArrayList<>();

            for (int j = 0; j < refPanels.size(); j++) {
                String rPanel = refPanels.get(j);
                String rpanelOutDirSummary = testTypeOutDir + File.separator + mixedCohort + "_for_" + rPanel + File.separator + "summary";

                String tmpTopHitsFileName = "tophits_" + testTypeName + "_" + mixedCohort + "_" + rPanel + ".txt.gz";
                GenericFile myTopHitsFile = new GenericFile(rpanelOutDirSummary, tmpTopHitsFileName, "compressed", "none");
                rpanelTopHitsFile.add(myTopHitsFile);
                
                String tmpCrossRangesFileName = "tophits_" + testTypeName + "_" + mixedCohort + "_" + rPanel + "_crossmodel_ranges.txt.gz";
                GenericFile myCrossRangesFile = new GenericFile(rpanelOutDirSummary, tmpCrossRangesFileName, "compressed", "none");
                rpanelCrossRangesFile.add(myCrossRangesFile);

                String tmpCorrectedPvaluesFileName = "corrected_pvalues_" + testTypeName + "_" + mixedCohort + "_" + rPanel + ".txt";
                GenericFile myCorrectedPvaluesFile = new GenericFile(rpanelOutDirSummary, tmpCorrectedPvaluesFileName, "uncompressed",
                        "none");
                rpanelCorrectedPvaluesFile.add(myCorrectedPvaluesFile);

                String tmpQqPlotPdfFileName = "QQplot_" + testTypeName + "_" + mixedCohort + "_" + rPanel + ".pdf";
                GenericFile myQqPlotPdfFile = new GenericFile(rpanelOutDirSummary, tmpQqPlotPdfFileName, "uncompressed", "none");
                rpanelQqPlotPdfFile.add(myQqPlotPdfFile);

                String tmpManhattanPdfFileName = "manhattan_" + testTypeName + "_" + mixedCohort + "_" + rPanel + ".pdf";
                GenericFile myManhattanPdfFile = new GenericFile(rpanelOutDirSummary, tmpManhattanPdfFileName, "uncompressed", "none");
                rpanelManhattanPdfFile.add(myManhattanPdfFile);

                String tmpQqPlotTiffFileName = "QQplot_" + testTypeName + "_" + mixedCohort + "_" + rPanel + ".tiff";
                GenericFile myQqPlotTiffFile = new GenericFile(rpanelOutDirSummary, tmpQqPlotTiffFileName, "uncompressed", "none");
                rpanelQqPlotTiffFile.add(myQqPlotTiffFile);

                String tmpManhattanTiffFileName = "manhattan_" + testTypeName + "_" + mixedCohort + "_" + rPanel + ".tiff";
                GenericFile myManhattanTiffFile = new GenericFile(rpanelOutDirSummary, tmpManhattanTiffFileName, "uncompressed", "none");
                rpanelManhattanTiffFile.add(myManhattanTiffFile);
            } // End of for panels

            // Now we have to build the list of reduced files for the type of Test. We store this list
            testTypeListOutDir.add(rpanelListOutDir);
            testTypeTopHitsFile.add(rpanelTopHitsFile);
            testTypeCrossRangesFile.add(rpanelCrossRangesFile);
            testTypeCorrectedPvaluesFile.add(rpanelCorrectedPvaluesFile);

            testTypeQqPlotPdfFile.add(rpanelQqPlotPdfFile);
            testTypeManhattanPdfFile.add(rpanelManhattanPdfFile);

            testTypeQqPlotTiffFile.add(rpanelQqPlotTiffFile);
            testTypeManhattanTiffFile.add(rpanelManhattanTiffFile);

        } // End of for test types
    }

    /**
     * Method to access correctedPvaluesFile information
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getListOutDir(int testTypeIndex, int rPanelIndex) {
        return this.testTypeListOutDir.get(testTypeIndex).get(rPanelIndex);
    }

    /**
     * Method to access topHitsFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getTopHitsFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeTopHitsFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }
    
    /**
     * Method to access crossRangesFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getCrossRangesFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeCrossRangesFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }

    /**
     * Method to access CorrectedPvaluesFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getCorrectedPvaluesFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeCorrectedPvaluesFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }

    /**
     * Method to access qqPlotPdfFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getQqPlotPdfFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeQqPlotPdfFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }

    /**
     * Method to access qqPlotTiffFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getQqPlotTiffFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeQqPlotTiffFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }

    /**
     * Method to access manhattanPdfFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getManhattanPdfFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeManhattanPdfFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }

    /**
     * Method to access manhattanTiffFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getManhattanTiffFile(int testTypeIndex, int rPanelIndex) {
        return this.testTypeManhattanTiffFile.get(testTypeIndex).get(rPanelIndex).getFullName();
    }

}
