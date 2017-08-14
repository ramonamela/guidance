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

import java.util.ArrayList;


public class MergeFiles {

    private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;

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

    private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedFilteredXFile = new ArrayList<>();
    // private ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>> testTypeCombinedReducedCondensedXFile = new
    // ArrayList<ArrayList<ArrayList<ArrayList<GenericFile>>>>();

    private int startChr = 0;
    private int endChr = 0;
    // private int endChrNormal=endChr;


    /**
     * New MergeFiles instance
     * 
     * @param parsingArgs
     * @param generalChromoInfo
     * @param baseOutDir
     * @param refPanels
     */
    public MergeFiles(ParseCmdLine parsingArgs, ChromoInfo generalChromoInfo, String baseOutDir, List<String> refPanels) {
        int i;
        int j;

        startChr = parsingArgs.getStart();
        endChr = parsingArgs.getEnd();

        // int endChrNormal =0;
        // for(i=startChr; i<=endChr;i++) {
        // if(i < 23 ) {
        // endChrNormal = i;
        // }
        // }

        int chunkSize = parsingArgs.getChunkSize();
        int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();

        /** We create the first directory name: the cohort directory */
        String mixedCohort = parsingArgs.getCohort();

        for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
            String testTypeName = parsingArgs.getTestTypeName(tt);
            String testTypeOutDir = baseOutDir + "/associations/" + testTypeName;
            String testTypeOutDir2 = baseOutDir + "/associations/" + testTypeName + "/" + mixedCohort + "_combined_panels";

            String rPanel = null;
            for (j = 0; j < refPanels.size(); j++) {
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

            ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedFilteredXFile = new ArrayList<>();
            // ArrayList<ArrayList<ArrayList<GenericFile>>> rpanelCombinedReducedCondensedXFile = new
            // ArrayList<ArrayList<ArrayList<GenericFile>>>();

            for (j = 0; j < refPanels.size(); j++) {
                rPanel = refPanels.get(j);
                String rpanelOutDir = testTypeOutDir + "/" + mixedCohort + "_for_" + rPanel;
                String rpanelOutDirSummary = testTypeOutDir + "/" + mixedCohort + "_for_" + rPanel + "/summary";

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

                ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedFilteredXFile = new ArrayList<>();
                // ArrayList<ArrayList<GenericFile>> chromoListCombinedReducedCondensedXFile = new
                // ArrayList<ArrayList<GenericFile>>();

                boolean addExtraCount = false;
                for (i = startChr; i <= endChr; i++) {
                    int chromo = i;

                    if (i != startChr) {

                        System.out.println("mon :: TRUE :: i = " + i);
                        addExtraCount = true;
                    }

                    String tmpChrDir = rpanelOutDir + "/" + "Chr_" + chromo;
                    chromoListOutDir.add(tmpChrDir);

                    int maxSize = generalChromoInfo.getMaxSize(chromo);
                    int total_chunks = maxSize / chunkSize;
                    int module = maxSize % chunkSize;
                    if (module != 0)
                        total_chunks++;
                    int lim1 = 1;
                    int lim2 = lim1 + chunkSize - 1;

                    ArrayList<GenericFile> listReducedFile = new ArrayList<>();
                    ArrayList<GenericFile> listReducedFilteredFile = new ArrayList<>();
                    ArrayList<GenericFile> listReducedCondensedFile = new ArrayList<>();

                    ArrayList<GenericFile> listCombinedReducedFilteredFile = new ArrayList<>();
                    ArrayList<GenericFile> listCombinedReducedCondensedFile = new ArrayList<>();

                    ArrayList<GenericFile> listCombinedReducedFilteredXFile = new ArrayList<>();
                    // ArrayList<GenericFile> listCombinedReducedCondensedXFile = new ArrayList<>();

                    // Now we have to create the rest of file names that will be used to reduce files
                    int counter = 0;

                    for (int deep = 0; deep < total_chunks - 1; deep++) {
                        String tmpReducedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_reduce_file_" + counter
                                + ".txt.gz";
                        String tmpReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_reduce_filtered_file_"
                                + counter + ".txt.gz";
                        String tmpReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_reduce_condensed_file_"
                                + counter + ".txt.gz";

                        String tmpReducedFile = tmpChrDir + "/" + tmpReducedFileName;
                        String tmpReducedFilteredFile = tmpChrDir + "/" + tmpReducedFilteredFileName;
                        String tmpReducedCondensedFile = tmpChrDir + "/" + tmpReducedCondensedFileName;

                        GenericFile myReducedFile = new GenericFile(tmpChrDir, tmpReducedFileName, "uncompressed", "none");
                        GenericFile myReducedFilteredFile = new GenericFile(tmpChrDir, tmpReducedFilteredFileName, "uncompressed", "none");
                        GenericFile myReducedCondensedFile = new GenericFile(tmpChrDir, tmpReducedCondensedFileName, "uncompressed",
                                "none");

                        listReducedFile.add(myReducedFile);
                        listReducedFilteredFile.add(myReducedFilteredFile);
                        listReducedCondensedFile.add(myReducedCondensedFile);

                        if (j == 0) {
                            if (i != 23) {
                                String tmpCombinedReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_reduce_filtered_file_"
                                        + counter + ".txt.gz";
                                // String tmpCombinedReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName +
                                // "_reduce_condensed_file_" + counter + ".txt.gz";
                                String tmpCombinedReducedFilteredFile = testTypeOutDir2 + "/" + tmpCombinedReducedFilteredFileName;
                                // String tmpCombinedReducedCondensedFile = testTypeOutDir2 + "/" +
                                // tmpCombinedReducedCondensedFileName;
                                GenericFile myCombinedReducedFilteredFile = new GenericFile(testTypeOutDir2,
                                        tmpCombinedReducedFilteredFileName, "uncompressed", "none");
                                // GenericFile myCombinedReducedCondensedFile = new GenericFile(testTypeOutDir2,
                                // tmpCombinedReducedCondensedFileName, "uncompressed", "none");

                                listCombinedReducedFilteredFile.add(myCombinedReducedFilteredFile);
                                // listCombinedReducedCondensedFile.add(myCombinedReducedCondensedFile);

                            } else if (i == 23) {
                                String tmpCombinedReducedFilteredXFileName = "chr_" + chromo + "_" + testTypeName + "_reduce_filtered_file_"
                                        + counter + ".txt.gz";
                                // String tmpCombinedReducedCondensedXFileName = "chr_" + chromo + "_" + testTypeName +
                                // "_reduce_condensed_file_" + counter + ".txt.gz";
                                String tmpCombinedReducedFilteredXFile = testTypeOutDir2 + "/" + tmpCombinedReducedFilteredXFileName;
                                // String tmpCombinedReducedCondensedXFile = testTypeOutDir2 + "/" +
                                // tmpCombinedReducedCondensedXFileName;
                                GenericFile myCombinedReducedFilteredXFile = new GenericFile(testTypeOutDir2,
                                        tmpCombinedReducedFilteredXFileName, "uncompressed", "none");
                                // GenericFile myCombinedReducedCondensedXFile = new GenericFile(testTypeOutDir2,
                                // tmpCombinedReducedCondensedXFileName, "uncompressed", "none");

                                listCombinedReducedFilteredXFile.add(myCombinedReducedFilteredXFile);
                                // listCombinedReducedCondensedXFile.add(myCombinedReducedCondensedXFile);
                            }

                            String tmpCombinedReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_reduce_condensed_file_"
                                    + counter + ".txt.gz";
                            String tmpCombinedReducedCondensedFile = testTypeOutDir2 + "/" + tmpCombinedReducedCondensedFileName;
                            GenericFile myCombinedReducedCondensedFile = new GenericFile(testTypeOutDir2,
                                    tmpCombinedReducedCondensedFileName, "uncompressed", "none");
                            listCombinedReducedCondensedFile.add(myCombinedReducedCondensedFile);

                        }

                        counter++;
                    }

                    if (i != 23)
                        System.out.println("mon :: ABANS F :: chr " + chromo + " ::  size = " + listCombinedReducedFilteredFile.size());
                    else
                        System.out.println("mon :: ABANS F :: chr " + chromo + " ::  size = " + listCombinedReducedFilteredXFile.size());
                    System.out.println("mon :: ABANS C :: chr " + chromo + " ::  size = " + listCombinedReducedCondensedFile.size());

                    if (addExtraCount) {

                        System.out.println("mon :: DINS EXTRA :: chr " + chromo + " ::  counter " + counter);

                        String tmpCombinedReducedFilteredFileName = "chr_" + chromo + "_" + testTypeName + "_reduce_filtered_file_"
                                + counter + ".txt.gz";
                        String tmpCombinedReducedCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_reduce_condensed_file_"
                                + counter + ".txt.gz";
                        String tmpCombinedReducedFilteredFile = testTypeOutDir2 + "/" + tmpCombinedReducedFilteredFileName;
                        String tmpCombinedReducedCondensedFile = testTypeOutDir2 + "/" + tmpCombinedReducedCondensedFileName;
                        GenericFile myCombinedReducedFilteredFile = new GenericFile(testTypeOutDir2, tmpCombinedReducedFilteredFileName,
                                "uncompressed", "none");
                        GenericFile myCombinedReducedCondensedFile = new GenericFile(testTypeOutDir2, tmpCombinedReducedCondensedFileName,
                                "uncompressed", "none");

                        listCombinedReducedFilteredFile.add(myCombinedReducedFilteredFile);
                        listCombinedReducedCondensedFile.add(myCombinedReducedCondensedFile);

                        addExtraCount = false;
                    }

                    if (i != 23)
                        System.out.println("mon :: DESPRES F :: chr " + chromo + " ::  size = " + listCombinedReducedFilteredFile.size());
                    else
                        System.out.println("mon :: DESPRES F :: chr " + chromo + " ::  size = " + listCombinedReducedFilteredXFile.size());
                    System.out.println("mon :: DESPRES C :: chr " + chromo + " ::  size = " + listCombinedReducedCondensedFile.size());

                    System.out.println("\n\n");

                    // Now we have built the list of reduce file for this chromosome. We store this list
                    chromoListReducedFile.add(listReducedFile);
                    chromoListReducedFilteredFile.add(listReducedFilteredFile);
                    chromoListReducedCondensedFile.add(listReducedCondensedFile);

                    chromoListCombinedReducedFilteredFile.add(listCombinedReducedFilteredFile);
                    chromoListCombinedReducedCondensedFile.add(listCombinedReducedCondensedFile);
                    chromoListCombinedReducedFilteredXFile.add(listCombinedReducedFilteredXFile);
                    // chromoListCombinedReducedCondensedXFile.add(listCombinedReducedCondensedXFile);

                    String tmpFilteredByAllFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel
                            + "_filtered_by_maf_info_hwe.txt.gz";
                    String tmpFilteredByAllFile = tmpChrDir + "/" + tmpFilteredByAllFileName;
                    GenericFile myFilteredByAllFile = new GenericFile(tmpChrDir, tmpFilteredByAllFileName, "uncompressed", "none");
                    chromoFilteredByAllFile.add(myFilteredByAllFile);

                    String tmpCondensedFileName = "chr_" + chromo + "_" + testTypeName + "_" + rPanel + "_condensed.txt.gz";
                    String tmpCondensedFile = tmpChrDir + "/" + tmpCondensedFileName;
                    GenericFile myCondensedFile = new GenericFile(tmpChrDir, tmpCondensedFileName, "uncompressed", "none");
                    chromoCondensedFile.add(myCondensedFile);

                } // end for(i=startChr; i<=endChr;i++)
                rpanelListOutDir.add(chromoListOutDir);
                rpanelReducedFile.add(chromoListReducedFile);
                rpanelReducedFilteredFile.add(chromoListReducedFilteredFile);
                rpanelReducedCondensedFile.add(chromoListReducedCondensedFile);

                // Only for all the panels combined
                if (j == 0) {
                    rpanelCombinedReducedFilteredFile.add(chromoListCombinedReducedFilteredFile);
                    rpanelCombinedReducedCondensedFile.add(chromoListCombinedReducedCondensedFile);
                    rpanelCombinedReducedFilteredXFile.add(chromoListCombinedReducedFilteredXFile);
                    // rpanelCombinedReducedCondensedXFile.add(chromoListCombinedReducedCondensedXFile);
                }

                rpanelFilteredByAllFile.add(chromoFilteredByAllFile);
                rpanelCondensedFile.add(chromoCondensedFile);

                // Here we have to create an additional list of condensed files that will be used when we execute
                // jointCondensedFiles Task,
                // for all chromosomes.
                // The number of additional files is the number of chromosomes minus 1.
                int addCondensed = 0;
                for (int deep = startChr; deep < endChr; deep++) {
                    String tmpAdditionalCondensedFileName = null;
                    // if(startChr == endChr) {
                    // tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_chr_" + startChr +
                    // ".txt.gz";
                    // } else
                    if (deep == (endChr - 1)) {
                        tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_chr_" + startChr + "_to_" + endChr
                                + ".txt.gz";
                    } else {
                        tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_" + addCondensed + ".txt.gz";
                    }
                    String tmpAdditionalCondensedFile = rpanelOutDirSummary + "/" + tmpAdditionalCondensedFileName;
                    GenericFile myAdditionalCondensedFile = new GenericFile(rpanelOutDirSummary, tmpAdditionalCondensedFileName,
                            "uncompressed", "none");
                    additionalCondensedFile.add(myAdditionalCondensedFile);
                    // System.out.println("\t[MergeFiles.java] " + tmpAdditionalCondensedFile);

                    addCondensed++;
                }

                if (startChr == endChr) {
                    String tmpAdditionalCondensedFileName = testTypeName + "_" + rPanel + "_condensed_chr_" + startChr + ".txt.gz";
                    String tmpAdditionalCondensedFile = rpanelOutDirSummary + "/" + tmpAdditionalCondensedFileName;
                    GenericFile myAdditionalCondensedFile = new GenericFile(rpanelOutDirSummary, tmpAdditionalCondensedFileName,
                            "uncompressed", "none");
                    additionalCondensedFile.add(myAdditionalCondensedFile);
                    // System.out.println("\t[MergeFiles.java] only " + tmpAdditionalCondensedFile);
                    addCondensed++;
                }

                rpanelAdditionalCondensedIndex.add(addCondensed);

                // Here we have to create an additional list of filteredByAll files that will be used when we execute
                // jointFilteredByAllFile task
                // for all chromosomes.
                // Unlike the previous case with condensed files, we can not include chromosome 23. (Chr 23 format fo
                // filteredByAllFile is different to
                // the rest of chromosomes (thanks to snptest).

                // The number of additional files is the number of chromosomes minus 1.
                int addFiltered = 0;
                int endChrNormal = endChr;
                if (startChr < endChr) {
                    if (endChr != 23) {
                        endChrNormal = endChr;
                    } else {
                        endChrNormal = endChr - 1;
                    }
                }

                for (int deep = startChr; deep < endChrNormal; deep++) {
                    String tmpAdditionalFilteredByAllFileName = null;
                    // if(startChr == endChrNormal) {
                    // tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + "_filteredByAll_chr_" +
                    // startChr + ".txt.gz";
                    // } else
                    if (deep == (endChrNormal - 1)) {
                        tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + "_filteredByAll_chr_" + startChr + "_to_"
                                + endChrNormal + ".txt.gz";
                    } else {
                        tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + "_filteredByAll_" + addFiltered + ".txt.gz";
                    }
                    String tmpAdditionalFilteredByAllFile = rpanelOutDirSummary + "/" + tmpAdditionalFilteredByAllFileName;
                    GenericFile myAdditionalFilteredByAllFile = new GenericFile(rpanelOutDirSummary, tmpAdditionalFilteredByAllFileName,
                            "uncompressed", "none");
                    additionalFilteredByAllFile.add(myAdditionalFilteredByAllFile);
                    // System.out.println("\t[MergeFiles.java] " + tmpAdditionalFilteredByAllFile);

                    addFiltered++;
                }

                if (startChr == endChrNormal) {
                    String tmpAdditionalFilteredByAllFileName = testTypeName + "_" + rPanel + "_filteredByAll_chr_" + startChr + ".txt.gz";
                    String tmpAdditionalFilteredByAllFile = rpanelOutDirSummary + "/" + tmpAdditionalFilteredByAllFileName;
                    GenericFile myAdditionalFilteredByAllFile = new GenericFile(rpanelOutDirSummary, tmpAdditionalFilteredByAllFileName,
                            "uncompressed", "none");
                    additionalFilteredByAllFile.add(myAdditionalFilteredByAllFile);
                    // System.out.println("\t[MergeFiles.java] only " + tmpAdditionalFilteredByAllFile);
                    addFiltered++;
                }

                rpanelAdditionalFilteredByAllIndex.add(addFiltered);

                rpanelAdditionalCondensedFile.add(additionalCondensedFile);
                rpanelAdditionalFilteredByAllFile.add(additionalFilteredByAllFile);

                // If there is chr 23:
                if (endChr == 23) {
                    String tmpAdditionalFilteredByAllXFileName = testTypeName + "_" + rPanel + "_filteredByAll_chr_" + endChr + ".txt.gz";
                    String tmpAdditionalFilteredByAllXFile = rpanelOutDirSummary + "/" + tmpAdditionalFilteredByAllXFileName;
                    GenericFile myAdditionalFilteredByAllXFile = new GenericFile(rpanelOutDirSummary, tmpAdditionalFilteredByAllXFileName,
                            "uncompressed", "none");

                    additionalFilteredByAllXFile.add(myAdditionalFilteredByAllXFile);
                    // System.out.println("\t[MergeFiles.java] " + tmpAdditionalFilteredByAllXFile);

                    rpanelAdditionalFilteredByAllXFile.add(additionalFilteredByAllXFile);
                }
            } // End of for(j=0; j<refPanels.size();j++)
              // Now we have to build the list of reduced files for the type of Test. We store this list
            testTypeReducedOutDir.add(rpanelListOutDir);
            testTypeReducedFile.add(rpanelReducedFile);
            testTypeReducedFilteredFile.add(rpanelReducedFilteredFile);
            testTypeReducedCondensedFile.add(rpanelReducedCondensedFile);

            testTypeCombinedReducedFilteredFile.add(rpanelCombinedReducedFilteredFile);
            testTypeCombinedReducedCondensedFile.add(rpanelCombinedReducedCondensedFile);
            testTypeCombinedReducedFilteredXFile.add(rpanelCombinedReducedFilteredXFile);
            // testTypeCombinedReducedCondensedXFile.add(rpanelCombinedReducedCondensedXFile);

            testTypeFilteredByAllFile.add(rpanelFilteredByAllFile);
            testTypeCondensedFile.add(rpanelCondensedFile);

            testTypeAdditionalCondensedFile.add(rpanelAdditionalCondensedFile);
            testTypeAdditionalFilteredByAllFile.add(rpanelAdditionalFilteredByAllFile);

            if (endChr == 23) {
                testTypeAdditionalFilteredByAllXFile.add(rpanelAdditionalFilteredByAllXFile);
            }

            testTypeAdditionalCondensedIndex.add(rpanelAdditionalCondensedIndex);
            testTypeAdditionalFilteredByAllIndex.add(rpanelAdditionalFilteredByAllIndex);

        } // end of for(int tt=0; tt< numberOfTestTypesNames; tt++)
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

        int i = chromo - startChr;
        // int i= chromo - 1 ;
        return testTypeReducedOutDir.get(testTypeIndex).get(rPanelIndex).get(i);
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFileName is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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
    public void setReducedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFilteredFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFilteredFileName is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFilteredFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFilteredFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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
    public void setReducedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFilteredFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFilteredFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFilteredFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFilteredFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedCondensedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedCondensedFileName is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedCondensedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedCondensedFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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
    public void setReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedCondensedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedCondensedFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedCondensedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedCondensedFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeCombinedReducedFilteredFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredFileName is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedFilteredFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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
    public void setCombinedReducedFilteredFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedFilteredFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedFilteredFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredFile is greater than the
        // existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedFilteredFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
    }

    /**
     * Method to access reducedCombinedReducedFilteredXFileName
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param index
     * @return
     */
    public String getCombinedReducedFilteredXFileName(int testTypeIndex, int rPanelIndex, int chromo, int index) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeCombinedReducedFilteredXFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredXFileName is greater
        // than the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedFilteredXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to access reducedCombinedReducedFilteredXFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param index
     * @return
     */
    public String getCombinedReducedFilteredXFile(int testTypeIndex, int rPanelIndex, int chromo, int index) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedFilteredXFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredXFile is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedFilteredXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the reducedCombinedReducedFilteredXFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param index
     * @param finalStatus
     */
    public void setCombinedReducedFilteredXFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedFilteredXFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredXFile is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeCombinedReducedFilteredXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the reducedCombinedReducedFilteredXFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param chromo
     * @param index
     * @return
     */
    public String getCombinedReducedFilteredXFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedFilteredXFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedFilteredXFile is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedFilteredXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeCombinedReducedCondensedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedFileName is greater
        // than the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedCondensedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedFile is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
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
    public void setCombinedReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index, String finalStatus) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedCondensedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedFile is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(finalStatus);
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
    public String getCombinedReducedCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int chromo, int index) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromo);

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeCombinedReducedCondensedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedFile is greater than
        // the existing in chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCombinedReducedCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus();
    }

    /*
     * // Method to access reducedCombinedReducedCondensedXFileName public String
     * getCombinedReducedCondensedXFileName(int testTypeIndex, int rPanelIndex,int chromo, int index) {
     * 
     * if ( (chromo < 1) || (chromo > maxNumberOfChromosomes) ) { System.err.println("[MergeFiles] Error, chromosome " +
     * chromo + "does not exist"); System.exit(1); } //int i = chromo -1; int i = chromo - startChr;
     * 
     * //TODO: IMPORTANT: Verify the index range!!!! //ArrayList<String> tmpList = new ArrayList<String>();
     * //testTypeCombinedReducedCondensedXFileName.get(rPanelIndex).get(i); //int lastIndex = tmpList.size() - 1;
     * 
     * //if(index > lastIndex) { // System.err.
     * println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedXFileName is greater than the existing in chromosome "
     * + chromo); // System.err.println("             index " + index + " > lastIndex = " + lastIndex); //
     * System.exit(1); //}
     * 
     * return testTypeCombinedReducedCondensedXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
     * }
     * 
     * 
     * // Method to access reducedCombinedReducedCondensedXFile public String getCombinedReducedCondensedXFile(int
     * testTypeIndex, int rPanelIndex,int chromo, int index) { if ( (chromo < 1) || (chromo > maxNumberOfChromosomes) )
     * { System.err.println("[MergeFiles] Error, chromosome " + chromo + "does not exist"); System.exit(1); } //int i =
     * chromo -1; int i = chromo - startChr; //TODO: IMPORTANT: Verify the index range!!!! //ArrayList<String> tmpList =
     * new ArrayList<String>(); //tmpList = (testTypeCombinedReducedCondensedXFileName.get(rPanelIndex, i)); //int
     * lastIndex = tmpList.size() - 1;
     * 
     * //if(index > lastIndex) { // System.err.
     * println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedXFile  is greater than the existing in chromosome "
     * + chromo); // System.err.println("             index " + index + " > lastIndex = " + lastIndex); //
     * System.exit(1); //}
     * 
     * return testTypeCombinedReducedCondensedXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFullName();
     * }
     * 
     * /** Method to set the finalStatus of the reducedCombinedReducedCondensedXFile / public void
     * setCombinedReducedCondensedXFileFinalStatus(int testTypeIndex, int rPanelIndex,int chromo,int index, String
     * finalStatus) { if ( (chromo < 1) || (chromo > maxNumberOfChromosomes) ) {
     * System.err.println("[MergeFiles] Error, chromosome " + chromo + "does not exist"); System.exit(1); } //int i =
     * chromo -1; int i = chromo - startChr; //TODO: IMPORTANT: Verify the index range!!!! //ArrayList<String> tmpList =
     * new ArrayList<String>(); //tmpList = (testTypeCombinedReducedCondensedXFileName.get(rPanelIndex, i)); //int
     * lastIndex = tmpList.size() - 1;
     * 
     * //if(index > lastIndex) { // System.err.
     * println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedXFile  is greater than the existing in chromosome "
     * + chromo); // System.err.println("             index " + index + " > lastIndex = " + lastIndex); //
     * System.exit(1); //}
     * 
     * testTypeCombinedReducedCondensedXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).setFinalStatus(
     * finalStatus); }
     * 
     * 
     * /** Method to access the finalStatus of the reducedCombinedReducedCondensedXFile / public String
     * getCombinedReducedCondensedXFileFinalStatus(int testTypeIndex, int rPanelIndex,int chromo,int index) { if (
     * (chromo < 1) || (chromo > maxNumberOfChromosomes) ) { System.err.println("[MergeFiles] Error, chromosome " +
     * chromo + "does not exist"); System.exit(1); } //int i = chromo -1; int i = chromo - startChr; //TODO: IMPORTANT:
     * Verify the index range!!!! //ArrayList<String> tmpList = new ArrayList<String>(); //tmpList =
     * (testTypeCombinedReducedCondensedXFileName.get(rPanelIndex, i)); //int lastIndex = tmpList.size() - 1;
     * 
     * //if(index > lastIndex) { // System.err.
     * println("[MergeFiles] Error, the number of testTypeCombinedReducedCondensedXFile  is greater than the existing in chromosome "
     * + chromo); // System.err.println("             index " + index + " > lastIndex = " + lastIndex); //
     * System.exit(1); //}
     * 
     * return
     * testTypeCombinedReducedCondensedXFile.get(testTypeIndex).get(rPanelIndex).get(i).get(index).getFinalStatus(); }
     */

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

        // int i = chromo -1;
        int i = chromo - startChr;

        // ArrayList<String> tmpList = new ArrayList<String>();
        int lastIndex = testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size() - 1;
        // System.out.println("[MergeFiles] lastIndex size = " + lastIndex);

        return testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).get(lastIndex).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // ArrayList<String> tmpList = new ArrayList<String>();
        int lastIndex = testTypeReducedFile.get(testTypeIndex).get(rPanelIndex).get(i).size() - 1;
        // System.out.println("[MergeFiles] lastIndex size = " + lastIndex);
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
        // System.out.println("\t[MergeFiles] theLastReducedFiles ARE:");
        for (int hh = 0; hh <= testTypeIndex; hh++) {
            for (int kk = 0; kk <= rPanelIndex; kk++) {
                for (int j = 0; j <= chromo - startChr; j++) {
                    // System.out.println("\t[MergeFiles] " + testTypeReducedFile.get(rPanelIndex).get(kk).get(j));
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFileName is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).setFinalStatus(finalStatus);
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(i).getFinalStatus();
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

        // int i = chromo -1;
        int i = chromo - startChr;

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFileName is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFullName();
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).setFinalStatus(finalStatus);
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

        // int i = chromo -1;
        int i = chromo - startChr;
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeCondensedFile.get(testTypeIndex).get(rPanelIndex).get(i).getFinalStatus();
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

        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFileName is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the additionalCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param index
     * @param finalStatus
     */
    public void setAdditionalCondensedFileFinalStatus(int testTypeIndex, int rPanelIndex, int index, String finalStatus) {
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).setFinalStatus(finalStatus);
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(index).getFinalStatus();
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFileName is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the additionalCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param index
     * @param finalStatus
     */
    public void setAdditionalFilteredByAllFileFinalStatus(int testTypeIndex, int rPanelIndex, int index, String finalStatus) {
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).setFinalStatus(finalStatus);
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(index).getFinalStatus();
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // testTypeReducedFileName.get(rPanelIndex).get(i);
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFileName is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the additionalCondensedFile
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @param index
     * @param finalStatus
     */
    public void setAdditionalFilteredByAllXFileFinalStatus(int testTypeIndex, int rPanelIndex, int index, String finalStatus) {
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index).setFinalStatus(finalStatus);
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
        // TODO: IMPORTANT: Verify the index range!!!!
        // ArrayList<String> tmpList = new ArrayList<String>();
        // tmpList = (testTypeReducedFileName.get(rPanelIndex, i));
        // int lastIndex = tmpList.size() - 1;

        // if(index > lastIndex) {
        // System.err.println("[MergeFiles] Error, the number of testTypeReducedFile is greater than the existing in
        // chromosome " + chromo);
        // System.err.println(" index " + index + " > lastIndex = " + lastIndex);
        // System.exit(1);
        // }

        return testTypeAdditionalFilteredByAllXFile.get(testTypeIndex).get(rPanelIndex).get(index).getFinalStatus();
    }

    /**
     * Method to access the last CondensedFile of each testType and rPanel
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getFinalCondensedFile(int testTypeIndex, int rPanelIndex) {
        int lastIndex = testTypeAdditionalCondensedIndex.get(testTypeIndex).get(rPanelIndex);
        return testTypeAdditionalCondensedFile.get(testTypeIndex).get(rPanelIndex).get(lastIndex - 1).getFullName();
    }

    /**
     * Method to access the last FilteredByAllFile of each testType and rPanel
     * 
     * @param testTypeIndex
     * @param rPanelIndex
     * @return
     */
    public String getFinalFilteredByAllFile(int testTypeIndex, int rPanelIndex) {
        int lastIndex = testTypeAdditionalFilteredByAllIndex.get(testTypeIndex).get(rPanelIndex);
        return testTypeAdditionalFilteredByAllFile.get(testTypeIndex).get(rPanelIndex).get(lastIndex - 1).getFullName();
    }

    private void checkChromoIndex(int chromo) {
        if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
            System.err.println("[MergeFiles] Error, chromosome " + chromo + "does not exist");
            System.exit(1);
        }
    }

}
