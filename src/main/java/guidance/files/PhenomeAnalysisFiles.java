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

import java.util.ArrayList;


public class PhenomeAnalysisFiles {

    private String phenomeAnalysisOutDir = null;

    private ArrayList<GenericFile> phenomeAnalysisFile = new ArrayList<>();
    private ArrayList<GenericFile> phenomeTmpFile = new ArrayList<>();
    private ArrayList<GenericFile> phenomeFinalFile = new ArrayList<>();


    /**
     * New PhenomeAnalysisFiles instance
     * 
     * @param parsingArgs
     * @param baseOutDir
     * @param refPanels
     */
    public PhenomeAnalysisFiles(ParseCmdLine parsingArgs, String baseOutDir, List<String> refPanels) {
        String testTypeName = null;
        String rPanel = null;
        int numberOfTestTypesNames = parsingArgs.getNumberOfTestTypeName();
        int numberOfRpanels = refPanels.size();

        int totalCombinations = numberOfTestTypesNames * numberOfRpanels - 1;

        // We create the first directory name: the cohort directory
        String mixedCohort = parsingArgs.getCohort();

        phenomeAnalysisOutDir = baseOutDir + "/associations/pheno_analysis/" + mixedCohort;
        // for(int tt=0; tt< numberOfTestTypesNames; tt++) {
        // testTypeName = parsingArgs.getTestTypeName(tt);
        // phenomeAnalysisOutDir = phenomeAnalysisOutDir + "_" + testTypeName;
        // }

        String basePrefixPhenomeName = "phenome_analysis_" + mixedCohort;
        String prefixPhenomeName = basePrefixPhenomeName;
        String prefixPhenomeTmpName = basePrefixPhenomeName;

        // Phantom name as the first position of phenomeFinalFile
        int counter = 0;
        String prefixPhenomeFinalName = prefixPhenomeTmpName + "_" + counter + "_final";
        String tmpPhenomeFinalFileName = prefixPhenomeFinalName + ".txt.gz";
        String tmpPhenomeFinalFile = phenomeAnalysisOutDir + "/" + tmpPhenomeFinalFileName;
        GenericFile myPhenomeFinalFile = new GenericFile(phenomeAnalysisOutDir, tmpPhenomeFinalFileName, "uncompressed", "none");
        phenomeFinalFile.add(myPhenomeFinalFile);

        counter = 0;
        for (int tt = 0; tt < numberOfTestTypesNames; tt++) {
            testTypeName = parsingArgs.getTestTypeName(tt);
            // prefixPhenomeName = basePrefixPhenomeName + "_" + testTypeName;

            int startj = 0;
            /*
             * if(tt==0) { rPanel = refPanels.get(0); prefixPhenomeName = prefixPhenomeName + "_" + rPanel; startj=1; }
             */
            for (int j = startj; j < refPanels.size(); j++) {
                rPanel = refPanels.get(j);

                String tmpPhenomeTmpFileName = prefixPhenomeTmpName + "_" + counter + ".txt.gz";
                String tmpPhenomeTmpFile = phenomeAnalysisOutDir + "/" + tmpPhenomeTmpFileName;
                GenericFile myPhenomeTmpFile = new GenericFile(phenomeAnalysisOutDir, tmpPhenomeTmpFileName, "uncompressed", "none");
                phenomeTmpFile.add(myPhenomeTmpFile);

                prefixPhenomeName = basePrefixPhenomeName + "_" + testTypeName + "_" + rPanel;
                String tmpPhenomeAnalysisFileName = prefixPhenomeName + ".txt.gz";
                String tmpPhenomeAnalysisFile = phenomeAnalysisOutDir + "/" + tmpPhenomeAnalysisFileName;
                GenericFile myPhenomeAnalysisFile = new GenericFile(phenomeAnalysisOutDir, tmpPhenomeAnalysisFileName, "uncompressed",
                        "none");
                phenomeAnalysisFile.add(myPhenomeAnalysisFile);

                prefixPhenomeFinalName = prefixPhenomeTmpName + "_" + counter + "_final";
                tmpPhenomeFinalFileName = prefixPhenomeFinalName + ".txt.gz";
                tmpPhenomeFinalFile = phenomeAnalysisOutDir + "/" + tmpPhenomeFinalFileName;
                myPhenomeFinalFile = new GenericFile(phenomeAnalysisOutDir, tmpPhenomeFinalFileName, "uncompressed", "none");
                phenomeFinalFile.add(myPhenomeFinalFile);

                // System.out.println("[PhenomeAnalysis] " + phenomeAnalysisOutDir + "/" + tmpPhenomeAnalysisFileName);
                counter++;

            } // End of for(j=1; j<refPanels.size();j++)
              // Now we have to build the list of reduced files for the type of Test. We store this list

        } // end of for(int tt=0; tt< numberOfTestTypesNames; tt++)
    }

    /**
     * Method to access correctedPvaluesFile information
     * 
     * @return
     */
    public String getPhenomeAnalysisOutDir() {
        return phenomeAnalysisOutDir;
    }

    /**
     * Method to access phenomeAnalysisFileName
     * 
     * @param index
     * @return
     */
    public String getPhenotypeIntermediateFileName(int index) {
        return phenomeTmpFile.get(index).getName();
    }

    /**
     * Method to access phenomeAnalysisFile
     * 
     * @param index
     * @return
     */
    public String getPhenotypeIntermediateFile(int index) {
        return phenomeTmpFile.get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the phenomeAnalysisFile
     * 
     * @param index
     * @param finalStatus
     */
    public void setPhenotypeIntermediateFileFinalStatus(int index, String finalStatus) {
        phenomeTmpFile.get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the phenomeAnalysisFile
     * 
     * @param index
     * @return
     */
    public String getPhenotypeIntermediateFileFinalStatus(int index) {
        return phenomeTmpFile.get(index).getFinalStatus();
    }

    /**
     * Method to access phenomeAnalysisFileName
     * 
     * @param index
     * @return
     */
    public String getPhenotypeFileName(int index) {
        return phenomeAnalysisFile.get(index).getName();
    }

    /**
     * Method to access phenomeAnalysisFile
     * 
     * @param index
     * @return
     */
    public String getPhenotypeFile(int index) {
        return phenomeAnalysisFile.get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the phenomeAnalysisFile
     * 
     * @param index
     * @param finalStatus
     */
    public void setPhenotypeFileFinalStatus(int index, String finalStatus) {
        phenomeAnalysisFile.get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the phenomeAnalysisFile
     * 
     * @param index
     * @return
     */
    public String getPhenotypeFileFinalStatus(int index) {
        return phenomeAnalysisFile.get(index).getFinalStatus();
    }

    /**
     * Method to access phenomeFinalFileName
     * 
     * @param index
     * @return
     */
    public String getPhenotypeFinalFileName(int index) {
        return phenomeFinalFile.get(index).getName();
    }

    /**
     * Method to access phenomeFinalFile
     * 
     * @param index
     * @return
     */
    public String getPhenotypeFinalFile(int index) {
        return phenomeFinalFile.get(index).getFullName();
    }

    /**
     * Method to set the finalStatus of the phenomeFinalFile
     * 
     * @param index
     * @param finalStatus
     */
    public void setPhenotypeFinalFileFinalStatus(int index, String finalStatus) {
        phenomeFinalFile.get(index).setFinalStatus(finalStatus);
    }

    /**
     * Method to access the finalStatus of the phenomeFinalFile
     * 
     * @param index
     * @return
     */
    public String getPhenotypeFinalFileFinalStatus(int index) {
        return phenomeFinalFile.get(index).getFinalStatus();
    }

}
