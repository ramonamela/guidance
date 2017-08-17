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

package guidance.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.File;
import java.util.HashMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;


public class ParseCmdLine {

    // Logger
    private static final Logger LOGGER = LogManager.getLogger("Console");

    // Error messages
    private static final String CLASS_HEADER = "[ParseCmdLine.java]";
    private static final String ERROR_PARAM_ORDER = CLASS_HEADER + " Error in the order of parameters, in parameter: ";
    private static final String ERROR_SYNTAX = CLASS_HEADER + " Error of sintax in ";
    private static final String ERROR_SYNTAX_SUFFIX = ", in parameter: ";

    // Chromosome sizes
    private static final int MINIMUMCHUNKSIZE = 1_000;
    private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;

    private String gwasConfigFile = null;
    private ArrayList<String> argumentsArray = new ArrayList<>();

    private String mixedCohort = null;
    private String mixedBedDir = null;
    private String mixedFamFileName = null;
    private String mixedBimFileName = null;
    private String mixedBedFileName = null;

    private String mixedChrDir = null;
    private ArrayList<String> mixedGenFileName = new ArrayList<>();
    private ArrayList<String> mixedGenFile = new ArrayList<>();

    private String mixedSampleDir = null;
    private String mixedSampleFileName = null;
    private String mixedSampleFile = null;

    private String gmapDir = null;
    private ArrayList<String> gmapFileName = new ArrayList<>();

    private String exclCgatSnp = null;
    private String exclSVSnp = "NO";
    private String imputationTool = null;

    // chunkSize is always rewritten below.
    private int chunkSize = 1_000_000;

    private int refPanelNumber = 0;
    // By default we do not combine panels
    private boolean refPanelCombine = false;
    private ArrayList<String> rpanelTypes = new ArrayList<>();
    private ArrayList<String> rpanelDir = new ArrayList<>();

    private ArrayList<ArrayList<String>> rpanelHapFileName = new ArrayList<>();
    private ArrayList<ArrayList<String>> rpanelLegFileName = new ArrayList<>();

    // testTypes will be organized as follows:
    // Each string will have: ["test_type_name","response_variable", "covariables"]
    // covariables will be a string like: "se,sex,bmi". That is values separated by coloms.
    private ArrayList<String> testTypesNames = new ArrayList<>();
    private ArrayList<String> responseVars = new ArrayList<>();
    private ArrayList<String> covariables = new ArrayList<>();

    // String rpanelTypes = "";
    private String outDir = null;
    private int start = 0;
    private int end = 0;
    private int endNormal = 0;

    private String wfDeepRequired = null;
    private HashMap<String, Integer> wfPossibleDeeps = new HashMap<>();
    private HashMap<String, Integer> wfAllStages = new HashMap<>();

    private Double mafThreshold = 0.000;
    private Double infoThreshold = 0.000;
    private Double hweCohortThreshold = 0.000;
    private Double hweCasesThreshold = 0.000;
    private Double hweControlsThreshold = 0.000;

    private String listOfStagesFile = "list_of_stages_default.txt";
    private String removeTemporalFiles = "NO";
    private String compressFiles = "NO";

    private String inputFormat = null;


    /**
     * Parse CMD Args into internal values
     * 
     * @param args
     */
    public ParseCmdLine(String[] args) {
        final Integer EMPTY_MASK = 0x0000000;
        wfPossibleDeeps.put("until_convertFromBedToBed", EMPTY_MASK);
        wfPossibleDeeps.put("until_phasing", EMPTY_MASK);
        wfPossibleDeeps.put("until_imputation", EMPTY_MASK);
        wfPossibleDeeps.put("until_association", EMPTY_MASK);
        wfPossibleDeeps.put("until_filterByAll", EMPTY_MASK);
        wfPossibleDeeps.put("until_summary", EMPTY_MASK);
        wfPossibleDeeps.put("whole_workflow", EMPTY_MASK);
        wfPossibleDeeps.put("from_phasing", EMPTY_MASK);
        wfPossibleDeeps.put("from_phasing_to_summary", EMPTY_MASK);
        wfPossibleDeeps.put("from_phasing_to_filterByAll", EMPTY_MASK);
        wfPossibleDeeps.put("from_phasing_to_association", EMPTY_MASK);
        wfPossibleDeeps.put("from_phasing_to_imputation", EMPTY_MASK);
        wfPossibleDeeps.put("from_imputation", EMPTY_MASK);
        wfPossibleDeeps.put("from_imputation_to_summary", EMPTY_MASK);
        wfPossibleDeeps.put("from_imputation_to_filterByAll", EMPTY_MASK);
        wfPossibleDeeps.put("from_imputation_to_association", EMPTY_MASK);
        wfPossibleDeeps.put("from_imputation_to_filterByInfo", EMPTY_MASK);
        wfPossibleDeeps.put("from_filterByInfo_to_qctoolS", EMPTY_MASK);
        wfPossibleDeeps.put("from_qctoolS_to_association", EMPTY_MASK);
        wfPossibleDeeps.put("from_association", EMPTY_MASK);
        wfPossibleDeeps.put("from_association_to_filterByAll", EMPTY_MASK);
        wfPossibleDeeps.put("from_association_to_summary", EMPTY_MASK);
        wfPossibleDeeps.put("from_filterByAll", EMPTY_MASK);
        wfPossibleDeeps.put("from_filterByAll_to_summary", EMPTY_MASK);
        wfPossibleDeeps.put("from_summary", EMPTY_MASK);

        // Step 1: We read the file with the configuration and
        // clean the lines from (spaces and comments
        // Then, we put the parameters in an array.
        if ((args.length < 1) && (args.length > 2)) {
            LOGGER.fatal(CLASS_HEADER + " Error in the command line parameters.");
            LOGGER.fatal(CLASS_HEADER + " Usage -1-: -config_file configuration_file.txt");
            System.exit(1);
        }

        if (args[0].equals("-config_file")) {
            gwasConfigFile = args[1];
        }

        try (FileReader fr = new FileReader(gwasConfigFile); BufferedReader br = new BufferedReader(fr);) {
            String line = null;
            while ((line = br.readLine()) != null) {
                char firstChar = line.charAt(0);
                if (firstChar != '#') {
                    String myLine = line.replaceAll(" ", "");
                    myLine = myLine.replaceAll("\t", "");
                    argumentsArray.add(myLine);
                }
                // Process the line
            }
        } catch (IOException ioe) {
            LOGGER.fatal(CLASS_HEADER + " Error opening/reading " + gwasConfigFile);
            LOGGER.fatal(ioe.getMessage());
        }

        // Now, we load the parameters of the execution.
        // There is a strict order in which parameters should be put in the configuration input file.
        // We follow this order.
        int i = 0;
        String tmpArg = argumentsArray.get(i++);
        String[] myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("wfDeep")) {
                wfDeepRequired = myArgument[1];
                boolean validKey = wfPossibleDeeps.containsKey(wfDeepRequired);
                if (!validKey) {
                    LOGGER.fatal(CLASS_HEADER + " Error, wfDeep parameter " + wfDeepRequired + " is not accepted");
                    LOGGER.fatal(CLASS_HEADER + "        The only accepted values are:");
                    for (String key : wfPossibleDeeps.keySet()) {
                        LOGGER.fatal(CLASS_HEADER + "- " + key);
                    }
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + i);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + i);
            System.exit(1);
        }

        int tmpStart = 0;
        int tmpEnd = 0;
        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("init_chromosome")) {
                tmpStart = Integer.parseInt(myArgument[1]);
                if (tmpStart < 1 || tmpStart > MAX_NUMBER_OF_CHROMOSOMES) {
                    LOGGER.fatal(CLASS_HEADER + " Error, init_chromosome = " + tmpStart);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: should be: 0<init_chromosome<=23");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("end_chromosome")) {
                tmpEnd = Integer.parseInt(myArgument[1]);
                if (tmpEnd < 1 || tmpEnd > MAX_NUMBER_OF_CHROMOSOMES || tmpEnd < tmpStart) {
                    LOGGER.fatal(CLASS_HEADER + " Error, end_chromosome = " + tmpEnd);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: should be: 0<init_chromosome<23 and >= init_chromosome");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        // Here we have to analyse if we have to include the chromosome 23 or not. If we have to include it, then we
        // enable the doChr23 variable.
        start = tmpStart;
        end = tmpEnd;
        for (int counter = start; counter <= end; counter++) {
            if (counter < 23) {
                endNormal = counter;
            }
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("maf_threshold")) {
                mafThreshold = Double.parseDouble(myArgument[1]);
                if (mafThreshold < 0) {
                    LOGGER.fatal(CLASS_HEADER + " Error, maf_threshold = " + mafThreshold);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: should be: > 0");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("info_threshold")) {
                infoThreshold = Double.parseDouble(myArgument[1]);
                if (infoThreshold < 0) {
                    LOGGER.fatal(CLASS_HEADER + " Error, info_threshold = " + infoThreshold);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: should be: > 0");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("hwe_cohort_threshold")) {
                hweCohortThreshold = Double.parseDouble(myArgument[1]);
                // if( hweCohortThreshold < 0 ) {
                // LOGGER.fatal(CLASS_HEADER + " Error, hwe_cohort_threshold = " + hweCohortThreshold);
                // LOGGER.fatal(" It should be: should be: > 0");
                // System.exit(1);
                // }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("hwe_cases_threshold")) {
                hweCasesThreshold = Double.parseDouble(myArgument[1]);
                // if( hweCasesThreshold < 0 ) {
                // LOGGER.fatal("Error, hwe_cases_threshold = " + hweCasesThreshold);
                // LOGGER.fatal(" It should be: should be: > 0");
                // System.exit(1);
                // }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("hwe_controls_threshold")) {
                hweControlsThreshold = Double.parseDouble(myArgument[1]);
                // if( hweControlsThreshold < 0 ) {
                // LOGGER.fatal("Error, hwe_controls_threshold = " + hweControlsThreshold);
                // LOGGER.fatal(" It should be: should be: > 0");
                // System.exit(1);
                // }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("exclude_cgat_snps")) {
                exclCgatSnp = myArgument[1].toUpperCase();
                if (!exclCgatSnp.equals("YES") && !exclCgatSnp.equals("NO")) {
                    LOGGER.fatal(CLASS_HEADER + " Error, exclude_cgat_snps = " + exclCgatSnp);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: should be: YES or NOT");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("imputation_tool")) {
                imputationTool = myArgument[1];
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("test_types")) {
                String tmpTestTypes = myArgument[1];
                String[] tmpTestTypesArray = tmpTestTypes.split(",");
                for (int kk = 0; kk < tmpTestTypesArray.length; kk++) {
                    testTypesNames.add(tmpTestTypesArray[kk]);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        for (int kk = 0; kk < testTypesNames.size(); kk++) {
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals(testTypesNames.get(kk))) {
                    String[] tmpFields = myArgument[1].split(":");
                    responseVars.add(tmpFields[0]);
                    covariables.add(tmpFields[1]);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("chunk_size_analysis")) {
                chunkSize = Integer.parseInt(myArgument[1]);
                if (chunkSize < MINIMUMCHUNKSIZE) {
                    LOGGER.fatal(CLASS_HEADER + " Error, the value for chunk_size_analysis parameter should not be less than "
                            + MINIMUMCHUNKSIZE);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("file_name_for_list_of_stages")) {
                listOfStagesFile = myArgument[1];
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("remove_temporal_files")) {
                removeTemporalFiles = myArgument[1];
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("compress_files")) {
                compressFiles = myArgument[1];
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("input_format")) {
                inputFormat = myArgument[1].toUpperCase();
                if (!inputFormat.equals("GEN") && !inputFormat.equals("BED")) {
                    LOGGER.fatal(CLASS_HEADER + " Error, input_format is incorrect = " + inputFormat);
                    LOGGER.fatal(CLASS_HEADER + "        It should be GEN or BED");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("mixed_cohort")) {
                mixedCohort = myArgument[1];
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        if (inputFormat.equals("GEN")) {
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("mixed_gen_file_dir")) {
                    mixedChrDir = myArgument[1];
                    checkExistence(mixedChrDir);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }

            for (int kk = start; kk <= end; kk++) {
                tmpArg = argumentsArray.get(i++);
                myArgument = tmpArg.split("=");
                if ((myArgument.length > 0) && (myArgument.length < 3)) {
                    String chromo = Integer.toString(kk);
                    String tmpfile = "mixed_gen_file_chr_" + chromo;
                    if (myArgument[0].equals(tmpfile)) {
                        mixedGenFileName.add(myArgument[1]);
                        checkExistence(mixedChrDir + File.separator + myArgument[1]);
                    } else {
                        LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                        System.exit(1);
                    }
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                    System.exit(1);
                }
            }
        } else if (inputFormat.equals("BED")) {
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("mixed_bed_file_dir")) {
                    mixedBedDir = myArgument[1];
                    checkExistence(mixedBedDir);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }

            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("mixed_bed_file")) {
                    mixedBedFileName = myArgument[1];
                    checkExistence(mixedBedDir + "/" + mixedBedFileName);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }

            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("mixed_bim_file")) {
                    mixedBimFileName = myArgument[1];
                    checkExistence(mixedBedDir + "/" + mixedBimFileName);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("mixed_fam_file")) {
                    mixedFamFileName = myArgument[1];
                    checkExistence(mixedBedDir + "/" + mixedFamFileName);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + " Error of in " + gwasConfigFile + ", the " + inputFormat + " input format is not supported");
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("mixed_sample_file_dir")) {
                mixedSampleDir = myArgument[1];
                checkExistence(mixedSampleDir);
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            String tmpfile = "mixed_sample_file";
            if (myArgument[0].equals(tmpfile)) {
                mixedSampleFileName = myArgument[1];
                checkExistence(mixedSampleDir + File.separator + mixedSampleFileName);
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("genmap_file_dir")) {
                gmapDir = myArgument[1];
                checkExistence(gmapDir);
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        for (int kk = start; kk <= end; kk++) {
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                String chromo = Integer.toString(kk);
                String tmpfile = "genmap_file_chr_" + chromo;
                if (myArgument[0].equals(tmpfile)) {
                    gmapFileName.add(myArgument[1]);
                    checkExistence(gmapDir + File.separator + myArgument[1]);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("refpanel_number")) {
                refPanelNumber = Integer.parseInt(myArgument[1]);
                if (refPanelNumber < 1) {
                    LOGGER.fatal(CLASS_HEADER + " Error, refpanel_number = " + refPanelNumber);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: should be: >0");
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("refpanel_combine")) {
                String myTempCombine = myArgument[1];
                myTempCombine = myTempCombine.toUpperCase();
                if (!myTempCombine.equals("YES") && !myTempCombine.equals("NO")) {
                    LOGGER.fatal(CLASS_HEADER + " Error, refpanel_combine = " + myTempCombine);
                    LOGGER.fatal(CLASS_HEADER + "        It should be: YES/NO");
                    System.exit(1);
                }
                // Now, if myTempCombine is "YES" and the refPanelNumber > 1, then refPanelCombine = true
                if (myTempCombine.equals("YES") && refPanelNumber > 1) {
                    refPanelCombine = true;
                } else {
                    refPanelCombine = false;
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        for (int kk = 0; kk < refPanelNumber; kk++) {
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("refpanel_type")) {
                    rpanelTypes.add(myArgument[1]);
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }
            tmpArg = argumentsArray.get(i++);
            myArgument = tmpArg.split("=");
            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                if (myArgument[0].equals("refpanel_file_dir")) {
                    rpanelDir.add(myArgument[1]);
                    checkExistence(myArgument[1]);
                    String tmpRpanelDir = myArgument[1];
                    ArrayList<String> chromoListRpanelHapFileName = new ArrayList<>();
                    for (int j = start; j <= end; j++) {
                        tmpArg = argumentsArray.get(i++);
                        myArgument = tmpArg.split("=");
                        if ((myArgument.length > 0) && (myArgument.length < 3)) {
                            String chromo = Integer.toString(j);
                            String tmpfile = "refpanel_hap_file_chr_" + chromo;
                            if (myArgument[0].equals(tmpfile)) {
                                chromoListRpanelHapFileName.add(myArgument[1]);
                                checkExistence(tmpRpanelDir + File.separator + myArgument[1]);
                            } else {
                                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                                System.exit(1);
                            }
                        } else {
                            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                            System.exit(1);
                        }
                    }
                    rpanelHapFileName.add(chromoListRpanelHapFileName);

                    // We have to know which is the imputaton tool, because:
                    // if we are going to use impute, then we need the legendfiles,
                    // if we are going to use minimac, we do not need legendfiles
                    if (imputationTool.equals("impute")) {
                        exclSVSnp = "NO";
                        LOGGER.info(CLASS_HEADER + " We are going to use 'impute' tool for imputation stage... ");
                        ArrayList<String> chromoListRpanelLegFileName = new ArrayList<>();
                        for (int j = start; j <= end; j++) {
                            tmpArg = argumentsArray.get(i++);
                            myArgument = tmpArg.split("=");
                            if ((myArgument.length > 0) && (myArgument.length < 3)) {
                                String chromo = Integer.toString(j);
                                String tmpfile = "refpanel_leg_file_chr_" + chromo;
                                if (myArgument[0].equals(tmpfile)) {
                                    chromoListRpanelLegFileName.add(myArgument[1]);
                                    checkExistence(tmpRpanelDir + File.separator + myArgument[1]);
                                } else {
                                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                                    System.exit(1);
                                }
                            } else {
                                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                                System.exit(1);
                            }
                        }
                        rpanelLegFileName.add(chromoListRpanelLegFileName);
                    } else if (imputationTool.equals("minimac")) {
                        exclSVSnp = "YES";
                        LOGGER.info(CLASS_HEADER + " We are going to use 'minimac' tool for imputation stage... ");
                    } else {
                        LOGGER.fatal(CLASS_HEADER
                                + " Sorry, Only 'impute' or 'minimac' Tools are supported right now for the imputation process... ");
                        LOGGER.fatal(CLASS_HEADER + " Future version will support other tools... ");
                        System.exit(1);
                    }
                } else {
                    LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                    System.exit(1);
                }
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
                System.exit(1);
            }
        }
        tmpArg = argumentsArray.get(i++);
        myArgument = tmpArg.split("=");
        if ((myArgument.length > 0) && (myArgument.length < 3)) {
            if (myArgument[0].equals("outputdir")) {
                outDir = myArgument[1];
                checkExistence(outDir);
            } else {
                LOGGER.fatal(CLASS_HEADER + ERROR_PARAM_ORDER + myArgument[0]);
                System.exit(1);
            }
        } else {
            LOGGER.fatal(CLASS_HEADER + ERROR_SYNTAX + gwasConfigFile + ERROR_SYNTAX_SUFFIX + myArgument[0]);
            System.exit(1);
        }

        // printInputCmd(inputFormat);

        // Finally we define which combination of GWAS analysis we are doing:
        // Options are: MIXED with GEN files
        // MIXED with BED files
        // SEPARATED with GEN files
        // SEPARATED with BED files

        // Now, depending on the wfDeepRequired, we have to activate the corresponding stages.
        activateStages(wfDeepRequired, imputationTool);

    } // End of ParseCmdLine class

    /**
     * Method to get the workflow deep the user wants to run
     * 
     * @return
     */
    public String getWfDeepRequired() {
        return this.wfDeepRequired;
    }

    /**
     * Method to get the workflow deep value the user wants to run
     * 
     * @param wfDeepRequired
     * @return
     */
    public int getWfDeepValue(String wfDeepRequired) {
        return this.wfPossibleDeeps.get(wfDeepRequired);
    }

    /**
     * Method to get start gen information
     * 
     * @return
     */
    public int getStart() {
        return this.start;
    }

    /**
     * Method to get end gen information
     * 
     * @return
     */
    public int getEnd() {
        return this.end;
    }

    /**
     * Method to get end gen information
     * 
     * @return
     */
    public int getEndNormal() {
        return this.endNormal;
    }

    /**
     * Method to get mafThreshold flag
     * 
     * @return
     */
    public Double getMafThreshold() {
        return this.mafThreshold;
    }

    /**
     * Method to get infoThreshold flag
     * 
     * @return
     */
    public Double getInfoThreshold() {
        return this.infoThreshold;
    }

    /**
     * Method to get hweThreshold flag
     * 
     * @return
     */
    public Double getHweCohortThreshold() {
        return this.hweCohortThreshold;
    }

    /**
     * Method to get hweCasesThreshold flag
     * 
     * @return
     */
    public Double getHweCasesThreshold() {
        return this.hweCasesThreshold;
    }

    /**
     * Method to get hweControlshreshold flag
     * 
     * @return
     */
    public Double getHweControlsThreshold() {
        return this.hweControlsThreshold;
    }

    /**
     * Method to get exclCgatSnp flag
     * 
     * @return
     */
    public String getExclCgatSnp() {
        return this.exclCgatSnp;
    }

    /**
     * Method to get exclSVSnp flag
     * 
     * @return
     */
    public String getExclSVSnp() {
        return this.exclSVSnp;
    }

    /**
     * Method to get imputationTool flag
     * 
     * @return
     */
    public String getImputationTool() {
        return this.imputationTool;
    }

    /**
     * Method to get testTypesNames flag
     * 
     * @param testNameIndex
     * @return
     */
    public String getTestTypeName(int testNameIndex) {
        return this.testTypesNames.get(testNameIndex);
    }

    /**
     * Method to get the number of testTypes
     * 
     * @return
     */
    public int getNumberOfTestTypeName() {
        return this.testTypesNames.size();
    }

    /**
     * Method to get the name of covariables to included in the snptest
     * 
     * @param testNameIndex
     * @return
     */
    public String getCovariables(int testNameIndex) {
        return this.covariables.get(testNameIndex);
    }

    /**
     * Method to get the name of responseVars to included in the snptest
     * 
     * @param testNameIndex
     * @return
     */
    public String getResponseVar(int testNameIndex) {
        return this.responseVars.get(testNameIndex);
    }

    /**
     * Method to get the chunk size to be used in the analysis
     * 
     * @return
     */
    public int getChunkSize() {
        return this.chunkSize;
    }

    /**
     * Method to get the list of stages
     * 
     * @return
     */
    public String getListOfStagesFile() {
        return this.listOfStagesFile;
    }

    /**
     * Method to get the input format
     * 
     * @return
     */
    public String getInputFormat() {
        return this.inputFormat;
    }

    /**
     * Method to get the value of the removeTemporalFiles variable It can be YES or NO only
     * 
     * @return
     */
    public String getRemoveTemporalFiles() {
        return this.removeTemporalFiles;
    }

    /**
     * Method to get the value of the compressFiles variable It can be YES or NO only
     * 
     * @return
     */
    public String getCompressFiles() {
        return this.compressFiles;
    }

    /**
     * Method to get bedFile information
     * 
     * @return
     */
    public String getBedFileName() {
        return this.mixedBedFileName;
    }

    /**
     * Method to get bimFile information
     * 
     * @return
     */
    public String getBimFileName() {
        return this.mixedBimFileName;
    }

    /**
     * Method to get famFile information
     * 
     * @return
     */
    public String getFamFileName() {
        return this.mixedFamFileName;
    }

    /**
     * Method to get Cohort information
     * 
     * @return
     */
    public String getCohort() {
        return this.mixedCohort;
    }

    // Cases information
    /**
     * Method to get the mixedBedDir information
     * 
     * @return
     */
    public String getBedDir() {
        return this.mixedBedDir;
    }

    // Cases information
    /**
     * Method to get the mixedChrDir information
     * 
     * @return
     */
    public String getChrDir() {
        return this.mixedChrDir;
    }

    /**
     * Method to get the xxxxxGenFileName
     * 
     * @param chromo
     * @return
     */
    public String getGenFileName(int chromo) {
        int index = chromo - getStart();
        if (index < 0) {
            LOGGER.fatal(CLASS_HEADER + " Error, the chromosome number should be > " + getStart());
            System.exit(1);
            return "none";
        }

        return this.mixedGenFileName.get(index);
    }

    /**
     * Method to get the xxxxxGenFile
     * 
     * @param chromo
     * @return
     */
    public String getGenFile(int chromo) {
        int index = chromo - getStart();

        if (index < 0) {
            LOGGER.fatal(CLASS_HEADER + " Error, the chromosome number should be > " + getStart());
            System.exit(1);
            return "none";
        }

        return this.mixedGenFile.get(index);
    }

    /**
     * Method to get sampleDir information
     * 
     * @return
     */
    public String getSampleDir() {
        return this.mixedSampleDir;
    }

    /**
     * Method to get the xxxxxSampleFileName
     * 
     * @return
     */
    public String getSampleFileName() {
        return this.mixedSampleFileName;
    }

    /**
     * Method to get the xxxxxSampleFile
     * 
     * @return
     */
    public String getSampleFile() {
        return this.mixedSampleFile;
    }

    // General information
    /**
     * Method to get gmapDir information
     * 
     * @return
     */
    public String getGmapDir() {
        return this.gmapDir;
    }

    /**
     * Method to get the gmapFileName
     * 
     * @param chromo
     * @return
     */
    public String getGmapFileName(int chromo) {
        int index = chromo - getStart();
        if (index < 0) {
            LOGGER.fatal(CLASS_HEADER + " Error, the chromosome number should be > " + getStart());
            System.exit(1);
            return "none";
        }
        return this.gmapFileName.get(index);
    }

    /**
     * Method to get refPanelNumber information
     * 
     * @return
     */
    public int getRpanelNumber() {
        return this.refPanelNumber;
    }

    /**
     * Method to get refPanelCombine variable
     * 
     * @return
     */
    public boolean getRefPanelCombine() {
        return this.refPanelCombine;
    }

    /**
     * Method to get refPanelType information
     * 
     * @param indexRpanel
     * @return
     */
    public String getRpanelType(int indexRpanel) {
        return this.rpanelTypes.get(indexRpanel);
    }

    /**
     * Method to get rpanelDir information
     * 
     * @param indexRpanel
     * @return
     */
    public String getRpanelDir(int indexRpanel) {
        return this.rpanelDir.get(indexRpanel);

    }

    /**
     * Method to get ALL the rpanelTypes in an ArrayList
     * 
     * @return
     */
    public ArrayList<String> getRpanelTypes() {
        return this.rpanelTypes;
    }

    /**
     * Method to get the rPanelHapFileName
     * 
     * @param indexRpanel
     * @param chromo
     * @return
     */
    public String getRpanelHapFileName(int indexRpanel, int chromo) {
        int index = chromo - getStart();
        if (chromo < getStart() || chromo > getEnd()) {
            LOGGER.fatal(CLASS_HEADER + " Error, the chromosome number should be " + getStart() + " <= chromo <= " + getEnd());
            System.exit(1);
            return "none";
        }

        if ((indexRpanel < 0) || (indexRpanel > rpanelTypes.size())) {
            LOGGER.fatal(CLASS_HEADER + " Error, the indexRpanel should be  0<= indexRpanel<=" + rpanelTypes.size());
            System.exit(1);
            return "none";
        }

        return this.rpanelHapFileName.get(indexRpanel).get(index);
    }

    /**
     * Method to get the rPanelLegFileName
     * 
     * @param indexRpanel
     * @param chromo
     * @return
     */
    public String getRpanelLegFileName(int indexRpanel, int chromo) {
        int index = chromo - getStart();
        if (chromo < getStart() || chromo > getEnd()) {
            LOGGER.fatal(CLASS_HEADER + " Error, the chromosome number should be " + getStart() + " <= chromo <= " + getEnd());
            System.exit(1);
            return "none";
        }

        if ((indexRpanel < 0) || (indexRpanel > rpanelTypes.size())) {
            LOGGER.fatal(CLASS_HEADER + " Error, the indexRpanel should be  0<= indexRpanel<=" + rpanelTypes.size());
            System.exit(1);
            return "none";
        }

        return this.rpanelLegFileName.get(indexRpanel).get(index);
    }

    /**
     * Method to get outDir information
     * 
     * @return
     */
    public String getOutDir() {
        return this.outDir;
    }

    /*
     * // Method for printing the help public void printHelp() { LOGGER.info("Usage: Guidance [options]");
     * LOGGER.info("All Possible options:");
     * LOGGER.info("   -gmapdir          : Directroy where the genetic map file is stored.");
     * LOGGER.info("   -rpaneldir        : Directroy where the directories of the reference panels are stored");
     * LOGGER.info("   -rpaneltypes      : The current reference panels supported are:");
     * LOGGER.info("                       1kg, dceg, hapmap, 1kg_xxx");
     * LOGGER.info("   -exclude_cgat_snps: Flag to indicate if we want to have CGAT exclusions:");
     * LOGGER.info("                     : YES/NO");
     * LOGGER.info("   -outdir           : Directory where the GWAS results will be stored.");
     * LOGGER.info("   -start            : First chromosome to process, from 1 to 23.");
     * LOGGER.info("   -end              : Last chromosome to process, from 1 to 23.");
     * LOGGER.info("   -totalprocs       : Total number of processors to use in the execution.");
     * LOGGER.info("   -chrdir           : Directory where the chromosomes for cases are stores");
     * LOGGER.info("   -sampledir        : Directory where the sample files for cases are stored.");
     * LOGGER.info("   -cohort           : The current supported cohorts for cases are:");
     * LOGGER.info("                       58C, NBS, T2D");
     * LOGGER.info("   -gmapdir          : Directroy where the genetic map files are stored.");
     * LOGGER.info("   -gmaptype         : The current genetic maps supported are:");
     * LOGGER.info("                       1kg, dceg, hapmap");
     * LOGGER.info("   -exclude_cgat_snps: Flag to indicate if we want to have CGAT exclusions:");
     * LOGGER.info("                     : YES/NO");
     * LOGGER.info("   -outdir           : Directory where the GWAS results will be stored.");
     * LOGGER.info("   -start            : First chromosome to process, from 1 to 23.");
     * LOGGER.info("   -end              : Last chromosome to process, from 1 to 23.");
     * LOGGER.info("---------------------------------"); }
     */

    /**
     * Method for printing the input command line
     * 
     * @param inputFormat
     */
    public void printInputCmd(String inputFormat) {
        LOGGER.info("[ParseCmdLine] Execution will be done with the following input parameters:");
        LOGGER.info("init_chromosome              = " + start);
        LOGGER.info("end_chromosome               = " + end);
        LOGGER.info("maf_threshold                = " + mafThreshold);
        LOGGER.info("info_threshold               = " + infoThreshold);
        LOGGER.info("hwe_cohort_threshold         = " + hweCohortThreshold);
        LOGGER.info("hwe_cases_threshold          = " + hweCasesThreshold);
        LOGGER.info("hwe_controls_threshold       = " + hweControlsThreshold);
        LOGGER.info("exclude_cgat_snps            = " + exclCgatSnp);
        int number_of_tests = testTypesNames.size();

        for (int kk = 0; kk < number_of_tests; kk++) {
            String tmp_test_type = testTypesNames.get(kk);
            String tmp_responseVar = responseVars.get(kk);
            String tmp_covariables = covariables.get(kk);
            LOGGER.info(tmp_test_type + " = " + tmp_responseVar + ":" + tmp_covariables);
        }
        LOGGER.info("imputation_tool              = " + imputationTool);
        LOGGER.info("names_of_covariables         = " + covariables);
        LOGGER.info("chunk_size_analysis          = " + chunkSize);
        LOGGER.info("file_name_for_list_of_stages = " + listOfStagesFile);

        LOGGER.info("mixed_cohort                 = " + mixedCohort);

        if (inputFormat.equals("GEN")) {
            LOGGER.info("mixed_gen_file_dir           = " + mixedBedDir);
            for (int kk = getStart(); kk <= getEnd(); kk++) {
                int index = kk - getStart();
                LOGGER.info("\tmixed_gen_file_chr_" + kk + " = " + mixedGenFileName.get(index));
            }

            LOGGER.info("mixed_sample_file_dir = " + mixedSampleDir);
            LOGGER.info("\tmixed_sample_file = " + mixedSampleFileName);
        } else if (inputFormat.equals("BED")) {
            LOGGER.info("mixed_bed_file_dir           = " + mixedBedDir);
        } else {
            LOGGER.fatal(CLASS_HEADER + " Error, input format " + inputFormat + " is not supported");
            System.exit(1);
        }

        LOGGER.info("genmap_file_dir = " + gmapDir);
        for (int kk = getStart(); kk <= getEnd(); kk++) {
            int index = kk - getStart();
            LOGGER.info("\tgenmap_file_chr_" + kk + " = " + gmapFileName.get(index));
        }

        LOGGER.info("refpanel_number = " + refPanelNumber);
        for (int kk = 0; kk < rpanelTypes.size(); kk++) {
            LOGGER.info("\trefpanel_type = " + rpanelTypes.get(kk));
            LOGGER.info("\trefpanel_file_dir = " + rpanelDir.get(kk));
            for (int jj = getStart(); jj <= getEnd(); jj++) {
                int index = jj - getStart();
                LOGGER.info("\t\trefpanel_hap_file_chr_" + jj + " = " + rpanelHapFileName.get(kk).get(index));
            }
            for (int jj = getStart(); jj <= getEnd(); jj++) {
                int index = jj - getStart();
                LOGGER.info("\t\trefpanel_leg_file_chr_" + jj + " = " + rpanelLegFileName.get(kk).get(index));
            }
        }
        LOGGER.info("\toutputdir = " + outDir);
        LOGGER.info("------------------------------------");
    }

    /**
     * Method to activate stages of the workflow depending on the wfDeepRequired string given by the user.
     * 
     * @param wfDeepRequired
     * @param imputationTool
     */
    public void activateStages(String wfDeepRequired, String imputationTool) {

        // Lets create the complete list of stages of the workflow
        // With the value 0, meaning initially they are not active.
        // They will be activated when activateStages is called, at the end of the class.

        final Integer DISABLED_MASK = 0;
        wfAllStages.put("convertFromBedToBed", DISABLED_MASK);
        wfAllStages.put("createRsIdList", DISABLED_MASK);
        wfAllStages.put("phasingBed", DISABLED_MASK);
        wfAllStages.put("phasing", DISABLED_MASK);
        wfAllStages.put("createListOfExcludedSnps", DISABLED_MASK);
        wfAllStages.put("filterHaplotypes", DISABLED_MASK);
        wfAllStages.put("imputeWithImpute", DISABLED_MASK);
        wfAllStages.put("imputeWithMinimac", DISABLED_MASK);
        wfAllStages.put("filterByInfo", DISABLED_MASK);
        wfAllStages.put("qctoolS", DISABLED_MASK);
        wfAllStages.put("snptest", DISABLED_MASK);
        wfAllStages.put("collectSummary", DISABLED_MASK);
        wfAllStages.put("mergeTwoChunks", DISABLED_MASK);
        wfAllStages.put("filterByAll", DISABLED_MASK);
        wfAllStages.put("jointCondensedFiles", DISABLED_MASK);
        wfAllStages.put("jointFilteredByAllFiles", DISABLED_MASK);
        wfAllStages.put("generateTopHits", DISABLED_MASK);
        wfAllStages.put("generateQQManhattanPlots", DISABLED_MASK);
        wfAllStages.put("combinePanelsComplex", DISABLED_MASK);
        wfAllStages.put("combineCondensedFiles", DISABLED_MASK);
        wfAllStages.put("initPhenoMatrix", DISABLED_MASK);
        wfAllStages.put("addToPhenoMatrix", DISABLED_MASK);
        wfAllStages.put("filloutPhenoMatrix", DISABLED_MASK);
        wfAllStages.put("finalizePhenoMatrix", DISABLED_MASK);
        wfAllStages.put("taskx", DISABLED_MASK);
        wfAllStages.put("tasky", DISABLED_MASK);
        wfAllStages.put("taskz", DISABLED_MASK);

        // First of all, we activate the correct PossibleDeeps depending on the kind of imputation tool that is used:
        // Important: The order of the bits in the value are related to the order of the stages in the wfAllStages
        // Hastable below. If you are going to modify the list of stages you should fix the new binary value in
        // wfPossibleDeeps.
        if (imputationTool.equals("impute")) {
            wfPossibleDeeps.put("until_convertFromBedToBed", 0x6000000);
            wfPossibleDeeps.put("until_phasing", 0x7800000);
            wfPossibleDeeps.put("until_imputation", 0x7900000);
            wfPossibleDeeps.put("until_association", 0x7970000);
            wfPossibleDeeps.put("until_filterByAll", 0x797E000);
            wfPossibleDeeps.put("until_summary", 0x797FF80);
            wfPossibleDeeps.put("whole_workflow", 0x797FFF8);
            wfPossibleDeeps.put("from_phasing", 0x017FFF8);
            wfPossibleDeeps.put("from_phasing_to_summary", 0x017FF80);
            wfPossibleDeeps.put("from_phasing_to_filterByAll", 0x017E000);
            wfPossibleDeeps.put("from_phasing_to_association", 0x0170000);
            wfPossibleDeeps.put("from_phasing_to_imputation", 0x0100000);
            wfPossibleDeeps.put("from_imputation", 0x007FFF8);
            wfPossibleDeeps.put("from_imputation_to_summary", 0x007FF80);
            wfPossibleDeeps.put("from_imputation_to_filterByAll", 0x007E000);
            wfPossibleDeeps.put("from_imputation_to_association", 0x0070000);
            wfPossibleDeeps.put("from_imputation_to_filterByInfo", 0x0040000);
            wfPossibleDeeps.put("from_filterByInfo_to_qctoolS", 0x0020000);
            wfPossibleDeeps.put("from_qctoolS_to_association", 0x0010000);
            wfPossibleDeeps.put("from_association", 0x000FFF8);
            wfPossibleDeeps.put("from_association_to_filterByAll", 0x000E000);
            wfPossibleDeeps.put("from_association_to_summary", 0x000FF80);
            wfPossibleDeeps.put("from_filterByAll", 0x0001FF8);
            wfPossibleDeeps.put("from_filterByAll_to_summary", 0x0001F80);
            wfPossibleDeeps.put("from_summary", 0x0000078);
        } else if (imputationTool.equals("minimac")) {
            wfPossibleDeeps.put("until_convertFromBedToBed", 0x6000000);
            wfPossibleDeeps.put("until_phasing", 0x7E00000);
            wfPossibleDeeps.put("until_imputation", 0x7E80000);
            wfPossibleDeeps.put("until_association", 0x7EF0000);
            wfPossibleDeeps.put("until_filterByAll", 0x7EFE000);
            wfPossibleDeeps.put("until_summary", 0x7EFFF80);
            wfPossibleDeeps.put("whole_workflow", 0x7EFFFF8);
            wfPossibleDeeps.put("from_phasing", 0x00FFFF8);
            wfPossibleDeeps.put("from_phasing_to_summary", 0x00FFF80);
            wfPossibleDeeps.put("from_phasing_to_filterByAll", 0x00FE000);
            wfPossibleDeeps.put("from_phasing_to_association", 0x00F0000);
            wfPossibleDeeps.put("from_phasing_to_imputation", 0x0080000);
            wfPossibleDeeps.put("from_imputation", 0x007FFF8);
            wfPossibleDeeps.put("from_imputation_to_summary", 0x007FF80);
            wfPossibleDeeps.put("from_imputation_to_filterByAll", 0x007E000);
            wfPossibleDeeps.put("from_imputation_to_association", 0x0070000);
            wfPossibleDeeps.put("from_imputation_to_filterByInfo", 0x0040000);
            wfPossibleDeeps.put("from_filterByInfo_to_qctoolS", 0x0020000);
            wfPossibleDeeps.put("from_qctoolS_to_association", 0x0010000);
            wfPossibleDeeps.put("from_association", 0x000FFF8);
            wfPossibleDeeps.put("from_association_to_filterByAll", 0x000E000);
            wfPossibleDeeps.put("from_association_to_summary", 0x000FF80);
            wfPossibleDeeps.put("from_filterByAll", 0x0001FF8);
            wfPossibleDeeps.put("from_filterByAll_to_summary", 0x0001F80);
            wfPossibleDeeps.put("from_summary", 0x0000078);
        } else {
            LOGGER.fatal(CLASS_HEADER + " Error, the imputation tool: " + imputationTool + " is not supported in this version...");
            System.exit(1);
        }

        final Integer MASK1 = 0x00001;
        int stageNumber = 0;
        // Shift 1 and Mask1
        int tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("taskz", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("tasky", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("taskx", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("finalizePhenoMatrix", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("filloutPhenoMatrix", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("addToPhenoMatrix", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("initPhenoMatrix", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("combineCondensedFiles", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("combinePanelsComplex", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("generateQQManhattanPlots", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("generateTopHits", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("jointFilteredByAllFiles", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("jointCondensedFiles", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("filterByAll", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("mergeTwoChunks", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("collectSummary", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("snptest", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("qctoolS", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("filterByInfo", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("imputeWithMinimac", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("imputeWithImpute", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("filterHaplotypes", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("createListOfExcludedSnps", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("phasing", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("phasingBed", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("createRsIdList", tmpVar);
        stageNumber++;

        tmpVar = (wfPossibleDeeps.get(wfDeepRequired) >> stageNumber) & MASK1;
        wfAllStages.put("convertFromBedToBed", tmpVar);
    }

    /**
     * Method to get the status of a stage. That is if it is active (1) or unactive (0)
     * 
     * @param myStage
     * @return
     */
    public int getStageStatus(String myStage) {
        return this.wfAllStages.get(myStage);
    }

    /**
     * Method to check the existence of a file or a directory defined in the configuration file
     * 
     * @param dirOrFileName
     */
    public void checkExistence(String dirOrFileName) {
        File theDir = new File(dirOrFileName);
        if (!theDir.exists()) {
            LOGGER.fatal(CLASS_HEADER + " Error, " + dirOrFileName + " does not exist!");
            LOGGER.fatal(CLASS_HEADER + "        Please verify the existence of all your input data set.");
            System.exit(1);
        }
    }

}
