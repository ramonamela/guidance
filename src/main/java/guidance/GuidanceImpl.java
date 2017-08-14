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

package guidance;

import java.util.zip.GZIPOutputStream;

import guidance.exceptions.GuidanceTaskException;

import java.io.File;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;

import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Hashtable;

import java.util.Set;

import java.nio.channels.FileChannel;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import static java.nio.file.StandardCopyOption.*;
import java.io.RandomAccessFile;
import java.util.zip.GZIPInputStream;

import java.util.Scanner;


public class GuidanceImpl {

    // Debug
    private static final boolean DEBUG = true;

    private static final String HEADER_MIXED = "chr	position	rs_id_all	info_all	certainty_all	alleleA	alleleB	index	average_maximum_posterior_call	info	cohort_1_AA	cohort_1_AB	cohort_1_BB	cohort_1_NULL	all_AA	all_AB	all_BB	all_NULL	all_total	cases_AA	cases_AB	cases_BB	cases_NULL	cases_total	controls_AA	controls_AB	controls_BB	controls_NULL	controls_total	all_maf	cases_maf	controls_maf	missing_data_proportion	cohort_1_hwe	cases_hwe	controls_hwe	het_OR	het_OR_lower	het_OR_upper	hom_OR	hom_OR_lower	hom_OR_upper	all_OR	all_OR_lower	all_OR_upper	frequentist_add_pvalue	frequentist_add_info	frequentist_add_beta_1	frequentist_add_se_1	frequentist_dom_pvalue	frequentist_dom_info	frequentist_dom_beta_1	frequentist_dom_se_1	frequentist_rec_pvalue	frequentist_rec_info	frequentist_rec_beta_1	frequentist_rec_se_1	frequentist_gen_pvalue	frequentist_gen_info	frequentist_gen_beta_1	frequentist_gen_se_1	frequentist_gen_beta_2	frequentist_gen_se_2	frequentist_het_pvalue	frequentist_het_info	frequentist_het_beta_1	frequentist_het_se_1	comment";
    private static final String HEADER_MIXED_X = "chr	position	rs_id_all	info_all	certainty_all	alleleA	alleleB	all_A	all_B	all_AA	all_AB	all_BB	all_NULL	all_total	all_maf	all_info	all_impute_info	cases_A	cases_B	cases_AA	cases_AB	cases_BB	cases_NULL	cases_total	cases_maf	cases_info	cases_impute_info	controls_A	controls_B	controls_AA	controls_AB	controls_BB	controls_NULL	controls_total	controls_maf	controls_info	controls_impute_info	sex=1_A	sex=1_B	sex=1_AA	sex=1_AB	sex=1_BB	sex=1_NULL	sex=1_total	sex=1_maf	sex=1_info	sex=1_impute_info	sex=2_A	sex=2_B	sex=2_AA	sex=2_AB	sex=2_BB	sex=2_NULL	sex=2_total	sex=2_maf	sex=2_info	sex=2_impute_info	frequentist_add_null_ll	frequentist_add_alternative_ll	frequentist_add_beta_1:genotype/sex=1	frequentist_add_beta_2:genotype/sex=2	frequentist_add_se_1:genotype/sex=1	frequentist_add_se_2:genotype/sex=2	frequentist_add_degrees_of_freedom	frequentist_add_pvalue	comment";


    /**
     * Method to perform the conversion from Bed to Ped Format file
     * 
     * @param bedFile
     * @param bimFile
     * @param famFile
     * @param newBedFile
     * @param newBimFile
     * @param newFamFile
     * @param logFile
     * @param chromo
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void convertFromBedToBed(String bedFile, String bimFile, String famFile, String newBedFile, String newBimFile,
            String newFamFile, String logFile, String chromo, String cmdToStore) throws GuidanceTaskException {

        String plinkBinary = System.getenv("PLINKBINARY");
        if (plinkBinary == null) {
            throw new GuidanceTaskException("[convertFromBedToBed] Error, PLINKBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running convertFromBedToBed with parameters:");
            System.out.println("[DEBUG] \t- Input bedFile   : " + bedFile);
            System.out.println("[DEBUG] \t- Input bimFile   : " + bimFile);
            System.out.println("[DEBUG] \t- Input famFile   : " + famFile);
            System.out.println("[DEBUG] \t- Output newBedFile  : " + newBedFile);
            System.out.println("[DEBUG] \t- Output newBimFile  : " + newBimFile);
            System.out.println("[DEBUG] \t- Output newFamFile  : " + newFamFile);
            System.out.println("[DEBUG] \t- Output logFile  : " + logFile);
            System.out.println("[DEBUG] \t- Chromosome      : " + chromo);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Workers:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            System.out.println("--------------------------------------");
        }

        long startTime = System.currentTimeMillis();
        String cmd = null;

        String tmpBedFile = newBedFile + ".bed";
        String tmpBimFile = newBedFile + ".bim";
        String tmpFamFile = newBedFile + ".fam";
        String tmpLogFile = newBedFile + ".log";

        cmd = plinkBinary + " --noweb --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile + " --chr " + chromo
                + " --recode --out " + newBedFile + " --make-bed";

        if (DEBUG) {
            System.out.println("[convertFromBedToBed] Command: " + cmd);
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process convertFromBedToBedProc = null;
        try {
            pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(convertFromBedToBedProc.getInputStream(), newBedFile + ".stdout", convertFromBedToBedProc.getErrorStream(),
                    newBedFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = convertFromBedToBedProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException("[convertFromBedToBed] Error executing convertFromBedToBed job, exit value is: " + exitValue);
        }

        // File (or directory) with old name
        File file = new File(newBedFile + ".bed");
        File file2 = new File(newBedFile);
        // if(file2.exists()) throw new java.io.IOException("file exists");
        // Rename file (or directory)
        boolean success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToBed] Error, the file " + newBedFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        file = new File(newBedFile + ".bim");
        file2 = new File(newBimFile);
        // Rename file (or directory)
        success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToBed] Error, the file " + newBimFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        file = new File(newBedFile + ".fam");
        file2 = new File(newFamFile);
        // Rename file (or directory)
        success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToBed] Error, the file " + newFamFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        file = new File(newBedFile + ".log");
        file2 = new File(logFile);
        // Rename file (or directory)
        success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToBed] Error, the file " + logFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        // If there is not output in the convertFromBedToBed process. Then we have to create some empty outputs.
        File fa = new File(newBedFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToBed] The file " + newBedFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        fa = new File(newBimFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToBed] The file " + newBimFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        fa = new File(newFamFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToBed] The file " + newFamFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        fa = new File(logFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToBed] The file " + logFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] convertFromBedToBed startTime: " + startTime);
            System.out.println("\n[DEBUG] convertFromBedToBed endTime: " + stopTime);
            System.out.println("\n[DEBUG] convertFromBedToBed elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of convertFromBedToBed.");
        }
    }

    /**
     * Method to perform the conversion from Bed to Ped Format file
     * 
     * @param bedPrefix
     * @param bedFile
     * @param bimFile
     * @param famFile
     * @param pedFile
     * @param mapFile
     * @param logFile
     * @param chromo
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void convertFromBedToPed(String bedPrefix, String bedFile, String bimFile, String famFile, String pedFile, String mapFile,
            String logFile, String chromo, String cmdToStore) throws GuidanceTaskException {

        String plinkBinary = System.getenv("PLINKBINARY");
        if (plinkBinary == null) {
            throw new GuidanceTaskException("[convertFromBedToPed] Error, PLINKBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running convertFromBedToPed with parameters:");
            System.out.println("[DEBUG] \t- Input bedPrefix : " + bedPrefix);
            System.out.println("[DEBUG] \t- Input bedFile   : " + bedFile);
            System.out.println("[DEBUG] \t- Input bimFile   : " + bimFile);
            System.out.println("[DEBUG] \t- Input famFile   : " + famFile);
            System.out.println("[DEBUG] \t- Output pedFile  : " + pedFile);
            System.out.println("[DEBUG] \t- Output mapFile  : " + mapFile);
            System.out.println("[DEBUG] \t- Output logFile  : " + logFile);
            System.out.println("[DEBUG] \t- Chromosome      : " + chromo);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- command : " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Workers:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            System.out.println("--------------------------------------");
        }

        long startTime = System.currentTimeMillis();
        String cmd = null;
        cmd = plinkBinary + " --noweb --bfile " + bedPrefix + " --chr " + chromo + " --recode --out " + pedFile;

        if (DEBUG) {
            System.out.println("[convertFromBedToPed] Command: " + cmd);
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process convertFromBedToPedProc;
        try {
            convertFromBedToPedProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(convertFromBedToPedProc.getInputStream(), pedFile + ".stdout", convertFromBedToPedProc.getErrorStream(),
                    pedFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = convertFromBedToPedProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException("[convertFromBedToPed] Error executing convertFromBedToPed job, exit value is: " + exitValue);
        }

        // File (or directory) with old name
        File file = new File(pedFile + ".ped");
        File file2 = new File(pedFile);
        // if(file2.exists()) throw new java.io.IOException("file exists");

        // Rename file (or directory)
        boolean success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToPed] Error, the file " + pedFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        file = new File(pedFile + ".map");
        file2 = new File(mapFile);

        // Rename file (or directory)
        success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToPed] Error, the file " + mapFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        file = new File(pedFile + ".log");
        file2 = new File(logFile);

        // Rename file (or directory)
        success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromBedToPed] Error, the file " + logFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        // If there is not output in the convertFromBedToPed process. Then we have to create some empty outputs.
        File fa = new File(mapFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToPed] The file " + mapFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        fa = new File(pedFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToPed] The file " + pedFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        fa = new File(logFile);
        if (!fa.exists()) {
            System.out.println("[convertFromBedToPed] The file " + logFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] convertFromBedToPed startTime: " + startTime);
            System.out.println("\n[DEBUG] convertFromBedToPed endTime: " + stopTime);
            System.out.println("\n[DEBUG] convertFromBedToPed elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of convertFromBedToPed.");
        }

    }

    /**
     * Method to perform the conversion from Ped to Gen Format file
     * 
     * @param pedFile
     * @param mapFile
     * @param genFile
     * @param sampleFile
     * @param logFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void convertFromPedToGen(String pedFile, String mapFile, String genFile, String sampleFile, String logFile,
            String cmdToStore) throws GuidanceTaskException {

        String gtoolBinary = System.getenv("GTOOLBINARY");
        if (gtoolBinary == null) {
            throw new GuidanceTaskException("[convertFromPedToGen] Error, GTOOLBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running convertFromPedToGen with parameters:");
            System.out.println("[DEBUG] \t- Input pedFile      : " + pedFile);
            System.out.println("[DEBUG] \t- Input mapFile      : " + mapFile);
            System.out.println("[DEBUG] \t- Output genFile     : " + genFile);
            System.out.println("[DEBUG] \t- Output sampleFile  : " + sampleFile);
            System.out.println("[DEBUG] \t- Output logFile     : " + logFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            System.out.println("--------------------------------------");
        }

        long startTime = System.currentTimeMillis();

        String cmd = null;
        cmd = gtoolBinary + " -P --ped " + pedFile + " --map " + mapFile + " --og " + genFile + " --os " + sampleFile
                + " --binary_phenotype --order --log " + logFile;

        if (DEBUG) {
            System.out.println("[convertFromPedToGen] Command: " + cmd);
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process convertFromPedToGenProc;
        try {
            convertFromPedToGenProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(convertFromPedToGenProc.getInputStream(), genFile + ".stdout", convertFromPedToGenProc.getErrorStream(),
                    genFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = convertFromPedToGenProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException("[convertFromPedToGen] Error executing convertFromPedToGed job, exit value is: " + exitValue);
        }

        // If there is not output in the convertFromPedToGen process. Then we have to create some empty outputs.
        File fa = new File(genFile);
        if (!fa.exists()) {
            System.out.println("[convertFromPedToGen] The file " + genFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        fa = new File(sampleFile);
        if (!fa.exists()) {
            System.out.println("[convertFromPedToGen] The file " + sampleFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        // Now we have to change the missing/unaffection/affection coding in the sample file from -9/1/2 to NA/0/1
        // Then, we go to the sample file and change the last column.

        // File changedSampleFile = new File(sampleFile + ".changed");
        // File changedSampleFile = new File("sample_file.changed");

        // changedSampleFile.createNewFile();

        File changedSampleFile = new File(sampleFile + ".changed");
        try {
            changedSampleFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // We read each line of th genFile and look for the rsID into the newPosMap
        try (FileReader fr = new FileReader(sampleFile);
                BufferedReader br = new BufferedReader(fr);
                BufferedWriter writerPos = new BufferedWriter(new FileWriter(changedSampleFile))) {

            String line = null;
            // Read the first to lines which are the headers
            line = br.readLine();
            writerPos.write(line);
            writerPos.newLine();

            line = br.readLine();
            writerPos.write(line);
            writerPos.newLine();

            while ((line = br.readLine()) != null) {
                int length = line.length();
                // System.out.println("[DEBUG]: Original line: " + line);
                String subline = line.substring(length - 1, length);
                // System.out.println("[DEBUG]: El substring es: |" + subline + "|");
                String myNewLine = null;
                if (subline.equals("0")) {
                    myNewLine = line.substring(0, length - 2) + " NA";
                    // System.out.println("[DEBUG]: NA New line : " + myNewLine);
                } else if (subline.equals("1")) {
                    myNewLine = line.substring(0, length - 2) + " 0";
                    // System.out.println("[DEBUG]: 0 New line : " + myNewLine);
                } else if (subline.equals("2")) {
                    myNewLine = line.substring(0, length - 2) + " 1";
                    // System.out.println("[DEBUG]: 1 New line : " + myNewLine);
                } else {
                    System.out.println(
                            "[convertFromPedToGen]: Error changing the sample file. Invalid affection coding in line " + myNewLine);
                    // throw new Exception("Error changing the sample file. Invalid affection coding in line " +
                    // counter);
                }
                writerPos.write(myNewLine);
                writerPos.newLine();
            }

            writerPos.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        File file = new File(sampleFile + ".changed");
        File file2 = new File(sampleFile);
        // if(file2.exists()) throw new java.io.IOException("file exists");

        // Rename file (or directory)
        boolean success = file.renameTo(file2);
        if (!success) {
            throw new GuidanceTaskException("[convertFromPedToGen] Error, the file " + sampleFile + " was not succesfully renamed");
            // File was not successfully renamed
        }

        fa = new File(logFile);
        if (!fa.exists()) {
            System.out.println("[convertFromPedToGen] The file " + logFile + " does not exist, then we create an empty file");
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] convertFromPedToGen startTime: " + startTime);
            System.out.println("\n[DEBUG] convertFromPedToGen endTime: " + stopTime);
            System.out.println("\n[DEBUG] convertFromPedToGen elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of convertFromPedToGen.");
        }

    }

    /**
     * Method to create a file that contains the SNPs positions of the input genFile
     * 
     * @param genOrBimFile
     * @param exclCgatFlag
     * @param pairsFile
     * @param inputFormat
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void createRsIdList(String genOrBimFile, String exclCgatFlag, String pairsFile, String inputFormat, String cmdToStore)
            throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running createRsIdList with parameters:");
            System.out.println("[DEBUG] \t- Input genOrBimFile : " + genOrBimFile);
            System.out.println("[DEBUG] \t- Input exclCgatFlag : " + exclCgatFlag);
            System.out.println("[DEBUG] \t- Output pairsFile   : " + pairsFile);
            System.out.println("[DEBUG] \t- InputFormat        : " + inputFormat);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // long freeMemory = Runtime.getRuntime().freeMemory()/1048576;
            // long totalMemory = Runtime.getRuntime().totalMemory()/1048576;
            // long maxMemory = Runtime.getRuntime().maxMemory()/1048576;

            // System.out.println("JVM freeMemory: " + freeMemory);
            // System.out.println("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
            // System.out.println("JVM maxMemory also equals to maximum heap size of JVM : " + maxMemory);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Worker:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            System.out.println("--------------------------------------");
        }

        // ArrayList objects = new ArrayList();

        // for (int ii = 0; ii < 10; ii++) {
        // objects.add(("" + 10 * 2710));
        // }

        // freeMemory = Runtime.getRuntime().freeMemory()/1048576;
        // totalMemory = Runtime.getRuntime().totalMemory()/1048576;
        // maxMemory = Runtime.getRuntime().maxMemory()/1048576;
        // System.out.println("Used Memory in JVM: " + (maxMemory - freeMemory));
        // System.out.println("freeMemory in JVM: " + freeMemory);
        // System.out.println("totalMemory in JVM shows current size of java heap : " + totalMemory);
        // System.out.println("maxMemory in JVM: " + maxMemory);

        long startTime = System.currentTimeMillis();

        // Check if the file is gzip or not.
        // This is done by checking the magic number of gzip files, which is 0x1f8b (the first two bytes)
        File tmpInputFile = new File(genOrBimFile);
        long n;
        try (RandomAccessFile raf = new RandomAccessFile(tmpInputFile, "r")) {
            n = raf.readInt();
            n = n & 0xFFFF0000;
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Create output file
        File outPairsFile = new File(pairsFile);
        try {
            outPairsFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        boolean thisIsGz = (n == 0x1f8b0000);

        if (thisIsGz) {
            System.out.println(
                    "[CreateRsIdList] It seems the file " + genOrBimFile + " is a gzip file. Magic Number is " + String.format("%x", n));

            try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(genOrBimFile));
                    Reader decoder = new InputStreamReader(inputGz);
                    BufferedReader br = new BufferedReader(decoder)) {

                writeOutPairsFile(br, outPairsFile, inputFormat, exclCgatFlag, pairsFile);
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        } else {
            try (FileReader fr = new FileReader(genOrBimFile); BufferedReader br = new BufferedReader(fr)) {
                writeOutPairsFile(br, outPairsFile, inputFormat, exclCgatFlag, pairsFile);
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] createRsIdList StartTime: " + startTime);
            System.out.println("\n[DEBUG] createRsIdList endTime: " + stopTime);
            System.out.println("\n[DEBUG] createRsIdList elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of createRsIdList.");
        }

    }

    private static void writeOutPairsFile(BufferedReader br, File outPairsFile, String inputFormat, String exclCgatFlag, String pairsFile)
            throws GuidanceTaskException {

        try (BufferedWriter writerPairs = new BufferedWriter(new FileWriter(outPairsFile))) {
            if (inputFormat.equals("BED")) {
                String line = null;
                while ((line = br.readLine()) != null) {
                    String[] splittedLine = line.split("\t");// delimiter I assume single space.
                    String allele = splittedLine[4] + splittedLine[5]; // store Allele (AC,GC,AG, GT,..., etc.)

                    // Store rsID of the SNP which its allele is AT or TA or GC or CG into the .pairs file
                    if (exclCgatFlag.equals("YES")) {
                        // Then we have to see if allele is AT TA GC CG to put the rsID into the .pairs file.
                        if (allele.equals("AT") || allele.equals("TA") || allele.equals("GC") || allele.equals("CG")) {
                            writerPairs.write(splittedLine[1]);
                            writerPairs.newLine();
                        }
                    }
                }
            } else if (inputFormat.equals("GEN")) {
                String line = null;
                while ((line = br.readLine()) != null) {
                    String[] splittedLine = line.split(" ");// delimiter I assume single space.
                    String allele = splittedLine[3] + splittedLine[4]; // store Allele (AC,GC,AG, GT,..., etc.)

                    // Store rsID of the SNP which its allele is AT or TA or GC or CG into the .pairs file
                    if (exclCgatFlag.equals("YES")) {
                        // Then we have to see if allele is AT TA GC CG to put the rsID into the .pairs file.
                        if (allele.equals("AT") || allele.equals("TA") || allele.equals("GC") || allele.equals("CG")) {
                            writerPairs.write(splittedLine[1]);
                            writerPairs.newLine();
                        }
                    }
                }
            } else {
                throw new GuidanceTaskException("[createRsIdList] Error, It was not possible to generate " + pairsFile + ". The "
                        + inputFormat + " is not valid!!");
            }

            writerPairs.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
    }

    /**
     * Method to execute gtool with option -S
     * 
     * @param newGenFile
     * @param modSampleFile
     * @param gtoolGenFile
     * @param gtoolSampleFile
     * @param sampleExclFile
     * @param snpWtccFile
     * @param gtoolLogFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void gtoolS(String newGenFile, String modSampleFile, String gtoolGenFile, String gtoolSampleFile, String sampleExclFile,
            String snpWtccFile, String gtoolLogFile, String cmdToStore) throws GuidanceTaskException {

        String gtoolBinary = System.getenv("GTOOLBINARY");
        if (gtoolBinary == null) {
            throw new GuidanceTaskException("[gtoolS] Error, GTOOLBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running gtoolS with parameters:");
            System.out.println("[DEBUG] \t- gtoolBinary            : " + gtoolBinary);
            System.out.println("[DEBUG] \t- Input newGenFile       : " + newGenFile);
            System.out.println("[DEBUG] \t- Input modSampleFile    : " + modSampleFile);
            System.out.println("[DEBUG] \t- Output gtoolGenFile    : " + gtoolGenFile);
            System.out.println("[DEBUG] \t- Output gtoolSampleFile : " + gtoolSampleFile);
            System.out.println("[DEBUG] \t- Input sampleExclFile   : " + sampleExclFile);
            System.out.println("[DEBUG] \t- Input snpWtccFile      : " + snpWtccFile);
            System.out.println("[DEBUG] \t- Output gtoolLogFile    : " + gtoolLogFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND            : " + cmdToStore);
            System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        String cmd = null;

        cmd = gtoolBinary + " -S --g " + newGenFile + " --s " + modSampleFile + " --og " + gtoolGenFile + " --os " + gtoolSampleFile
                + " --sample_excl " + sampleExclFile + " --exclusion " + snpWtccFile + " --log " + gtoolLogFile;

        if (DEBUG) {
            System.out.println("\n[DEBUG] Cmd -> " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process gtoolSProc;
        try {
            gtoolSProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(gtoolSProc.getInputStream(), gtoolGenFile + ".stdout", gtoolSProc.getErrorStream(),
                    gtoolGenFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = gtoolSProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException("[gtoolS] Error executing gtoolSProc job, exit value is: " + exitValue);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] gtoolS startTime: " + startTime);
            System.out.println("\n[DEBUG] gtoolS endTime: " + stopTime);
            System.out.println("\n[DEBUG] gtoolS elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of gtoolS.");
        }

    }

    /**
     * Method to execute qctool
     * 
     * @param genFile
     * @param sampleFile
     * @param qctoolGenFile
     * @param qctoolSampleFile
     * @param qctoolLogFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void qctool(String genFile, String sampleFile, String qctoolGenFile, String qctoolSampleFile, String qctoolLogFile,
            String cmdToStore) throws GuidanceTaskException {

        String qctoolBinary = System.getenv("QCTOOLBINARY");
        if (qctoolBinary == null) {
            throw new GuidanceTaskException("[qctool] Error, QCTOOLBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running qctool with parameters:");
            System.out.println("[DEBUG] \t- Input genFile          : " + genFile);
            System.out.println("[DEBUG] \t- Input sampleFile       : " + sampleFile);
            System.out.println("[DEBUG] \t- Output qctoolGenFile   : " + qctoolGenFile);
            System.out.println("[DEBUG] \t- Output qctoolSampleFile: " + qctoolSampleFile);
            System.out.println("[DEBUG] \t- Output qctoolLogFile   : " + qctoolLogFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND            : " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        String cmd = null;
        cmd = qctoolBinary + " -g " + genFile + " -s " + sampleFile + " -og " + qctoolGenFile + " -os " + qctoolSampleFile
                + " -omit-chromosome -sort -log " + qctoolLogFile;

        if (DEBUG) {
            System.out.println("\n[DEBUG] Cmd -> " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process qctoolProc;
        try {
            qctoolProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(qctoolProc.getInputStream(), qctoolGenFile + ".stdout", qctoolProc.getErrorStream(),
                    qctoolGenFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = qctoolProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException("[qctool] Error executing qctoolProc job, exit value is: " + exitValue);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] qctool startTime: " + startTime);
            System.out.println("\n[DEBUG] qctool endTime: " + stopTime);
            System.out.println("\n[DEBUG] qctool elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of qctool.");
        }
    }

    /**
     * Method to execute qctoolS
     * 
     * @param imputeFile
     * @param inclusionRsIdFile
     * @param mafThresholdS
     * @param filteredFile
     * @param filteredLogFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void qctoolS(String imputeFile, String inclusionRsIdFile, String mafThresholdS, String filteredFileGz,
            String filteredLogFile, String cmdToStore) throws GuidanceTaskException {

        String qctoolBinary = System.getenv("QCTOOLBINARY");
        if (qctoolBinary == null) {
            throw new GuidanceTaskException("[qctoolS] Error, QCTOOLBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\nRunning qctoolS for generation a subset of rsids with parameters:");
            System.out.println("\t- qctoolBinary               : " + qctoolBinary);
            System.out.println("\t- Input imputeFile           : " + imputeFile);
            System.out.println("\t- Input inclusionRsIdFile    : " + inclusionRsIdFile);
            System.out.println("\t- Input mafThreshold         : " + mafThresholdS);
            System.out.println("\t- Output filteredFile        : " + filteredFileGz);
            System.out.println("\t- Output filteredLogFile     : " + filteredLogFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND            : " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // We have to make sure whether we are using renamed files of the original gz files.
        // We detect this situation by scanning the last three characters:
        String extension = imputeFile.substring(Math.max(0, imputeFile.length() - 3));
        // System.out.println("DEBUG \t The file extension is: " + extension + " and the file is " + imputeFile);

        String imputeFileGz = null;
        if (extension.equals(".gz")) {
            imputeFileGz = imputeFile;
        } else {
            // If imputeFile exists, then imputeFileGz exists also.
            // We reused the imputFileGz
            imputeFileGz = imputeFile + ".gz";
            // String imputeFileGz = imputeFile + ".gz";
        }

        // Now we create filteredFileGz
        String filteredFile = filteredFileGz.substring(0, filteredFileGz.length() - 3);

        String cmd = null;
        cmd = qctoolBinary + " -g " + imputeFileGz + " -og " + filteredFile + " -incl-rsids " + inclusionRsIdFile
                + " -omit-chromosome -force -log " + filteredLogFile + " -maf " + mafThresholdS + " 1";

        if (DEBUG) {
            System.out.println("\n[DEBUG] Command: " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process qctoolSProc;
        try {
            qctoolSProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(qctoolSProc.getInputStream(), filteredFile + ".stdout", qctoolSProc.getErrorStream(),
                    filteredFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = qctoolSProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException("[qctoolS] Error executing qctoolSProc job, exit value is: " + exitValue);
        }

        // This tool always creates a filteredFile output (empty or not)
        // Now we need to GZ it and move to the output parameter
        gzipFile(filteredFile, filteredFileGz);

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] qctoolS startTime: " + startTime);
            System.out.println("\n[DEBUG] qctoolS endTime: " + stopTime);
            System.out.println("\n[DEBUG] qctoolS elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of qctoolS.");
        }
    }

    /**
     * Method to execute createListOfExcludedSnps
     * 
     * @param shapeitHapsFile
     * @param excludedSnpsFile
     * @param exclCgatFlag
     * @param exclSVFlag
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void createListOfExcludedSnps(String shapeitHapsFile, String excludedSnpsFile, String exclCgatFlag, String exclSVFlag,
            String cmdToStore) throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running createListOfExcludedSnps method:");
            System.out.println("[DEBUG] \t- Input shapeitHapsFile   : " + shapeitHapsFile);
            System.out.println("[DEBUG] \t- Output excludedSnpsFile : " + excludedSnpsFile);
            System.out.println("[DEBUG] \t- Input exclCgatFlag      : " + exclCgatFlag);
            System.out.println("[DEBUG] \t- Input exclSVFlag        : " + exclSVFlag);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // We have to check the haplotypes file and extract tha allele information. For the format file used by shapeit,
        // the
        // indices for alleles are:
        int rsIdIndex = 1;
        int posIndex = 2;
        int a1Index = 3;
        int a2Index = 4;
        // And the column separator is:
        String separator = " ";

        // We have to make sure whether we are using renamed files of the original gz files.
        // We detect this situation by scanning the last three characters of the file name:
        String extension = shapeitHapsFile.substring(Math.max(0, shapeitHapsFile.length() - 3));
        // System.out.println("DEBUG \t The file extension is: " + extension + " and the file is " + shapeitHapsFile);
        String shapeitHapsFileGz = null;

        if (extension.equals(".gz")) {
            shapeitHapsFileGz = shapeitHapsFile;
        } else {
            // If the renamed shapeitHapsFile exists, then an extended .gz version exists also.
            shapeitHapsFileGz = shapeitHapsFile + ".gz";
            // String imputeFileGz = imputeFile + ".gz";
        }

        // Array of string to store positions of SNPs to exclude
        ArrayList<String> excludeList = new ArrayList<>();

        // Then, we read the gz File line by line
        try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(shapeitHapsFileGz));
                Reader decoder = new InputStreamReader(inputGz);
                BufferedReader br = new BufferedReader(decoder)) {

            String line = "";
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split(separator);// delimiter defined previously.
                String allele1 = splitted[a1Index];
                String allele2 = splitted[a2Index];
                String positionS = splitted[posIndex];
                boolean excludeSNP = false;

                if (exclSVFlag.equals("YES")) {
                    if (!splitted[a1Index].equals("A") && !splitted[a1Index].equals("C") && !splitted[a1Index].equals("G")
                            && !splitted[a1Index].equals("T")) {
                        // This SNP is a SV because:
                        // 1) It has more than one point: e.g "AAA.." or "ACG..."
                        // 2) it is a deletion: e.g "-"
                        excludeSNP = true;
                    }
                    if (!splitted[a2Index].equals("A") && !splitted[a2Index].equals("C") && !splitted[a2Index].equals("G")
                            && !splitted[a2Index].equals("T")) {
                        // This SNP is a SV because:
                        // 1) It has more than one point: e.g "AAA.." or "ACG..."
                        // 2) it is a deletion: e.g "-"
                        excludeSNP = true;
                    }
                }

                String allele = splitted[a1Index] + splitted[a2Index];
                if (exclCgatFlag.equals("YES")) {
                    // Then we have to see if allele is AT TA GC CG to exclude it.
                    if (allele.equals("AT") || allele.equals("TA") || allele.equals("GC") || allele.equals("CG")) {
                        excludeSNP = true;
                    }
                }

                if (excludeSNP) {
                    // Store the positon in the excludeList
                    excludeList.add(positionS);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Finally we put the excludedList into the outputFile
        // We have to create the outputFile:
        // We verify that a file with the same name does not exist. (It should not exist!!)
        File outputFile = new File(excludedSnpsFile);
        try {
            outputFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            for (int i = 0; i < excludeList.size(); i++) {
                writer.write(excludeList.get(i) + "\n");
            }
            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] createListOfExcludedSnps startTime: " + startTime);
            System.out.println("\n[DEBUG] createListOfExcludedSnps endTime: " + stopTime);
            System.out.println("\n[DEBUG] createListOfExcludedSnps elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of createListOfExcludedSnps.");
        }
    }

    /**
     * Method to execute phasing task where input files are in BED format
     * 
     * @param chromo
     * @param bedFile
     * @param bimFile
     * @param famFile
     * @param gmapFile
     * @param shapeitHapsFile
     * @param shapeitSampleFile
     * @param shapeitLogFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void phasingBed(String chromo, String bedFile, String bimFile, String famFile, String gmapFile, String shapeitHapsFile,
            String shapeitSampleFile, String shapeitLogFile, String cmdToStore) throws GuidanceTaskException {

        String shapeitBinary = System.getenv("SHAPEITBINARY");
        if (shapeitBinary == null) {
            throw new GuidanceTaskException("[phasingBed] Error, SHAPEITBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running phasing with parameters:");
            System.out.println("[DEBUG] \t- shapeitBinary            : " + shapeitBinary);
            System.out.println("[DEBUG] \t- Input bedFile            : " + bedFile);
            System.out.println("[DEBUG] \t- Input bimFile            : " + bimFile);
            System.out.println("[DEBUG] \t- Input famFile            : " + famFile);
            System.out.println("[DEBUG] \t- Input gmapFile           : " + gmapFile);
            System.out.println("[DEBUG] \t- Output shapeitHapsFile   : " + shapeitHapsFile);
            System.out.println("[DEBUG] \t- Output shapeitSampleFile : " + shapeitSampleFile);
            System.out.println("[DEBUG] \t- Output shapeitLogFile    : " + shapeitLogFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            // System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // Now we create shapeitHapsFileGz
        String shapeitHapsFileGz = shapeitHapsFile + ".gz";

        String cmd = null;

        if (chromo.equals("23")) {
            cmd = shapeitBinary + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map " + gmapFile
                    + " --chrX --output-max " + shapeitHapsFileGz + " " + shapeitSampleFile
                    + " --thread 47 --effective-size 20000 --output-log " + shapeitLogFile;
        } else {
            cmd = shapeitBinary + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map " + gmapFile + " --output-max "
                    + shapeitHapsFileGz + " " + shapeitSampleFile + " --thread 47 --effective-size 20000 --output-log " + shapeitLogFile;
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Command: " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process shapeitProc;
        try {
            shapeitProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        // Ugly issue: If we run shapeit_v1, all the stdXXX is done stderr, and there is not stdout
        // Ugly issue: If we run shapeit_v2, all the stdXXX is done stdout, and there is not stderr
        try {
            readOutputAndError(shapeitProc.getInputStream(), shapeitHapsFile + ".stdout", shapeitProc.getErrorStream(),
                    shapeitHapsFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = shapeitProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            System.err.println("[phasingBed] Warning executing shapeitProc job, exit value is: " + exitValue);
            System.err.println("[phasingBed]                         (This warning is not fatal).");
            // throw new Exception("Error executing shapeitProc job, exit value is: " + exitValue);
        }

        // Ugly, because shapeit_v2 automatically puts the .log to the file.
        // If there is not output in the impute process. Then we have to create some empty outputs.
        String tmpFile = shapeitLogFile + ".log";
        File tmpShapeitLogFile = new File(tmpFile);
        if (tmpShapeitLogFile.exists()) {
            tmpShapeitLogFile.renameTo(new File(shapeitLogFile));
        }

        // Now we rename shapeitHapsFileGz to shapeitHapsFile
        /*
         * File myHapsGz = new File(shapeitHapsFileGz); File myHaps = new File(shapeitHapsFile); if(myHaps.exists()) {
         * throw new java.io.IOException("[phasingBed] Error, file " + myHaps + " exist!!"); } boolean success =
         * myHapsGz.renameTo(myHaps); if (!success) { throw new java.io.IOException("[phasingBed] Error, file " + myHaps
         * + " could not be renamed to + " + myHapsGz); }
         */

        // Now we rename shapeitHapsFileGz to shapeitHapsFile
        File source = new File(shapeitHapsFileGz);
        File dest = new File(shapeitHapsFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] phasing startTime: " + startTime);
            System.out.println("\n[DEBUG] phasing endTime: " + stopTime);
            System.out.println("\n[DEBUG] phasing elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of phasing.");
        }
    }

    /**
     * Method to execute phasing where input files are in GEN format
     * 
     * @param chromo
     * @param inputGenFile
     * @param inputSampleFile
     * @param gmapFile
     * @param shapeitHapsFile
     * @param shapeitSampleFile
     * @param shapeitLogFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void phasing(String chromo, String inputGenFile, String inputSampleFile, String gmapFile, String shapeitHapsFile,
            String shapeitSampleFile, String shapeitLogFile, String cmdToStore) throws GuidanceTaskException {

        String shapeitBinary = System.getenv("SHAPEITBINARY");
        if (shapeitBinary == null) {
            throw new GuidanceTaskException("[phasing] Error, SHAPEITBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running phasing with parameters:");
            System.out.println("[DEBUG] \t- shapeitBinary            : " + shapeitBinary);
            System.out.println("[DEBUG] \t- Input inputGenFile       : " + inputGenFile);
            System.out.println("[DEBUG] \t- Input inputSampleFile    : " + inputSampleFile);
            System.out.println("[DEBUG] \t- Input gmapFile           : " + gmapFile);
            System.out.println("[DEBUG] \t- Output shapeitHapsFile   : " + shapeitHapsFile);
            System.out.println("[DEBUG] \t- Output shapeitSampleFile : " + shapeitSampleFile);
            System.out.println("[DEBUG] \t- Output shapeitLogFile    : " + shapeitLogFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            // System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // Check if the file is genfile is gzip or not.
        // This is done by checking the magic number of gzip files, which is 0x1f8b (the first two bytes)

        /*
         * String tmpGtoolGenFile = null; tmpGtoolGenFile = inputGenFile;
         * 
         * File tmpInputFile = new File(inputGenFile); RandomAccessFile raf = new RandomAccessFile(tmpInputFile, "r");
         * long n = raf.readInt(); n = n & 0xFFFF0000; raf.close();
         * 
         * if (n == 0x1f8b0000) { System.out.println("[Shapeit] It seems the file " + inputGenFile +
         * " is a gzip file. Magic Number is " + String.format("%x", n)); // Then,use a renamed file tmpGtoolGenFile =
         * inputGenFile + ".gz"; //Now we rename shapeitHapsFileGz to shapeitHapsFile File source = new
         * File(inputGenFile); File dest = new File(tmpGtoolGenFile); copyFile(source, dest); }
         */

        // Now we creat shapeitHapsFileGz
        String shapeitHapsFileGz = shapeitHapsFile + ".gz";

        String cmd = null;
        if (chromo.equals("23")) {
            cmd = shapeitBinary + " --input-gen " + inputGenFile + " " + inputSampleFile + " --input-map " + gmapFile
                    + " --chrX --output-max " + shapeitHapsFileGz + " " + shapeitSampleFile
                    + " --thread 47 --effective-size 20000 --output-log " + shapeitLogFile;
        } else {
            cmd = shapeitBinary + " --input-gen " + inputGenFile + " " + inputSampleFile + " --input-map " + gmapFile + " --output-max "
                    + shapeitHapsFileGz + " " + shapeitSampleFile + " --thread 47 --effective-size 20000 --output-log " + shapeitLogFile;
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Command: " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process shapeitProc;
        try {
            shapeitProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(shapeitProc.getInputStream(), shapeitHapsFile + ".stdout", shapeitProc.getErrorStream(),
                    shapeitHapsFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = shapeitProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            System.err.println("[phasing] Warning executing shapeitProc job, exit value is: " + exitValue);
            System.err.println("[phasing]                         (This warning is not fatal).");
            // throw new Exception("Error executing shapeitProc job, exit value is: " + exitValue);
        }

        // Ugly, because shapeit_v2 automatically puts the .log to the file.
        // If there is not output in the impute process. Then we have to create some empty outputs.
        String tmpFile = shapeitLogFile + ".log";
        File tmpShapeitLogFile = new File(tmpFile);
        if (tmpShapeitLogFile.exists()) {
            tmpShapeitLogFile.renameTo(new File(shapeitLogFile));
        }

        // Now we rename shapeitHapsFileGz to shapeitHapsFile
        File source = new File(shapeitHapsFileGz);
        File dest = new File(shapeitHapsFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] phasing startTime: " + startTime);
            System.out.println("\n[DEBUG] phasing endTime: " + stopTime);
            System.out.println("\n[DEBUG] phasing elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of phasing.");
        }
    }

    /**
     * Method to execute filterHaplotypes where input files are in GEN format
     * 
     * @param hapsFile
     * @param sampleFile
     * @param excludedSnpsFile
     * @param filteredHapsFile
     * @param filteredSampleFile
     * @param filteredLogFile
     * @param filteredHapsVcfFile
     * @param listOfSnpsFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void filterHaplotypes(String hapsFile, String sampleFile, String excludedSnpsFile, String filteredHapsFile,
            String filteredSampleFile, String filteredLogFile, String filteredHapsVcfFile, String listOfSnpsFile, String cmdToStore)
            throws GuidanceTaskException {

        String shapeitBinary = System.getenv("SHAPEITBINARY");
        if (shapeitBinary == null) {
            throw new GuidanceTaskException("[phasing] Error, SHAPEITBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running filterHaplotypes with parameters:");
            System.out.println("[DEBUG] \t- shapeitBinary             : " + shapeitBinary);
            System.out.println("[DEBUG] \t- Input hapsFile            : " + hapsFile);
            System.out.println("[DEBUG] \t- Input sampleFile          : " + sampleFile);
            System.out.println("[DEBUG] \t- Input excludedSnpsFile    : " + excludedSnpsFile);
            System.out.println("[DEBUG] \t- Output filteredHapsFile   : " + filteredHapsFile);
            System.out.println("[DEBUG] \t- Output filteredSampleFile : " + filteredSampleFile);
            System.out.println("[DEBUG] \t- Output filteredLogFile     : " + filteredLogFile);
            System.out.println("[DEBUG] \t- Output filteredHapsVcfFile : " + filteredHapsVcfFile);
            System.out.println("[DEBUG] \t- Output lisfOfSnpsFile     : " + listOfSnpsFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            // System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // We have to make sure whether we are using renamed files of the original gz files.
        // We detect this situation by scanning the last three characters:
        String extension = hapsFile.substring(Math.max(0, hapsFile.length() - 3));
        // System.out.println("DEBUG \t The file extension is: " + extension + " and the file is " + shapeitHapsFile);

        String hapsFileGz = null;
        if (extension.equals(".gz")) {
            hapsFileGz = hapsFile;
        } else {
            // If hapsFile exists, then hapsFileGz exists also.
            hapsFileGz = hapsFile + ".gz";
            // String imputeFileGz = imputeFile + ".gz";
        }

        // Check if the file is genfile is gzip or not.
        // This is done by checking the magic number of gzip files, which is 0x1f8b (the first two bytes)

        /*
         * String tmpGtoolGenFile = null; tmpGtoolGenFile = inputGenFile;
         * 
         * File tmpInputFile = new File(inputGenFile); RandomAccessFile raf = new RandomAccessFile(tmpInputFile, "r");
         * long n = raf.readInt(); n = n & 0xFFFF0000; raf.close();
         * 
         * if (n == 0x1f8b0000) { System.out.println("[Shapeit] It seems the file " + inputGenFile +
         * " is a gzip file. Magic Number is " + String.format("%x", n)); // Then,use a renamed file tmpGtoolGenFile =
         * inputGenFile + ".gz"; //Now we rename shapeitHapsFileGz to shapeitHapsFile File source = new
         * File(inputGenFile); File dest = new File(tmpGtoolGenFile); copyFile(source, dest); }
         */

        // Now we creat hapsFileGz
        String filteredHapsFileGz = filteredHapsFile + ".gz";
        String filteredHapsVcfFileGz = filteredHapsVcfFile + ".gz";

        String cmd = null;
        cmd = shapeitBinary + " -convert --input-haps " + hapsFileGz + " " + sampleFile + " --exclude-snp " + excludedSnpsFile + " "
                + " --output-haps " + filteredHapsFileGz + " " + filteredSampleFile + " --output-log " + filteredLogFile + " --output-vcf "
                + filteredHapsVcfFileGz;

        if (DEBUG) {
            System.out.println("\n[DEBUG] Command: " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process shapeitProc;
        try {
            shapeitProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        // Ugly issue: If we run shapeit_v1, all the stdXXX is done stderr, and there is not stdout
        // Ugly issue: If we run shapeit_v2, all the stdXXX is done stdout, and there is not stderr
        try {
            readOutputAndError(shapeitProc.getInputStream(), filteredHapsFile + ".stdout", shapeitProc.getErrorStream(),
                    filteredHapsFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = shapeitProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            System.err.println("[filterHaplotypes] Warning executing shapeitProc job in mode -convert, exit value is: " + exitValue);
            System.err.println("[filterHaplotypes]                         (This warning is not fatal).");
            // throw new Exception("Error executing shapeitProc job, exit value is: " + exitValue);
        }

        // Ugly, because shapeit_v2 automatically puts the .log to the file.
        // If there is not output in the impute process. Then we have to create some empty outputs.
        String tmpFile = filteredLogFile + ".log";
        File tmpFilteredLogFile = new File(tmpFile);
        if (tmpFilteredLogFile.exists()) {
            tmpFilteredLogFile.renameTo(new File(filteredLogFile));
        }

        System.err.println("[filterHaplotypes] Filtering haplotypes OK. Now we create the listofSnps...");

        // Now we have to create the list of snps and write them into the output file.
        // Taking into account that shapeit had generated a gziped file, then we have to use GZIPInputStream.
        File outFilteredFile = new File(listOfSnpsFile);

        // Trys to create the file
        boolean bool = false;
        try {
            bool = outFilteredFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about de existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + listOfSnpsFile + " was succesfuly created? " + bool);

        try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(filteredHapsFileGz));
                Reader decoder = new InputStreamReader(inputGz);
                BufferedReader br = new BufferedReader(decoder);
                BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(outFilteredFile))) {
            String line = null;
            while ((line = br.readLine()) != null) {
                String[] splittedLine = line.split(" ");
                String rsId = splittedLine[1];

                writerFiltered.write(rsId);
                writerFiltered.newLine();
            }

            writerFiltered.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Now we rename filteredHapsFileGz to filteredHapsFile
        File source = new File(filteredHapsFileGz);
        File dest = new File(filteredHapsFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        source = new File(filteredHapsVcfFileGz);
        dest = new File(filteredHapsVcfFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] filterHaplotypes startTime: " + startTime);
            System.out.println("\n[DEBUG] filterHaplotypes endTime: " + stopTime);
            System.out.println("\n[DEBUG] filterHaplotypes elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of filterHaplotypes.");
        }
    }

    /**
     * Method to impute with impute
     * 
     * @param gmapFile
     * @param knownHapFile
     * @param legendFile
     * @param shapeitHapsFile
     * @param shapeitSampleFile
     * @param lim1S
     * @param lim2S
     * @param pairsFile
     * @param imputeFile
     * @param imputeFileInfo
     * @param imputeFileSummary
     * @param imputeFileWarnings
     * @param theChromo
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void imputeWithImpute(String gmapFile, String knownHapFile, String legendFile, String shapeitHapsFile,
            String shapeitSampleFile, String lim1S, String lim2S, String pairsFile, String imputeFile, String imputeFileInfo,
            String imputeFileSummary, String imputeFileWarnings, String theChromo, String cmdToStore) throws GuidanceTaskException {

        String impute2Binary = System.getenv("IMPUTE2BINARY");
        if (impute2Binary == null) {
            throw new GuidanceTaskException("[impute] Error, IMPUTE2BINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running impute with parameters:");
            System.out.println("[DEBUG] \t- impute2Binary             : " + impute2Binary);
            System.out.println("[DEBUG] \t- Input gmapFile            : " + gmapFile);
            System.out.println("[DEBUG] \t- Input knownHapFile        : " + knownHapFile);
            System.out.println("[DEBUG] \t- Input legendHapFile       : " + legendFile);
            System.out.println("[DEBUG] \t- Input shapeitHapsFile     : " + shapeitHapsFile);
            System.out.println("[DEBUG] \t- Input shapeitSampleFile   : " + shapeitSampleFile);
            System.out.println("[DEBUG] \t- Input lim1S               : " + lim1S);
            System.out.println("[DEBUG] \t- Input lim2S               : " + lim2S);
            System.out.println("[DEBUG] \t- Input pairsFile           : " + pairsFile);
            System.out.println("[DEBUG] \t- Output imputeFile         : " + imputeFile);
            System.out.println("[DEBUG] \t- Output imputeFileInfo     : " + imputeFileInfo);
            System.out.println("[DEBUG] \t- Output imputeFileSummary  : " + imputeFileSummary);
            System.out.println("[DEBUG] \t- Output imputeFileWarnings : " + imputeFileWarnings);
            System.out.println("[DEBUG] \t- Input  theChromo          : " + theChromo);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            // System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // We have to make sure whether we are using renamed files of the original gz files.
        // We detect this situation by scanning the last three characters:
        String extension = shapeitHapsFile.substring(Math.max(0, shapeitHapsFile.length() - 3));
        // System.out.println("DEBUG \t The file extension is: " + extension + " and the file is " + shapeitHapsFile);

        String shapeitHapsFileGz = null;
        if (extension.equals(".gz")) {
            shapeitHapsFileGz = shapeitHapsFile;
        } else {
            // If shapeitHapsFile exists, then shapeitHapsFileGz exists also.
            shapeitHapsFileGz = shapeitHapsFile + ".gz";
            // String imputeFileGz = imputeFile + ".gz";
        }

        String cmd = null;

        if (theChromo.equals("23")) {
            cmd = impute2Binary + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l " + legendFile + " -known_haps_g "
                    + shapeitHapsFileGz + " -sample_g " + shapeitSampleFile + " -int " + lim1S + " " + lim2S + "  -chrX -exclude_snps_g "
                    + pairsFile + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r " + imputeFileSummary
                    + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz";
        } else {
            cmd = impute2Binary + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l " + legendFile + " -known_haps_g "
                    + shapeitHapsFileGz + " -int " + lim1S + " " + lim2S + " -exclude_snps_g " + pairsFile
                    + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r " + imputeFileSummary + " -w "
                    + imputeFileWarnings + " -no_sample_qc_info -o_gz";
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Command: " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process imputeProc;
        try {
            imputeProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(imputeProc.getInputStream(), imputeFile + ".stdout", imputeProc.getErrorStream(), imputeFile + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = imputeProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            System.err.println("[impute] Warning executing imputeProc job, exit value is: " + exitValue);
            System.err.println("                        (This warning is not fatal).");
            // throw new Exception("Error executing imputeProc job, exit value is: " + exitValue);
        }

        // With the -o_gz option in the comand, the outputs are imputeFile.gz
        // If there is not output in the impute process. Then we have to create some empty outputs.
        File fImputeGz = new File(imputeFile + ".gz");
        File fImpute = new File(imputeFile);
        if (!fImputeGz.exists()) {
            // System.out.println("[impute] The file " + imputeFile + ".gz" + "does not exist, then, we create it.. " +
            // imputeFile);
            try {
                fImpute.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            gzipFile(imputeFile, imputeFile + ".gz");
        }
        File source = new File(imputeFile + ".gz");
        File dest = new File(imputeFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        File fImputeInfo = new File(imputeFileInfo);
        if (!fImputeInfo.exists()) {
            try {
                fImputeInfo.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            /* do something */
        }

        File fImputeSummary = new File(imputeFileSummary);
        if (!fImputeSummary.exists()) {
            try {
                fImputeSummary.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            /* do something */
        }

        File fImputeWarnings = new File(imputeFileWarnings);
        if (!fImputeWarnings.exists()) {
            try {
                fImputeWarnings.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            /* do something */
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] imputeWithImpute startTime: " + startTime);
            System.out.println("\n[DEBUG] imputeWithImpute endTime: " + stopTime);
            System.out.println("\n[DEBUG] imputeWithImpute elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of imputeWithImpute with parameters:");
        }

    }

    /**
     * Method to impute with minimac
     * 
     * @param knownHapFile
     * @param filteredHapsFile
     * @param filteredSampleFile
     * @param filteredListOfSnpsFile
     * @param imputedMMFileName
     * @param imputedMMInfoFile
     * @param imputedMMErateFile
     * @param imputedMMRecFile
     * @param imputedMMDoseFile
     * @param imputedMMLogFile
     * @param theChromo
     * @param lim1S
     * @param lim2S
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void imputeWithMinimac(String knownHapFile, String filteredHapsFile, String filteredSampleFile,
            String filteredListOfSnpsFile, String imputedMMFileName, String imputedMMInfoFile,
            // String imputedMMDraftFile,
            String imputedMMErateFile, String imputedMMRecFile, String imputedMMDoseFile, String imputedMMLogFile, String theChromo,
            String lim1S, String lim2S, String cmdToStore) throws GuidanceTaskException {

        String minimacBinary = System.getenv("MINIMACBINARY");
        if (minimacBinary == null) {
            throw new GuidanceTaskException("[impute] Error, MINIMACBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running imputation with parameters:");
            System.out.println("[DEBUG] \t- minimacBinary                : " + minimacBinary);
            System.out.println("[DEBUG] \t- Input knownHapFile           : " + knownHapFile);
            System.out.println("[DEBUG] \t- Input filteredHapsFile       : " + filteredHapsFile);
            System.out.println("[DEBUG] \t- Input filteredSampleFile     : " + filteredSampleFile);
            System.out.println("[DEBUG] \t- Input filteredListOfSnpsFile : " + filteredListOfSnpsFile);
            System.out.println("[DEBUG] \t- Output imputedMMFileName     : " + imputedMMFileName);
            System.out.println("[DEBUG] \t- Output imputedMMInfoFile     : " + imputedMMInfoFile);
            // System.out.println("[DEBUG] \t- Output imputedMMDraftFile : " + imputedMMDraftFile);
            System.out.println("[DEBUG] \t- Output imputedMMErateFile    : " + imputedMMErateFile);
            System.out.println("[DEBUG] \t- Output imputedMMRecFile      : " + imputedMMRecFile);
            System.out.println("[DEBUG] \t- Output imputedMMDoseFile     : " + imputedMMDoseFile);
            System.out.println("[DEBUG] \t- Output imputedMMLogFile      : " + imputedMMLogFile);
            System.out.println("[DEBUG] \t- Input  theChromo             : " + theChromo);
            System.out.println("[DEBUG] \t- Input lim1S                  : " + lim1S);
            System.out.println("[DEBUG] \t- Input lim2S                  : " + lim2S);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            // System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // We have to make sure whether we are using renamed files of the original gz files.
        // We detect this situation by scanning the last three characters:
        String extension = filteredHapsFile.substring(Math.max(0, filteredHapsFile.length() - 3));
        System.out.println("DEBUG \t The file extension is: " + extension + " and the file is " + filteredHapsFile);

        String filteredHapsFileGz = null;
        if (extension.equals(".gz")) {
            filteredHapsFileGz = filteredHapsFile;
        } else {
            // If filteredHapsFile exists, then filteredHapsFileGz exists also.
            filteredHapsFileGz = filteredHapsFile + ".gz";
        }

        String cmd = null;

        if (theChromo.equals("23")) {
            cmd = minimacBinary + " --vcfReference --refHaps " + knownHapFile + " --snps " + filteredListOfSnpsFile + " --shape_haps "
                    + filteredHapsFileGz + " --sample " + filteredSampleFile + " --vcfstart " + lim1S + " --vcfend " + lim2S + " --chr "
                    + theChromo + " --vcfwindow 250000 --rounds 5 --states 200 --prefix " + imputedMMFileName + " --gzip";
        } else {
            cmd = minimacBinary + " --vcfReference --refHaps " + knownHapFile + " --snps " + filteredListOfSnpsFile + " --shape_haps "
                    + filteredHapsFileGz + " --sample " + filteredSampleFile + " --vcfstart " + lim1S + " --vcfend " + lim2S + " --chr "
                    + theChromo + " --vcfwindow 250000 --rounds 5 --states 200 --prefix " + imputedMMFileName + " --gzip";
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Command: " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process minimacProc;
        try {
            minimacProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Handling the streams so that dead lock situation never occurs
        try {
            readOutputAndError(minimacProc.getInputStream(), imputedMMFileName + ".stdout", minimacProc.getErrorStream(),
                    imputedMMFileName + ".stderr");
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = minimacProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            System.err.println("[impute] Warning executing minimacProc job, exit value is: " + exitValue);
            System.err.println("                        (This warning is not fatal).");
            // throw new Exception("Error executing imputeProc job, exit value is: " + exitValue);
        }

        // With the -o_gz option in the comand, the outputs are imputeFile.gz
        // If there is not output in the impute process. Then we have to create some empty outputs.
        File infoFileGz = new File(imputedMMFileName + ".info.gz");
        File infoFile = new File(imputedMMInfoFile);
        if (!infoFileGz.exists()) {
            // System.out.println("[impute] The file " + imputeFile + ".gz" + "does not exist, then, we create it.. " +
            // imputeFile);
            try {
                infoFile.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            gzipFile(imputedMMInfoFile, imputedMMInfoFile + ".info.gz");
        }
        File source = new File(imputedMMFileName + ".info.gz");
        File dest = new File(imputedMMInfoFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // File draftFileGz = new File(imputedMMFileName + ".info.draft.gz");
        // File draftFile = new File(imputedMMDraftFile);
        // if(!draftFileGz.exists()) {
        // //System.out.println("[impute] The file " + imputeFile + ".gz" + "does not exist, then, we create it.. " +
        // imputeFile);
        // draftFile.createNewFile();
        // gzipFile(imputedMMDraftFile, imputedMMDraftFile + ".info.draft.gz");
        // }
        // source = new File(imputedMMFileName + ".info.draft.gz");
        // dest = new File(imputedMMDraftFile);
        // copyFile(source, dest);

        File erateFileGz = new File(imputedMMFileName + ".erate.gz");
        File erateFile = new File(imputedMMErateFile);
        if (!erateFileGz.exists()) {
            // System.out.println("[impute] The file " + imputeFile + ".gz" + "does not exist, then, we create it.. " +
            // imputeFile);
            try {
                erateFile.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            gzipFile(imputedMMErateFile, imputedMMErateFile + ".erate.gz");
        }
        source = new File(imputedMMFileName + ".erate.gz");
        dest = new File(imputedMMErateFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        File recFileGz = new File(imputedMMFileName + ".rec.gz");
        File recFile = new File(imputedMMRecFile);
        if (!recFileGz.exists()) {
            // System.out.println("[impute] The file " + imputeFile + ".gz" + "does not exist, then, we create it.. " +
            // imputeFile);
            try {
                recFile.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            gzipFile(imputedMMRecFile, imputedMMRecFile + ".rec.gz");
        }
        source = new File(imputedMMFileName + ".rec.gz");
        dest = new File(imputedMMRecFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        File doseFileGz = new File(imputedMMFileName + ".dose.gz");
        File doseFile = new File(imputedMMDoseFile);
        if (!doseFileGz.exists()) {
            // System.out.println("[impute] The file " + imputeFile + ".gz" + "does not exist, then, we create it.. " +
            // imputeFile);
            try {
                doseFile.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            gzipFile(imputedMMDoseFile, imputedMMDoseFile + ".dose.gz");
        }
        source = new File(imputedMMFileName + ".dose.gz");
        dest = new File(imputedMMDoseFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        source = new File(imputedMMFileName + ".stdout");
        dest = new File(imputedMMLogFile);
        try {
            copyFile(source, dest);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] imputeWithMinimac startTime: " + startTime);
            System.out.println("\n[DEBUG] imputeWithMinimac endTime: " + stopTime);
            System.out.println("\n[DEBUG] imputeWithMinimac elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of imputeWithMinimac with parameters:");
        }

    }

    /**
     * Method to filter by info
     * 
     * @param imputeFileInfo
     * @param filteredFile
     * @param threshold
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void filterByInfo(String imputeFileInfo, String filteredFile, String threshold, String cmdToStore)
            throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running filterByInfo with parameters:");
            System.out.println("[DEBUG] \t- Input inputeFileInfo   : " + imputeFileInfo);
            System.out.println("[DEBUG] \t- Output filteredFile    : " + filteredFile);
            System.out.println("[DEBUG] \t- Input threshold        : " + threshold);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // The position of info and rsId values in the imputeFileInfo
        int infoIndex = 6;
        int rsIdIndex = 1;

        // Convert threshold string into thresholdDouble
        Double thresholdDouble = Double.parseDouble(threshold); // store info value in Double format

        File outFilteredFile = new File(filteredFile);
        // Trys to create the file
        boolean bool = false;
        try {
            bool = outFilteredFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about de existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + filteredFile + " was succesfuly created? " + bool);

        try (FileReader fr = new FileReader(imputeFileInfo);
                BufferedReader br = new BufferedReader(fr);
                BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(outFilteredFile))) {

            // We read each line of the imputeFileInfo and put them into string.
            // I read the header
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] splittedLine = line.split(" ");// delimiter I assume single space.
                Double info = Double.parseDouble(splittedLine[infoIndex]); // store info value in Double format

                // Store rsID into filteredFile if info >= threshold
                int retval = Double.compare(info, thresholdDouble);
                if (retval >= 0) {
                    // The info value is greater or equal to the threshold, then store the rsID into the output file.
                    writerFiltered.write(splittedLine[rsIdIndex]);
                    writerFiltered.newLine();
                }
            }

            writerFiltered.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] filterByInfo startTime: " + startTime);
            System.out.println("\n[DEBUG] filterByInfo endTime: " + stopTime);
            System.out.println("\n[DEBUG] filterByInfo elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of filterByInfo");
        }

    }

    /**
     * Method to filter by all
     * 
     * @param inputFile
     * @param outputFile
     * @param outputCondensedFile
     * @param mafThresholdS
     * @param infoThresholdS
     * @param hweCohortThresholdS
     * @param hweCasesThresholdS
     * @param hweControlsThresholdS
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void filterByAll(String inputFile, String outputFile, String outputCondensedFile, String mafThresholdS,
            String infoThresholdS, String hweCohortThresholdS, String hweCasesThresholdS, String hweControlsThresholdS, String cmdToStore)
            throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running filterByAll with parameters:");
            System.out.println("[DEBUG] \t- Input summaryFile             : " + inputFile);
            System.out.println("[DEBUG] \t- Output outputFile             : " + outputFile);
            System.out.println("[DEBUG] \t- Output outputCondensedFile    : " + outputCondensedFile);
            System.out.println("[DEBUG] \t- Input maf threshold           : " + mafThresholdS);
            System.out.println("[DEBUG] \t- Input info threshold          : " + infoThresholdS);
            System.out.println("[DEBUG] \t- Input hwe cohort threshold    : " + hweCohortThresholdS);
            System.out.println("[DEBUG] \t- Input hwe controls threshold  : " + hweCasesThresholdS);
            System.out.println("[DEBUG] \t- Input hwe cases threshold     : " + hweControlsThresholdS);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // Convert threshold string into thresholdDouble
        Double mafThreshold = Double.parseDouble(mafThresholdS);
        Double infoThreshold = Double.parseDouble(infoThresholdS);
        Double hweCohortThreshold = Double.parseDouble(hweCohortThresholdS);
        Double hweCasesThreshold = Double.parseDouble(hweCasesThresholdS);
        Double hweControlsThreshold = Double.parseDouble(hweControlsThresholdS);

        File outFilteredFile = new File(outputFile);
        File outCondensedFile = new File(outputCondensedFile);

        // Trys to create the file
        boolean bool = false;
        try {
            outFilteredFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about de existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outputFile + " was succesfuly created? " + bool);

        // Print information about de existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outputCondensedFile + " was succesfuly created? " + bool);

        Hashtable<String, Integer> inputFileHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> inputFileHashTableIndexReversed = new Hashtable<>();

        try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(inputFile));
                Reader decoder = new InputStreamReader(inputGz);
                BufferedReader br = new BufferedReader(decoder);
                BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(outFilteredFile));
                BufferedWriter writerCondensed = new BufferedWriter(new FileWriter(outCondensedFile))) {

            // I read the header
            String line = br.readLine();
            // Put the header in the output file.
            writerFiltered.write(line);
            writerFiltered.newLine();

            inputFileHashTableIndex = createHashWithHeader(line, "\t");
            inputFileHashTableIndexReversed = createHashWithHeaderReversed(line, "\t");

            String headerCondensed = "chr\tposition\talleleA\talleleB\tpvalue\tinfo_all";
            writerCondensed.write(headerCondensed);
            writerCondensed.newLine();

            while ((line = br.readLine()) != null) {
                String[] splittedLine = line.split("\t");// delimiter I assume single space.

                String chromo = splittedLine[inputFileHashTableIndex.get("chr")];
                String infoS = splittedLine[inputFileHashTableIndex.get("info_all")];
                String alleleA = splittedLine[inputFileHashTableIndex.get("alleleA")];
                String alleleB = splittedLine[inputFileHashTableIndex.get("alleleB")];

                // We start with these values for hwe values just to allows the X chromosome to pass the if statement of
                // the
                // next lines
                // Just remember that hwe filtering when chromo X is being processed does not make sense.
                String hwe_cohortS = "1.0";
                String hwe_casesS = "1.0";
                String hwe_controlsS = "1.0";

                if (!chromo.equals("23")) {
                    hwe_cohortS = splittedLine[inputFileHashTableIndex.get("cohort_1_hwe")];
                    hwe_casesS = splittedLine[inputFileHashTableIndex.get("cases_hwe")];
                    hwe_controlsS = splittedLine[inputFileHashTableIndex.get("controls_hwe")];
                }

                String cases_mafS = splittedLine[inputFileHashTableIndex.get("cases_maf")];
                String controls_mafS = splittedLine[inputFileHashTableIndex.get("controls_maf")];

                String position = splittedLine[inputFileHashTableIndex.get("position")];
                // String beta = splittedLine[inputFileHashTableIndex.get("frequentist_add_beta_1")];
                // String se = splittedLine[inputFileHashTableIndex.get("frequentist_add_se_1")];
                String pva = splittedLine[inputFileHashTableIndex.get("frequentist_add_pvalue")];

                String chrbpb = chromo + "\t" + position + "\t" + alleleA + "\t" + alleleB + "\t" + pva + "\t" + infoS;

                if (!cases_mafS.equals("NA") && !controls_mafS.equals("NA") && !infoS.equals("NA") && !hwe_cohortS.equals("NA")
                        && !hwe_casesS.equals("NA") && !hwe_controlsS.equals("NA") && !pva.equals("NA")) {
                    Double cases_maf = Double.parseDouble(cases_mafS);
                    Double controls_maf = Double.parseDouble(controls_mafS);
                    Double info = Double.parseDouble(infoS);
                    Double hweCohort = 1.0;
                    Double hweCases = 1.0;
                    Double hweControls = 1.0;

                    if (!chromo.equals("23")) {
                        hweCohort = Double.parseDouble(hwe_cohortS);
                        hweCases = Double.parseDouble(hwe_casesS);
                        hweControls = Double.parseDouble(hwe_controlsS);
                    }

                    if (cases_maf >= mafThreshold && controls_maf >= mafThreshold && info >= infoThreshold && // VERIFICAR
                                                                                                              // LA
                                                                                                              // CONDICION
                            hweCohort >= hweCohortThreshold && hweCases >= hweCasesThreshold && hweControls >= hweControlsThreshold) {

                        writerFiltered.write(line);
                        writerFiltered.newLine();

                        writerCondensed.write(chrbpb);
                        writerCondensed.newLine();
                    }
                }
            }

            writerFiltered.flush();
            writerCondensed.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(outputFile, outputFile + ".gz");
        File fc = new File(outputFile);
        File fGz = new File(outputFile + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        gzipFile(outputCondensedFile, outputCondensedFile + ".gz");
        fc = new File(outputCondensedFile);
        fGz = new File(outputCondensedFile + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] filterByAll startTime: " + startTime);
            System.out.println("\n[DEBUG] filterByAll endTime: " + stopTime);
            System.out.println("\n[DEBUG] filterByAll elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of filterByAll");
        }

    }

    /**
     * Method to joint filtered by all files
     * 
     * @param filteredByAllA
     * @param filteredByAllB
     * @param filteredByAllC
     * @param rpanelName
     * @param rpanelFlag
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void jointFilteredByAllFiles(String filteredByAllA, String filteredByAllB, String filteredByAllC, String rpanelName,
            String rpanelFlag, String cmdToStore) throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running jointFilteredByAllFiles with parameters:");
            System.out.println("[DEBUG] \t- Input filteredByAllA          : " + filteredByAllA);
            System.out.println("[DEBUG] \t- Input filteredByAllB          : " + filteredByAllB);
            System.out.println("[DEBUG] \t- Output filteredByAllC         : " + filteredByAllC);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);

        }
        long startTime = System.currentTimeMillis();

        File outFilteredByAllFile = new File(filteredByAllC);

        // Try to create the file
        boolean bool = false;
        try {
            bool = outFilteredByAllFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about the existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outFilteredByAllFile + " was succesfuly created? " + bool);

        try (BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(outFilteredByAllFile))) {
            try (GZIPInputStream filteredByAllGz = new GZIPInputStream(new FileInputStream(filteredByAllA));
                    Reader decoder = new InputStreamReader(filteredByAllGz);
                    BufferedReader br = new BufferedReader(decoder)) {

                boolean putRefpanel = false;

                // I read the header
                String line = br.readLine();
                // I put the refpanel column in the header:
                // if( rpanelFlag.equals("YES") ) {
                String[] splittedHeader = line.split("\t");
                if (!splittedHeader[splittedHeader.length - 1].equals("refpanel")) {
                    line = line + "\trefpanel";
                    putRefpanel = true;
                }
                // }

                // Put the header in the output file.
                writerFiltered.write(line);
                writerFiltered.newLine();

                while ((line = br.readLine()) != null) {
                    // if( rpanelFlag.equals("YES") ) {
                    if (putRefpanel == true) {
                        line = line + "\t" + rpanelName;
                    }
                    // }
                    writerFiltered.write(line);
                    writerFiltered.newLine();
                }
            }

            // Do the same with the filteredByAllB file if this is different to the filteredByAllA file
            // OK, I explain now: The only way filteredByAllB = filteredByAllA is when there is only one chromosome to
            // process.
            // In that case, in the main program, we put the same file as filteredByAllA and filteredByAllB.
            if (!filteredByAllA.equals(filteredByAllB)) {
                try (GZIPInputStream filteredByAllGz = new GZIPInputStream(new FileInputStream(filteredByAllB));
                        Reader decoder = new InputStreamReader(filteredByAllGz);
                        BufferedReader br = new BufferedReader(decoder)) {

                    boolean putRefpanel = false;
                    // I read the header and skip it.
                    String line = br.readLine();
                    // if( rpanelFlag.equals("YES") ) {
                    String[] splittedHeader = line.split("\t");
                    if (!splittedHeader[splittedHeader.length - 1].equals("refpanel")) {
                        putRefpanel = true;
                    }

                    while ((line = br.readLine()) != null) {
                        // if( rpanelFlag.equals("YES") ) {
                        if (putRefpanel == true) {
                            line = line + "\t" + rpanelName;
                        }
                        // }
                        writerFiltered.write(line);
                        writerFiltered.newLine();
                    }
                }
            }

            writerFiltered.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(filteredByAllC, filteredByAllC + ".gz");
        File fc = new File(filteredByAllC);
        File fGz = new File(filteredByAllC + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] jointFilteredByAllFiles startTime: " + startTime);
            System.out.println("\n[DEBUG] jointFilteredByAllFiles endTime: " + stopTime);
            System.out.println("\n[DEBUG] jointFilteredByAllFiles elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of jointFilteredByAllFiles");
        }

    }

    /**
     * Method to joint condensed files
     * 
     * @param inputAFile
     * @param inputBFile
     * @param outputFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void jointCondensedFiles(String inputAFile, String inputBFile, String outputFile, String cmdToStore)
            throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running jointCondensedFiles with parameters:");
            System.out.println("[DEBUG] \t- InputAFile                    : " + inputAFile);
            System.out.println("[DEBUG] \t- InputBFile                    : " + inputBFile);
            System.out.println("[DEBUG] \t- Output outputCondensedFile    : " + outputFile);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }

        long startTime = System.currentTimeMillis();

        File outCondensedFile = new File(outputFile);

        // Try to create the file
        // bool = outCondensedFile.createNewFile();

        // Print information about de existence of the file
        // System.out.println("\n[DEBUG] \t- Output file " + outputFile + " was succesfuly created? " + bool);

        try (BufferedWriter writerCondensed = new BufferedWriter(new FileWriter(outCondensedFile))) {

            try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(inputAFile));
                    Reader decoder = new InputStreamReader(inputGz);
                    BufferedReader br = new BufferedReader(decoder)) {

                // I read the header
                String line = br.readLine();
                // Put the header in the output file.
                writerCondensed.write(line);
                writerCondensed.newLine();

                while ((line = br.readLine()) != null) {
                    writerCondensed.write(line);
                    writerCondensed.newLine();
                }
            }

            // Do the same with the inputB file if this is not null.
            // OK, I explain now: The only way inputAFile = inputBFile is when there is only one chromosome to process.
            // In that case, in the main program, we put he same file as inputAFile and inputBFile.
            if (!inputAFile.equals(inputBFile)) {
                // Now the next file: inputBFile
                try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(inputBFile));
                        Reader decoder = new InputStreamReader(inputGz);
                        BufferedReader br = new BufferedReader(decoder)) {

                    // I read the header and skip it.
                    String line = br.readLine();

                    while ((line = br.readLine()) != null) {
                        writerCondensed.write(line);
                        writerCondensed.newLine();
                    }
                }
            }

            writerCondensed.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(outputFile, outputFile + ".gz");
        File fc = new File(outputFile);
        File fGz = new File(outputFile + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] jointCondensedFiles startTime: " + startTime);
            System.out.println("\n[DEBUG] jointCondensedFiles endTime: " + stopTime);
            System.out.println("\n[DEBUG] jointCondensedFiles elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of jointCondensedFiles");
        }
    }

    /**
     * Method to combine panels
     * 
     * @param resultsPanelA
     * @param resultsPanelB
     * @param resultsPanelC
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void combinePanels(String resultsPanelA, String resultsPanelB, String resultsPanelC, String cmdToStore)
            throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running combinePanels with parameters:");
            System.out.println("[DEBUG] \t- resultsPanelA             : " + resultsPanelA);
            System.out.println("[DEBUG] \t- resultsPanelB             : " + resultsPanelB);
            System.out.println("[DEBUG] \t- resultsPanelC             : " + resultsPanelC);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // We read each line of the resultsPanelA and put them into the String
        // IMPORTANT: In that case we sort by position and not by position_rsID. So, maybe we
        // are going to lose some SNPs...
        TreeMap<Integer, String> fileTreeMapA = new TreeMap<>();
        String header = null;
        try (FileReader frA = new FileReader(resultsPanelA); BufferedReader brA = new BufferedReader(frA)) {
            String line = null;
            // First: read the header and avoid it
            header = brA.readLine();

            while ((line = brA.readLine()) != null) {
                String[] splitted = line.split("\t");
                int position = Integer.parseInt(splitted[0]);

                // Now, we put this String into the treemap with the key positionAndRsId
                fileTreeMapA.put(position, line);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // We read each line of the resultsPanelB and put them into fileAList array of Strings
        TreeMap<Integer, String> fileTreeMapB = new TreeMap<>();
        try (FileReader frB = new FileReader(resultsPanelB); BufferedReader brB = new BufferedReader(frB)) {
            // First: read the header and avoid it
            String line = brB.readLine();

            while ((line = brB.readLine()) != null) {
                String[] splitted = line.split("\t");
                int position = Integer.parseInt(splitted[0]);

                // Now, we put this String into the treemap with the key positionAndRsId
                fileTreeMapB.put(position, line);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // A place to store the results of this combining
        TreeMap<Integer, String> fileTreeMapC = new TreeMap<>();
        // We first iterate the fileTreeMapA
        Set<Entry<Integer, String>> mySet = fileTreeMapA.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<Integer, String>> iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<Integer, String> m = iter.next();
            // getKey is used to get key of Map
            // position=(Integer)m.getKey();
            int position = (Integer) m.getKey();
            String lineA = (String) m.getValue();

            if (fileTreeMapB.containsKey(position)) {
                // If the fileTreeMapB contains this key, then we have to choose
                // the ones that have a better info (greater info).
                String[] lineASplitted = lineA.split("\t");
                Double infoA = Double.parseDouble(lineASplitted[2]);

                String lineB = fileTreeMapB.get(position);
                String[] lineBSplitted = lineB.split("\t");
                Double infoB = Double.parseDouble(lineBSplitted[2]);

                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(position, lineA);
                } else {
                    fileTreeMapC.put(position, lineB);
                }
                // Now we can remove this value from the fileTreeMapB
                fileTreeMapB.remove(position);
            } else {
                // Then, the fileTreeMapB does not contain this key.
                // therefore we keep the one in fileTreeMapA.
                fileTreeMapC.put(position, lineA);
            }
        }

        // Now we have to put in fileTreeMapC the rest of values that remain in fileTreeMapB.
        // We iterate the fileTreeMapB (the rest of the...)
        mySet = fileTreeMapB.entrySet();
        // Move next key and value of Map by iterator
        iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<Integer, String> m = iter.next();
            // getKey is used to get key of Map
            int position = (Integer) m.getKey();
            String lineB = (String) m.getValue();
            // Then we have to choose the value in fileTreeMapC
            fileTreeMapC.put(position, lineB);
        }

        // Finally we put the fileTreeMapC into the outputFile
        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File outputFile = new File(resultsPanelC);
        try {
            outputFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            // We print the header which is the same always!.
            writer.write(header);
            writer.newLine();

            mySet = fileTreeMapC.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<Integer, String> m = iter.next();
                // getKey is used to get key of Map
                String myLine = (String) m.getValue();

                writer.write(myLine);
                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] combinePanels startTime: " + startTime);
            System.out.println("\n[DEBUG] combinePanels endTime: " + stopTime);
            System.out.println("\n[DEBUG] combinePanels elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of combinePanels");
        }
    }

    /**
     * Method to combine panels complex
     * 
     * @param resultsPanelA
     * @param resultsPanelB
     * @param resultsPanelC
     * @param lim1
     * @param lim2
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void combinePanelsComplex(String resultsPanelA, String resultsPanelB, String resultsPanelC, int lim1, int lim2,
            String cmdToStore) throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running combinePanelsComplex with parameters:");
            System.out.println("[DEBUG] \t- resultsPanelA             : " + resultsPanelA);
            System.out.println("[DEBUG] \t- resultsPanelB             : " + resultsPanelB);
            System.out.println("[DEBUG] \t- resultsPanelC             : " + resultsPanelC);
            System.out.println("[DEBUG] \t- lim1               : " + lim1);
            System.out.println("[DEBUG] \t- lim2                 : " + lim2);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // int chrStart = Integer.parseInt(chromoStart);
        // int chrEnd = Integer.parseInt(chromoEnd);

        /*
         * int posIdx = 0; int a1Idx = 4; int a2Idx = 5; int chrIdx = 35; int infoIdx = 2;
         */
        int posIdx = 0;
        int a1Idx = 0;
        int a2Idx = 0;
        int chrIdx = 0;
        int infoIdx = 0;

        // First, we uncompress the input files
        gunzipFile(resultsPanelA, resultsPanelA + ".temp");
        gunzipFile(resultsPanelB, resultsPanelB + ".temp");

        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File outputFile = new File(resultsPanelC);
        try {
            outputFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // We read each line of the resultsPanelA and put them into the String
        // IMPORTANT: In that case we sort by position and not by position_rsID. So, maybe we
        // are going to lose some SNPs...
        // for(int chromo=chrStart; chromo <=chrEnd; chromo++) {

        Hashtable<String, Integer> resultsHashTableIndex = new Hashtable<>();
        // String chromoS = Integer.toString(chromo);

        // Create the first treeMap for the chromo
        TreeMap<String, String> fileTreeMapA = new TreeMap<>();

        String header = null;
        String positionA1A2Chr = null;
        try (Scanner sc1 = new Scanner(new File(resultsPanelA + ".temp"))) {
            // Get the header
            header = sc1.nextLine();

            if (!header.equals("chr\tposition\trs_id_all\tinfo_all\tcertainty_all\t")) {
                resultsHashTableIndex = createHashWithHeader(header, "\t");
                chrIdx = resultsHashTableIndex.get("chr");
                posIdx = resultsHashTableIndex.get("position");
                a1Idx = resultsHashTableIndex.get("alleleA");
                a2Idx = resultsHashTableIndex.get("alleleB");
                infoIdx = resultsHashTableIndex.get("info_all");

                while (sc1.hasNextLine()) {
                    String line = sc1.nextLine();
                    String[] splitted = line.split("\t");

                    // if(splitted[chrIdx].equals(chromoS)) {
                    positionA1A2Chr = splitted[posIdx] + "_" + splitted[a1Idx] + "_" + splitted[a2Idx] + "_" + splitted[chrIdx];
                    // Now, we put this String into the treemap with the key positionA1A1Chr
                    fileTreeMapA.put(positionA1A2Chr, line);
                    // contador++;
                    // }
                }
                // note that Scanner suppresses exceptions
                if (sc1.ioException() != null) {
                    throw sc1.ioException();
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // System.out.println("i\n[DEBUG] We have read the chromo " + chromoS + " from first File. contador = " +
        // contador);

        // Create the second treeMap for the chromo
        TreeMap<String, String> fileTreeMapB = new TreeMap<>();
        // contador=0;
        try (Scanner sc = new Scanner(new File(resultsPanelB + ".temp"))) {
            // Get the header
            header = sc.nextLine();

            if (!header.equals("chr\tposition\trs_id_all\tinfo_all\tcertainty_all\t")) {
                chrIdx = resultsHashTableIndex.get("chr");
                posIdx = resultsHashTableIndex.get("position");
                a1Idx = resultsHashTableIndex.get("alleleA");
                a2Idx = resultsHashTableIndex.get("alleleB");
                infoIdx = resultsHashTableIndex.get("info_all");

                while (sc.hasNextLine()) {
                    String line = sc.nextLine();
                    String[] splitted = line.split("\t");

                    // if(splitted[chrIdx].equals(chromoS)) {
                    positionA1A2Chr = splitted[posIdx] + "_" + splitted[a1Idx] + "_" + splitted[a2Idx] + "_" + splitted[chrIdx];
                    // Now, we put this String into the treemap with the key positionA1A1Chr
                    fileTreeMapB.put(positionA1A2Chr, line);
                    // contador++;
                    // }
                }
                // note that Scanner suppresses exceptions
                if (sc.ioException() != null) {
                    throw sc.ioException();
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // System.out.println("\n[DEBUG] We have read the chromo " + chromoS + " from second File. contador = " +
        // contador);

        // A place to store the results of this combining
        String lineA = null;
        String lineB = null;

        String[] splittedA = null;
        String[] splittedB = null;

        Double infoA;
        Double infoB;

        String posAllelesEqual = null;
        String posAllelesReverse = null;
        String posAllelesComplement = null;
        String posAllelesComplementAndReverse = null;

        TreeMap<String, String> fileTreeMapC = new TreeMap<>();
        // contador=0;

        // We first iterate the fileTreeMapA
        Set<Entry<String, String>> mySet = fileTreeMapA.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<String, String>> iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, String> m = iter.next();
            positionA1A2Chr = (String) m.getKey();
            lineA = (String) m.getValue();
            splittedA = lineA.split("\t");
            infoA = Double.parseDouble(splittedA[infoIdx]);

            posAllelesEqual = positionA1A2Chr;
            posAllelesReverse = splittedA[posIdx] + "_" + getAllele(splittedA[a1Idx], splittedA[a2Idx], "reverse") + "_"
                    + splittedA[chrIdx];
            posAllelesComplement = splittedA[posIdx] + "_" + getAllele(splittedA[a1Idx], splittedA[a2Idx], "complement") + "_"
                    + splittedA[chrIdx];
            posAllelesComplementAndReverse = splittedA[posIdx] + "_" + getAllele(splittedA[a1Idx], splittedA[a2Idx], "complementAndReverse")
                    + "_" + splittedA[chrIdx];

            // System.out.println("[combinePanelsComplex] " + positionA1A2Chr + " " + posAllelesEqual + " " +
            // posAllelesReverse + " " + posAllelesComplement + " " + posAllelesComplementAndReverse);

            // The same: position, a1 and a2?
            if (fileTreeMapB.containsKey(positionA1A2Chr)) {
                // If the fileTreeMapB contains this positionA1A2Chr combination, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                lineB = fileTreeMapB.get(positionA1A2Chr);
                splittedB = lineB.split("\t");
                infoB = Double.parseDouble(splittedB[infoIdx]);
                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(positionA1A2Chr, lineB);
                }
                // System.out.println("WOW alelos iguales: " + positionA1A2Chr);

                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(positionA1A2Chr);
            } else if (fileTreeMapB.containsKey(posAllelesReverse)) {
                // If the fileTreeMapB contains this posAllelesReverse, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                lineB = fileTreeMapB.get(posAllelesReverse);
                splittedB = lineB.split("\t");
                infoB = Double.parseDouble(splittedB[infoIdx]);
                // Then we have to choose between A and B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(posAllelesReverse, lineB);
                }
                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(posAllelesReverse);
                // System.out.println("WOW alelos reversos: " + positionA1A2Chr + " " + posAllelesReverse);
            } else if (fileTreeMapB.containsKey(posAllelesComplement)) {
                // If the fileTreeMapB contains this posAllelesComplement, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                lineB = fileTreeMapB.get(posAllelesComplement);
                splittedB = lineB.split("\t");
                infoB = Double.parseDouble(splittedB[infoIdx]);
                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(posAllelesComplement, lineB);
                }
                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(posAllelesComplement);
                // System.out.println("WOW alelos complementarios: " + positionA1A2Chr + " " +
                // posAllelesComplement);
            } else if (fileTreeMapB.containsKey(posAllelesComplementAndReverse)) {
                // If the fileTreeMapB contains this posAllelesComplement, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                lineB = fileTreeMapB.get(posAllelesComplementAndReverse);
                splittedB = lineB.split("\t");
                infoB = Double.parseDouble(splittedB[infoIdx]);
                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(posAllelesComplementAndReverse, lineB);
                }
                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(posAllelesComplementAndReverse);
                // System.out.println("WOW alelos complementariosYreversos: " + positionA1A2Chr + " " +
                // posAllelesComplementAndReverse);
            } else {
                // Else means that fileTreeMapB does not contain this SNP or any of its variants.
                // Therefore, we keep the one in fileTreeMapA
                fileTreeMapC.put(positionA1A2Chr, lineA);
                // System.out.println("WOW fileTreeMapB does not contain this SNP: " + positionA1A2Chr);
            }
            // contador++;
        }

        fileTreeMapA.clear();

        // Now we have to put in fileTreeMapC the rest of values that remain in fileTreeMapB.
        // We iterate the fileTreeMapB (the rest of the...)
        mySet = fileTreeMapB.entrySet();
        // Move next key and value of Map by iterator
        iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, String> m = iter.next();
            // getKey is used to get key of Map
            positionA1A2Chr = (String) m.getKey();
            lineB = (String) m.getValue();
            // Then we have to store the value in fileTreeMapC
            fileTreeMapC.put(positionA1A2Chr, lineB);
            // contador++;
        }

        fileTreeMapB.clear();
        // System.out.println("\n[DEBUG] We have processed the chromosome " + chromoS + ". contador " + contador);

        // Finally we put the fileTreeMapC into the outputFile
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            // We print the header which is the same always!.
            // if(chromo == chrStart) {
            writer.write(header);
            writer.newLine();
            // }

            String myLine = null;

            mySet = fileTreeMapC.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, String> m = iter.next();
                // getKey is used to get key of Map
                myLine = (String) m.getValue();

                writer.write(myLine);
                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        fileTreeMapC.clear();
        // System.out.println("\n[DEBUG] We have stored snps from chromosome " + chromoS + " in the output file");

        // Then, we create the gz file and rename it
        gzipFile(resultsPanelC, resultsPanelC + ".gz");
        File fc = new File(resultsPanelC);
        File fGz = new File(resultsPanelC + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        File fA = new File(resultsPanelA + ".temp");
        fA.delete();
        File fB = new File(resultsPanelB + ".temp");
        fB.delete();

        System.out.println("\n[DEBUG] Finished all chromosomes");

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] combinePanelsComplex startTime: " + startTime);
            System.out.println("\n[DEBUG] combinePanelsComplex endTime: " + stopTime);
            System.out.println("\n[DEBUG] combinePanelsComplex elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of combinePanelsComplex");
        }

    }

    /**
     * Method to combine condensed files
     * 
     * @param filteredA
     * @param filteredX
     * @param combinedCondensedFile
     * @param mafThresholdS
     * @param infoThresholdS
     * @param hweCohortThresholdS
     * @param hweCasesThresholdS
     * @param hweControlsThresholdS
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void combineCondensedFiles(String filteredA, String filteredX, String combinedCondensedFile, String mafThresholdS,
            String infoThresholdS, String hweCohortThresholdS, String hweCasesThresholdS, String hweControlsThresholdS, String cmdToStore)
            throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running combineCondensedFiles with parameters:");
            System.out.println("[DEBUG] \t- Input filteredA             : " + filteredA);
            System.out.println("[DEBUG] \t- Input filteredX             : " + filteredX);
            System.out.println("[DEBUG] \t- Output combinedCondensedFile : " + combinedCondensedFile);
            System.out.println("[DEBUG] \t- Input maf threshold           : " + mafThresholdS);
            System.out.println("[DEBUG] \t- Input info threshold          : " + infoThresholdS);
            System.out.println("[DEBUG] \t- Input hwe cohort threshold    : " + hweCohortThresholdS);
            System.out.println("[DEBUG] \t- Input hwe controls threshold  : " + hweCasesThresholdS);
            System.out.println("[DEBUG] \t- Input hwe cases threshold     : " + hweControlsThresholdS);
            System.out.println("\n");

            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // Convert threshold string into thresholdDouble
        Double mafThreshold = Double.parseDouble(mafThresholdS);
        Double infoThreshold = Double.parseDouble(infoThresholdS);
        Double hweCohortThreshold = Double.parseDouble(hweCohortThresholdS);
        Double hweCasesThreshold = Double.parseDouble(hweCasesThresholdS);
        Double hweControlsThreshold = Double.parseDouble(hweControlsThresholdS);

        File outCondensedFile = new File(combinedCondensedFile);

        // Tries to create the file
        // Print information about the existence of the file
        // System.out.println("\n[DEBUG] \t- Output file " + combinedCondensedFile + " was succesfuly created? " +
        // bool);

        Hashtable<String, Integer> inputFileHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> inputFileHashTableIndexReversed = new Hashtable<>();

        try (BufferedWriter writerCondensed = new BufferedWriter(new FileWriter(outCondensedFile))) {

            try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(filteredA));
                    Reader decoder = new InputStreamReader(inputGz);
                    BufferedReader br = new BufferedReader(decoder)) {
                // I read the header
                String line = br.readLine();

                inputFileHashTableIndex = createHashWithHeader(line, "\t");
                inputFileHashTableIndexReversed = createHashWithHeaderReversed(line, "\t");

                String headerCondensed = "CHR\tBP\tP";
                writerCondensed.write(headerCondensed);
                writerCondensed.newLine();

                while ((line = br.readLine()) != null) {
                    String[] splittedLine = line.split("\t");// delimiter I assume single space.

                    String chromo = splittedLine[inputFileHashTableIndex.get("chr")];

                    String infoS = splittedLine[inputFileHashTableIndex.get("info_all")];

                    // We start with these values for hwe values just to allows the X chromosome to pass the if
                    // statement of
                    // the
                    // next lines
                    // Just remember that hwe filtering when chromo X is being processed does not make sense.
                    String hwe_cohortS = "1.0";
                    String hwe_casesS = "1.0";
                    String hwe_controlsS = "1.0";

                    if (!chromo.equals("23")) {
                        hwe_cohortS = splittedLine[inputFileHashTableIndex.get("cohort_1_hwe")];
                        hwe_casesS = splittedLine[inputFileHashTableIndex.get("cases_hwe")];
                        hwe_controlsS = splittedLine[inputFileHashTableIndex.get("controls_hwe")];
                    }

                    String cases_mafS = splittedLine[inputFileHashTableIndex.get("cases_maf")];
                    String controls_mafS = splittedLine[inputFileHashTableIndex.get("controls_maf")];

                    String position = splittedLine[inputFileHashTableIndex.get("position")];
                    String pva = splittedLine[inputFileHashTableIndex.get("frequentist_add_pvalue")];

                    String chrbpb = chromo + "\t" + position + "\t" + pva;

                    if (!cases_mafS.equals("NA") && !controls_mafS.equals("NA") && !infoS.equals("NA") && !hwe_cohortS.equals("NA")
                            && !hwe_casesS.equals("NA") && !hwe_controlsS.equals("NA") && !pva.equals("NA")) {
                        Double cases_maf = Double.parseDouble(cases_mafS);
                        Double controls_maf = Double.parseDouble(controls_mafS);
                        Double info = Double.parseDouble(infoS);
                        Double hweCohort = 1.0;
                        Double hweCases = 1.0;
                        Double hweControls = 1.0;

                        if (!chromo.equals("23")) {
                            hweCohort = Double.parseDouble(hwe_cohortS);
                            hweCases = Double.parseDouble(hwe_casesS);
                            hweControls = Double.parseDouble(hwe_controlsS);
                        }

                        if (cases_maf >= mafThreshold && controls_maf >= mafThreshold && info >= infoThreshold && // VERIFICAR
                                                                                                                  // LA
                                                                                                                  // CONDICION
                                hweCohort >= hweCohortThreshold && hweCases >= hweCasesThreshold && hweControls >= hweControlsThreshold) {

                            writerCondensed.write(chrbpb);
                            writerCondensed.newLine();
                        }
                    }
                }
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }

            // Now with crh 23
            // If filteredA != filteredX then there is chr23 file (filteredX, therefore we have to include it in the
            // results.
            // Other wise, there is nothing to do.
            if (!filteredA.equals(filteredX)) {
                try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(filteredX));
                        Reader decoder = new InputStreamReader(inputGz);
                        BufferedReader br = new BufferedReader(decoder);) {

                    inputFileHashTableIndex = new Hashtable<>();
                    inputFileHashTableIndexReversed = new Hashtable<>();

                    // I read the header
                    String line = br.readLine();

                    inputFileHashTableIndex = createHashWithHeader(line, "\t");
                    inputFileHashTableIndexReversed = createHashWithHeaderReversed(line, "\t");

                    while ((line = br.readLine()) != null) {
                        String[] splittedLine = line.split("\t");// delimiter I assume single space.

                        String chromo = splittedLine[inputFileHashTableIndex.get("chr")];

                        String infoS = splittedLine[inputFileHashTableIndex.get("info_all")];

                        // We start with these values for hwe values just to allows the X chromosome to pass the if
                        // statement of
                        // the next lines
                        // Just remember that hwe filtering when chromo X is being processed does not make sense.
                        String hwe_cohortS = "1.0";
                        String hwe_casesS = "1.0";
                        String hwe_controlsS = "1.0";

                        if (!chromo.equals("23")) {
                            hwe_cohortS = splittedLine[inputFileHashTableIndex.get("cohort_1_hwe")];
                            hwe_casesS = splittedLine[inputFileHashTableIndex.get("cases_hwe")];
                            hwe_controlsS = splittedLine[inputFileHashTableIndex.get("controls_hwe")];
                        }

                        String cases_mafS = splittedLine[inputFileHashTableIndex.get("cases_maf")];
                        String controls_mafS = splittedLine[inputFileHashTableIndex.get("controls_maf")];

                        String position = splittedLine[inputFileHashTableIndex.get("position")];
                        String pva = splittedLine[inputFileHashTableIndex.get("frequentist_add_pvalue")];

                        String chrbpb = chromo + "\t" + position + "\t" + pva;

                        if (!cases_mafS.equals("NA") && !controls_mafS.equals("NA") && !infoS.equals("NA") && !hwe_cohortS.equals("NA")
                                && !hwe_casesS.equals("NA") && !hwe_controlsS.equals("NA") && !pva.equals("NA")) {
                            Double cases_maf = Double.parseDouble(cases_mafS);
                            Double controls_maf = Double.parseDouble(controls_mafS);
                            Double info = Double.parseDouble(infoS);
                            Double hweCohort = 1.0;
                            Double hweCases = 1.0;
                            Double hweControls = 1.0;

                            if (!chromo.equals("23")) {
                                hweCohort = Double.parseDouble(hwe_cohortS);
                                hweCases = Double.parseDouble(hwe_casesS);
                                hweControls = Double.parseDouble(hwe_controlsS);
                            }

                            if (cases_maf >= mafThreshold && controls_maf >= mafThreshold && info >= infoThreshold && // VERIFICAR
                                                                                                                      // LA
                                                                                                                      // CONDICION
                                    hweCohort >= hweCohortThreshold && hweCases >= hweCasesThreshold
                                    && hweControls >= hweControlsThreshold) {

                                writerCondensed.write(chrbpb);
                                writerCondensed.newLine();
                            }
                        }
                    }
                } catch (IOException ioe) {
                    throw new GuidanceTaskException(ioe);
                }
            }

            writerCondensed.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        gzipFile(combinedCondensedFile, combinedCondensedFile + ".gz");
        File fc = new File(combinedCondensedFile);
        File fGz = new File(combinedCondensedFile + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] combineCondensedFiles startTime: " + startTime);
            System.out.println("\n[DEBUG] combineCondensedFiles endTime: " + stopTime);
            System.out.println("\n[DEBUG] combineCondensedFiles elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of combinedCondensedFiles");
        }
    }

    /**
     * Method to combine panels complex (version 2)
     * 
     * @param resultsPanelA
     * @param resultsPanelB
     * @param resultsPanelC
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void combinePanelsComplex1(String resultsPanelA, String resultsPanelB, String resultsPanelC, String cmdToStore)
            throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running combinePanelsComplex1 with parameters:");
            System.out.println("[DEBUG] \t- resultsPanelA             : " + resultsPanelA);
            System.out.println("[DEBUG] \t- resultsPanelB             : " + resultsPanelB);
            System.out.println("[DEBUG] \t- resultsPanelC             : " + resultsPanelC);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // We read each line of the resultsPanelA and put them into the String
        // IMPORTANT: In that case we sort by position and not by position_rsID. So, maybe we
        // are going to lose some SNPs...
        TreeMap<String, String> fileTreeMapA = new TreeMap<>();
        String header = null;
        try (FileReader frA = new FileReader(resultsPanelA); BufferedReader brA = new BufferedReader(frA)) {
            String line = null;
            // First: read the header and avoid it
            header = brA.readLine();

            while ((line = brA.readLine()) != null) {
                String[] splitted = line.split("\t");
                String positionA1A2Chr = splitted[0] + "_" + splitted[4] + "_" + splitted[5] + "_" + splitted[35];
                // Now, we put this String into the treemap with the key positionA1A1Chr
                fileTreeMapA.put(positionA1A2Chr, line);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // We read each line of the resultsPanelB and put them into fileAList array of Strings
        TreeMap<String, String> fileTreeMapB = new TreeMap<>();
        try (FileReader frB = new FileReader(resultsPanelB); BufferedReader brB = new BufferedReader(frB)) {
            // First: read the header and avoid it
            String line = brB.readLine();

            while ((line = brB.readLine()) != null) {
                String[] splitted = line.split("\t");
                String positionA1A2Chr = splitted[0] + "_" + splitted[4] + "_" + splitted[5] + "_" + splitted[35];
                // Now, we put this String into the treemap with the key positionA1A2Chr
                fileTreeMapB.put(positionA1A2Chr, line);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // A place to store the results of this combining
        TreeMap<String, String> fileTreeMapC = new TreeMap<>();
        // We first iterate the fileTreeMapA
        Set<Entry<String, String>> mySet = fileTreeMapA.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<String, String>> iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, String> m = iter.next();

            String positionA1A2Chr = (String) m.getKey();
            String lineA = (String) m.getValue();
            String[] splittedA = lineA.split("\t");
            Double infoA = Double.parseDouble(splittedA[2]);

            String posAllelesEqual = positionA1A2Chr;
            String posAllelesReverse = splittedA[0] + "_" + getAllele(splittedA[4], splittedA[5], "reverse") + "_" + splittedA[35];
            String posAllelesComplement = splittedA[0] + "_" + getAllele(splittedA[4], splittedA[5], "complement") + "_" + splittedA[35];
            String posAllelesComplementAndReverse = splittedA[0] + "_" + getAllele(splittedA[4], splittedA[5], "complementAndReverse") + "_"
                    + splittedA[35];

            // System.out.println("[combinePanelsComplex] " + positionA1A2Chr + " " + posAllelesEqual + " " +
            // posAllelesReverse + " " + posAllelesComplement + " " + posAllelesComplementAndReverse);

            // The same: position, a1 and a2?
            if (fileTreeMapB.containsKey(positionA1A2Chr)) {
                // If the fileTreeMapB contains this positionA1A2Chr combination, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                String lineB = fileTreeMapB.get(positionA1A2Chr);
                String[] splittedB = lineB.split("\t");
                Double infoB = Double.parseDouble(splittedB[2]);
                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(positionA1A2Chr, lineB);
                }
                // System.out.println("WOW alelos iguales: " + positionA1A2Chr);

                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(positionA1A2Chr);
            } else if (fileTreeMapB.containsKey(posAllelesReverse)) {
                // If the fileTreeMapB contains this posAllelesReverse, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                String lineB = fileTreeMapB.get(posAllelesReverse);
                String[] splittedB = lineB.split("\t");
                Double infoB = Double.parseDouble(splittedB[2]);
                // Then we have to choose between A and B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(posAllelesReverse, lineB);
                }
                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(posAllelesReverse);
                // System.out.println("WOW alelos reversos: " + positionA1A2Chr + " " + posAllelesReverse);
            } else if (fileTreeMapB.containsKey(posAllelesComplement)) {
                // If the fileTreeMapB contains this posAllelesComplement, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                String lineB = fileTreeMapB.get(posAllelesComplement);
                String[] splittedB = lineB.split("\t");
                Double infoB = Double.parseDouble(splittedB[2]);
                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(posAllelesComplement, lineB);
                }
                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(posAllelesComplement);
                // System.out.println("WOW alelos complementarios: " + positionA1A2Chr + " " + posAllelesComplement);
            } else if (fileTreeMapB.containsKey(posAllelesComplementAndReverse)) {
                // If the fileTreeMapB contains this posAllelesComplement, then we have to choose
                // the ones that has a better info (that is the ones with greater info).
                String lineB = fileTreeMapB.get(posAllelesComplementAndReverse);
                String[] splittedB = lineB.split("\t");
                Double infoB = Double.parseDouble(splittedB[2]);
                // Then we have to choose between A o B.
                if (infoA >= infoB) {
                    fileTreeMapC.put(positionA1A2Chr, lineA);
                } else {
                    fileTreeMapC.put(posAllelesComplementAndReverse, lineB);
                }
                // Now we remove this value from the fileTreeMapB
                fileTreeMapB.remove(posAllelesComplementAndReverse);
                // System.out.println("WOW alelos complementariosYreversos: " + positionA1A2Chr + " " +
                // posAllelesComplementAndReverse);
            } else {
                // Else means that fileTreeMapB does not contain this SNP or any of its variants.
                // Therefore, we keep the one in fileTreeMapA
                fileTreeMapC.put(positionA1A2Chr, lineA);
            }
        }

        // Now we have to put in fileTreeMapC the rest of values that remain in fileTreeMapB.
        // We iterate the fileTreeMapB (the rest of the...)
        mySet = fileTreeMapB.entrySet();
        // Move next key and value of Map by iterator
        iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, String> m = iter.next();
            // getKey is used to get key of Map
            String positionA1A2Chr = (String) m.getKey();
            String lineB = (String) m.getValue();
            // Then we have to store the value in fileTreeMapC
            fileTreeMapC.put(positionA1A2Chr, lineB);
        }

        // Finally we put the fileTreeMapC into the outputFile
        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File outputFile = new File(resultsPanelC);
        try {
            outputFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            // We print the header which is the same always!.
            writer.write(header);
            writer.newLine();

            mySet = fileTreeMapC.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, String> m = iter.next();
                // getKey is used to get key of Map
                String myLine = (String) m.getValue();

                writer.write(myLine);
                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] combinePanelsComplex1 startTime: " + startTime);
            System.out.println("\n[DEBUG] combinePanelsComplex1 endTime: " + stopTime);
            System.out.println("\n[DEBUG] combinePanelsComplex1 elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of combinePanelsComplex1");
        }

    }

    /**
     * Method to generate top hits
     * 
     * @param resultsFile
     * @param outputTopHitFile
     * @param pvaThreshold
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void generateTopHits(String resultsFile, String outputTopHitFile, String pvaThreshold, String cmdToStore)
            throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running generateTopHits with parameters:");
            System.out.println("[DEBUG] \t- resultsFile                : " + resultsFile);
            System.out.println("[DEBUG] \t- outputTopHitFile           : " + outputTopHitFile);
            System.out.println("[DEBUG] \t- pvaThreshold               : " + pvaThreshold);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();
        double pvaThres = Double.parseDouble(pvaThreshold);

        // double pvaThres = Double.parseDouble(splitted[17]);

        // We read each line of the resultsFile and put them into the String
        TreeMap<String, String> fileTreeMap = new TreeMap<>();

        String header = null;
        String newHeader = null;
        Hashtable<String, Integer> resultsFileHashTableIndex = new Hashtable<>();
        // First: read the header and avoid it
        try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(resultsFile));
                Reader decoder = new InputStreamReader(inputGz);
                BufferedReader br = new BufferedReader(decoder)) {
            header = br.readLine();
            resultsFileHashTableIndex = createHashWithHeader(header, "\t");

            int indexPosition = resultsFileHashTableIndex.get("position");
            int indexRsId = resultsFileHashTableIndex.get("rs_id_all");

            System.out.println("ANTES 3");

            int indexPvalue = resultsFileHashTableIndex.get("frequentist_add_pvalue");
            int indexChromo = resultsFileHashTableIndex.get("chr");
            int indexAllMaf = resultsFileHashTableIndex.get("all_maf");
            int indexAlleleA = resultsFileHashTableIndex.get("alleleA");
            int indexAlleleB = resultsFileHashTableIndex.get("alleleB");

            newHeader = "chr\tposition\trsid\tMAF\ta1\ta2\tpval_add";
            String line = null;
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split("\t");
                String positionAndRsId = splitted[indexPosition] + "_" + splitted[indexRsId];
                double myPva = Double.parseDouble(splitted[indexPvalue]);

                if (myPva <= pvaThres && myPva > 0.0) {
                    // Now, we put this String into the treemap with the key positionAndRsId
                    // reducedLine is chr;position;RSID_ALL;MAF;a1;a2;pval
                    String reducedLine = splitted[indexChromo] + "\t" + splitted[indexPosition] + "\t" + splitted[indexRsId] + "\t"
                            + splitted[indexAllMaf] + "\t" + splitted[indexAlleleA] + "\t" + splitted[indexAlleleB] + "\t"
                            + splitted[indexPvalue];
                    fileTreeMap.put(positionAndRsId, reducedLine);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Finally we put the fileTreeMap into the output file
        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File outputFile = new File(outputTopHitFile);
        try {
            outputFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            // We print the header which is the same always!.
            writer.write(newHeader);
            // writer.newLine();

            Set<Entry<String, String>> mySet = fileTreeMap.entrySet();
            // Move next key and value of Map by iterator
            Iterator<Entry<String, String>> iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, String> m = iter.next();
                // getKey is used to get key of Map
                String myLine = (String) m.getValue();

                writer.newLine();
                writer.write(myLine);
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(outputTopHitFile, outputTopHitFile + ".gz");
        File fc = new File(outputTopHitFile);
        File fGz = new File(outputTopHitFile + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] generateTopHits startTime:  " + startTime);
            System.out.println("\n[DEBUG] generateTopHits endTime:    " + stopTime);
            System.out.println("\n[DEBUG] generateTopHits elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of generateTopHits");
        }
    }

    /**
     * Method to generate all the top hits
     * 
     * @param resultsAFile
     * @param resultsBFile
     * @param outputTopHitFile
     * @param pvaThreshold
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void generateTopHitsAll(String resultsAFile, String resultsBFile, String outputTopHitFile, String pvaThreshold,
            String cmdToStore) throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running generateTopHits with parameters:");
            System.out.println("[DEBUG] \t- resultsAFile                : " + resultsAFile);
            System.out.println("[DEBUG] \t- resultsBFile                : " + resultsBFile);
            System.out.println("[DEBUG] \t- outputTopHitFile           : " + outputTopHitFile);
            System.out.println("[DEBUG] \t- pvaThreshold               : " + pvaThreshold);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();
        double pvaThres = Double.parseDouble(pvaThreshold);

        // double pvaThres = Double.parseDouble(splitted[17]);

        // We read each line of the resultsAFile and put them into the String
        TreeMap<String, String> fileATreeMap = new TreeMap<>();
        TreeMap<String, String> fileBTreeMap = new TreeMap<>();

        String header = null;
        String newHeader = null;
        Hashtable<String, Integer> resultsAFileHashTableIndex = new Hashtable<>();
        Hashtable<String, Integer> resultsBFileHashTableIndex = new Hashtable<>();

        // First: read the header and avoid it
        try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(resultsAFile));
                Reader decoder = new InputStreamReader(inputGz);
                BufferedReader br = new BufferedReader(decoder)) {

            header = br.readLine();
            resultsAFileHashTableIndex = createHashWithHeader(header, "\t");

            int indexPosition = resultsAFileHashTableIndex.get("position");
            int indexRsId = resultsAFileHashTableIndex.get("rs_id_all");

            int indexPvalue = resultsAFileHashTableIndex.get("frequentist_add_pvalue");
            int indexChromo = resultsAFileHashTableIndex.get("chr");
            int indexAllMaf = resultsAFileHashTableIndex.get("all_maf");
            int indexAlleleA = resultsAFileHashTableIndex.get("alleleA");
            int indexAlleleB = resultsAFileHashTableIndex.get("alleleB");

            newHeader = "chr\tposition\trsid\tMAF\ta1\ta2\tpval_add";
            String line = null;
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split("\t");
                String positionAndRsId = splitted[indexPosition] + "_" + splitted[indexRsId];
                double myPva = Double.parseDouble(splitted[indexPvalue]);

                if (myPva <= pvaThres && myPva > 0.0) {
                    // Now, we put this String into the treemap with the key positionAndRsId
                    // reducedLine is chr;position;RSID_ALL;MAF;a1;a2;pval
                    String reducedLine = splitted[indexChromo] + "\t" + splitted[indexPosition] + "\t" + splitted[indexRsId] + "\t"
                            + splitted[indexAllMaf] + "\t" + splitted[indexAlleleA] + "\t" + splitted[indexAlleleB] + "\t"
                            + splitted[indexPvalue];
                    fileATreeMap.put(positionAndRsId, reducedLine);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Now we have to see if we have to include results for Chr23 that come in the second input file (resultsBFile)
        // The way to know whether we have results for chr23 is by checking that resultsAFile is equal to resultsBFile.
        // If they are equal, the we do not have chr23. Other wise we have to results for chr23 y resultsBFile and we
        // have
        // to include it in the outputTopHitFile
        if (!resultsAFile.equals(resultsBFile)) {
            try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(resultsBFile));
                    Reader decoder = new InputStreamReader(inputGz);
                    BufferedReader br = new BufferedReader(decoder)) {
                resultsBFileHashTableIndex = new Hashtable<>();

                // First: read the header and avoid it
                header = br.readLine();
                resultsBFileHashTableIndex = createHashWithHeader(header, "\t");

                int indexPosition = resultsBFileHashTableIndex.get("position");
                int indexRsId = resultsBFileHashTableIndex.get("rs_id_all");

                System.out.println("ANTES X");

                int indexPvalue = resultsBFileHashTableIndex.get("frequentist_add_pvalue");
                int indexChromo = resultsBFileHashTableIndex.get("chr");
                int indexAllMaf = resultsBFileHashTableIndex.get("all_maf");
                int indexAlleleA = resultsBFileHashTableIndex.get("alleleA");
                int indexAlleleB = resultsBFileHashTableIndex.get("alleleB");

                String line = null;
                while ((line = br.readLine()) != null) {
                    String[] splitted = line.split("\t");
                    String positionAndRsId = splitted[indexPosition] + "_" + splitted[indexRsId];
                    double myPva = Double.parseDouble(splitted[indexPvalue]);

                    if (myPva <= pvaThres && myPva > 0.0) {
                        // Now, we put this String into the treemap with the key positionAndRsId
                        // reducedLine is chr;position;RSID_ALL;MAF;a1;a2;pval
                        String reducedLine = splitted[indexChromo] + "\t" + splitted[indexPosition] + "\t" + splitted[indexRsId] + "\t"
                                + splitted[indexAllMaf] + "\t" + splitted[indexAlleleA] + "\t" + splitted[indexAlleleB] + "\t"
                                + splitted[indexPvalue];
                        fileBTreeMap.put(positionAndRsId, reducedLine);
                    }
                }
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        // Then we put the fileATreeMap into the output file
        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File outputFile = new File(outputTopHitFile);
        try {
            outputFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
            // We print the header which is the same always!.
            writer.write(newHeader);
            // writer.newLine();

            Set<Entry<String, String>> mySet = fileATreeMap.entrySet();
            // Move next key and value of Map by iterator
            Iterator<Entry<String, String>> iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, String> m = iter.next();
                // getKey is used to get key of Map
                String myLine = (String) m.getValue();

                writer.newLine();
                writer.write(myLine);
            }

            // Finally we put the fileBTreeMap into the output file.
            // This should happens only if we really have a resultsAFile != resultsBFile
            if (!resultsAFile.equals(resultsBFile)) {
                mySet = fileBTreeMap.entrySet();
                // Move next key and value of Map by iterator
                iter = mySet.iterator();
                while (iter.hasNext()) {
                    // key=value separator this by Map.Entry to get key and value
                    Entry<String, String> m = iter.next();
                    // getKey is used to get key of Map
                    String myLine = (String) m.getValue();

                    writer.newLine();
                    writer.write(myLine);
                }
            }

            // Now we close the output buffer.
            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(outputTopHitFile, outputTopHitFile + ".gz");
        File fc = new File(outputTopHitFile);
        File fGz = new File(outputTopHitFile + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] generateTopHits startTime:  " + startTime);
            System.out.println("\n[DEBUG] generateTopHits endTime:    " + stopTime);
            System.out.println("\n[DEBUG] generateTopHits elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of generateTopHits");
        }
    }

    /**
     * Method to generate QQ Manhattan Plots
     * 
     * @param lastCondensedFile
     * @param qqPlotFile
     * @param manhattanPlotFile
     * @param qqPlotTiffFile
     * @param manhattanPlotTiffFile
     * @param correctedPvaluesFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void generateQQManhattanPlots(String lastCondensedFile, String qqPlotFile, String manhattanPlotFile,
            String qqPlotTiffFile, String manhattanPlotTiffFile, String correctedPvaluesFile, String cmdToStore)
            throws GuidanceTaskException {

        String rScriptBinDir = System.getenv("RSCRIPTBINDIR");
        if (rScriptBinDir == null) {
            throw new GuidanceTaskException(
                    "[generateQQManhattanPlots] Error, RSCRIPTBINDIR environment variable is not defined in .bashrc!!!");
        }

        String rScriptDir = System.getenv("RSCRIPTDIR");
        if (rScriptDir == null) {
            throw new GuidanceTaskException(
                    "[generateQQManhattanPlots] Error, RSCRIPTDIR environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running generateQQManhattanPlots with parameters:");
            System.out.println("[DEBUG] \t- lastCondensedFile             : " + lastCondensedFile);
            System.out.println("[DEBUG] \t- qqPlotFile                    : " + qqPlotFile);
            System.out.println("[DEBUG] \t- manhattanPlotFile             : " + manhattanPlotFile);
            System.out.println("[DEBUG] \t- qqPlotTiffFile                : " + qqPlotTiffFile);
            System.out.println("[DEBUG] \t- manhattanPlotTiffFile         : " + manhattanPlotTiffFile);
            System.out.println("[DEBUG] \t- Output outputCondensedFile    : " + correctedPvaluesFile);

            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND            : " + cmdToStore);

        }

        long startTime = System.currentTimeMillis();

        // First, we have to decompress the input file
        String theInputFile = lastCondensedFile + ".temp1";
        gunzipFile(lastCondensedFile, theInputFile);

        String cmd = null;
        cmd = rScriptBinDir + "/Rscript " + rScriptDir + "/qqplot_manhattan.R " + theInputFile + " " + qqPlotFile + " " + manhattanPlotFile
                + " " + qqPlotTiffFile + " " + manhattanPlotTiffFile + " " + correctedPvaluesFile;

        if (DEBUG) {
            System.out.println("\n[DEBUG] Cmd -> " + cmd);
            System.out.println(" ");
        }
        ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
        pb.environment().remove("LD_PRELOAD");
        Process generateQQManProc;
        try {
            generateQQManProc = pb.start();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Check the proper ending of the process
        int exitValue = -1;
        try {
            exitValue = generateQQManProc.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        if (exitValue != 0) {
            throw new GuidanceTaskException(
                    "[generateQQManhattanPlots] Error executing generateQQManProc job, exit value is: " + exitValue);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] generateQQManhattanPlots startTime: " + startTime);
            System.out.println("\n[DEBUG] generateQQManhattanPlots endTime: " + stopTime);
            System.out.println("\n[DEBUG] generateQQManhattanPlots elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of generate generateQQManhattanPlots.");
        }

    }

    /**
     * Method to perform the SNPTest
     * 
     * @param mergedGenFile
     * @param mergedSampleFile
     * @param snptestOutFile
     * @param snptestLogFile
     * @param responseVar
     * @param covariables
     * @param theChromo
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void snptest(String mergedGenFile, String mergedSampleFile, String snptestOutFile, String snptestLogFile,
            String responseVar, String covariables, String theChromo, String cmdToStore) throws GuidanceTaskException {

        String snptestBinary = System.getenv("SNPTESTBINARY");
        if (snptestBinary == null) {
            throw new GuidanceTaskException("[snptest] Error, SNPTESTBINARY environment variable is not defined in .bashrc!!!");
        }

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running snptest with parameters:");
            System.out.println("[DEBUG] \t- snptestBinary                    : " + snptestBinary);
            System.out.println("[DEBUG] \t- Input mergedGenFile              : " + mergedGenFile);
            System.out.println("[DEBUG] \t- Input mergedSampleFile           : " + mergedSampleFile);
            System.out.println("[DEBUG] \t- Output snptestOutFile            : " + snptestOutFile);
            System.out.println("[DEBUG] \t- Output snptestLogFile            : " + snptestLogFile);
            System.out.println("[DEBUG] \t- Input responseVar               : " + responseVar);
            System.out.println("[DEBUG] \t- Input covariables                : " + covariables);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND            : " + cmdToStore);
        }
        // replace commas in the string covariables
        String newStr = covariables.replace(',', ' ');
        if (DEBUG) {
            System.out.println("[DEBUG] \t- Changing covariable format. New covariables : " + newStr);
        }

        long startTime = System.currentTimeMillis();

        // We have to make sure whether we are using renamed files of the original gz files.
        // We detect this situation by scanning the last three characters:
        String extension = mergedGenFile.substring(Math.max(0, mergedGenFile.length() - 3));
        // System.out.println("DEBUG \t The file extension is: " + extension + " and the file is " + mergedGenFile);

        String mergedGenFileGz = null;
        if (extension.equals(".gz")) {
            mergedGenFileGz = mergedGenFile;
        } else {
            // If mergeGenFile exists, then mergedGenFileGz exists also.
            // We reused the mergedGenFileGz
            mergedGenFileGz = mergedGenFile + ".gz";
        }

        // Now we create the output file name: snptestOutFileGz
        String snptestOutFileGz = snptestOutFile + ".gz";

        String cmd = null;

        // Before executing snptest, I have to verify that the input mergedGenFile is not empty
        int nBytes = -1;
        try (FileInputStream fis = new FileInputStream(new File(mergedGenFile))) {
            nBytes = fis.read();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        if (nBytes != -1) {
            if (covariables.equals("none")) {
                cmd = snptestBinary + " -data " + mergedGenFileGz + " " + mergedSampleFile + " -o " + snptestOutFile + " -pheno "
                        + responseVar + " -hwe -log " + snptestLogFile;
            } else {
                cmd = snptestBinary + " -data " + mergedGenFileGz + " " + mergedSampleFile + " -o " + snptestOutFile + " -pheno "
                        + responseVar + " -cov_names " + newStr + " -hwe -log " + snptestLogFile;
            }

            // Different parameters for chromo 23 (X) and the rest.
            if (theChromo.equals("23")) {
                cmd = cmd + " -method newml -assume_chromosome X -stratify_on sex -frequentist 1 ";
            } else {
                cmd = cmd + " -method em -frequentist 1 2 3 4 5 ";
            }

            if (DEBUG) {
                System.out.println("\n[DEBUG] Cmd -> " + cmd);
                System.out.println(" ");
            }
            ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
            pb.environment().remove("LD_PRELOAD");
            Process snptestProc;
            try {
                snptestProc = pb.start();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }

            // Handling the streams so that dead lock situation never occurs
            try (BufferedInputStream bisInp = new BufferedInputStream(snptestProc.getInputStream());
                    BufferedOutputStream bosInp = new BufferedOutputStream(new FileOutputStream(snptestLogFile + ".stdout"))) {

                byte[] b = new byte[10240];
                int read;
                while ((read = bisInp.read(b)) >= 0) {
                    bosInp.write(b, 0, read);
                }
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            try (BufferedInputStream bisErr = new BufferedInputStream(snptestProc.getErrorStream());
                    BufferedOutputStream bosErr = new BufferedOutputStream(new FileOutputStream(snptestLogFile + ".stderr"))) {

                byte[] b = new byte[10240];
                int read;
                while ((read = bisErr.read(b)) >= 0) {
                    bosErr.write(b, 0, read);
                }
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }

            // Check the proper ending of the process
            int exitValue = -1;
            try {
                exitValue = snptestProc.waitFor();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
            if (exitValue != 0) {
                System.err.println("[snptest] Error executing snptestProc job, exit value is: " + exitValue);
                System.err.println("                         (This error is not fatal).");
                // throw new Exception("Error executing snptestProc job, exit value is: " + exitValue);
            }
        }

        // If there is not output in the snptest process. Then we have to create some empty outputs.
        File fa = new File(snptestOutFile);
        if (!fa.exists()) {
            try {
                fa.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        // Then, we create the gz file and rename it
        gzipFile(snptestOutFile, snptestOutFileGz);

        File fGz = new File(snptestOutFileGz);
        try {
            copyFile(fGz, fa);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        File fb = new File(snptestLogFile);
        if (!fb.exists()) {
            try {
                fb.createNewFile();
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] snptest startTime: " + startTime);
            System.out.println("\n[DEBUG] snptest endTime: " + stopTime);
            System.out.println("\n[DEBUG] snptest elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of snptest.");
        }
    }

    /*
     * public static void mergeTwoChunks(String reduceFileA, String reduceFileB, String reduceFileC) throws IOException,
     * InterruptedException, Exception { if (debug) { System.out.println("\nRunning mergeTwoChunks with parameters:");
     * System.out.println("\t- Input reduceFileA            : " + reduceFileA);
     * System.out.println("\t- Input reduceFileB            : " + reduceFileB);
     * System.out.println("\t- Output reduceFileC           : " + reduceFileC); }
     * 
     * String cmd = null; String workingDir = System.getenv("PWD");
     * 
     * // Before executing snptest, I have to verify that the input mergedGenFile is not empty. cmd = "perl " +
     * perlScriptDir + "/merge_two_chunks.pl " + reduceFileA + " " + reduceFileB + " " + reduceFileC;
     * 
     * if(debug){ System.out.println("\nCmd -> "+cmd); System.out.println(" "); } ProcessBuilder pb = new
     * ProcessBuilder(cmd); pb.environment().remove("LD_PRELOAD"); Process mergeTwoChunksProc = pb.start();
     * 
     * // Check the proper ending of the process int exitValue = mergeTwoChunksProc.waitFor(); if (exitValue != 0) {
     * System.err.println("Error executing mergeTwoChunksProc job, exit value is: " + exitValue); //throw new
     * Exception("Error executing snptestProc job, exit value is: " + exitValue); }
     * 
     * }
     */

    /**
     * Method to merge two chunks
     * 
     * @param reduceFileA
     * @param reduceFileB
     * @param reduceFileC
     * @param chrS
     * @param type
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void mergeTwoChunks(String reduceFileA, String reduceFileB, String reduceFileC, String chrS, String type,
            String cmdToStore) throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running mergeTwoChunks with parameters:");
            System.out.println("[DEBUG] \t- Input reduceFileA            : " + reduceFileA);
            System.out.println("[DEBUG] \t- Input reduceFileB            : " + reduceFileB);
            System.out.println("[DEBUG] \t- Output reduceFileC           : " + reduceFileC);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND         : " + cmdToStore);

            // long freeMemory = Runtime.getRuntime().freeMemory()/1048576;
            // long totalMemory = Runtime.getRuntime().totalMemory()/1048576;
            // long maxMemory = Runtime.getRuntime().maxMemory()/1048576;

            // System.out.println("JVM freeMemory: " + freeMemory);
            // System.out.println("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
            // System.out.println("JVM maxMemory also equals to maximum heap size of JVM : " + maxMemory);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
            // System.out.println("--------------------------------------");
        }
        long startTime = System.currentTimeMillis();

        // We read each line of th reducedFileA and put them into fileAList array of Strings
        TreeMap<String, ArrayList<String>> fileATreeMap = new TreeMap<>();

        Hashtable<String, Integer> reduceFileAHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> reduceFileAHashTableIndexReversed = new Hashtable<>();

        Hashtable<String, Integer> reduceFileBHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> reduceFileBHashTableIndexReversed = new Hashtable<>();

        try (GZIPInputStream reduceGz = new GZIPInputStream(new FileInputStream(reduceFileA));
                Reader decoder = new InputStreamReader(reduceGz);
                BufferedReader br = new BufferedReader(decoder)) {

            // FileReader fr = new FileReader(reduceFileA);
            // BufferedReader br = new BufferedReader(fr);
            String line = "";
            // First: read the header and avoid it
            line = br.readLine();

            // We do not use the previous line, instead, we use a predefined header

            if (type.equals("filtered")) {
                if (chrS.equals("23")) {
                    line = HEADER_MIXED_X;
                } else {
                    line = HEADER_MIXED;
                }
            } else if (type.equals("condensed")) {
                line = "chr\tposition\talleleA\talleleB\tpvalue\tinfo_all";
            }

            // By default values for indexes in the header
            // int indexRsId =2;
            int indexChr = 0;
            int indexPosition = 0;
            int indexAlleleA = 0;
            int indexAlleleB = 0;

            if (line != null && !line.isEmpty()) {
                reduceFileAHashTableIndex = createHashWithHeader(line, "\t");
                reduceFileAHashTableIndexReversed = createHashWithHeaderReversed(line, "\t");

                indexChr = reduceFileAHashTableIndex.get("chr");
                indexPosition = reduceFileAHashTableIndex.get("position");
                indexAlleleA = reduceFileAHashTableIndex.get("alleleA");
                indexAlleleB = reduceFileAHashTableIndex.get("alleleB");
            }

            while ((line = br.readLine()) != null) {
                ArrayList<String> fileAList = new ArrayList<>();
                // delimiter I assume a tap.
                String[] splitted = line.split("\t");

                // Store Position:Store rsIDCases:Store infoCases:Store certCases
                fileAList.add(line);
                String variantKey = splitted[indexChr] + "_" + splitted[indexPosition] + "_" + splitted[indexAlleleA] + "_"
                        + splitted[indexAlleleB];

                // We only store the ones that are not repeated.
                // Question for Siliva: what is the criteria?
                if (!fileATreeMap.containsKey(variantKey)) {
                    // Now, we put this casesList into the treemap with the key position
                    fileATreeMap.put(variantKey, fileAList);
                } else {
                    fileATreeMap.remove(variantKey);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // if (debug) {
        // long freeMemory = Runtime.getRuntime().freeMemory()/1048576;
        // long totalMemory = Runtime.getRuntime().totalMemory()/1048576;
        // long maxMemory = Runtime.getRuntime().maxMemory()/1048576;

        // System.out.println("JVM freeMemory: " + freeMemory);
        // System.out.println("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
        // System.out.println("JVM maxMemory also equals to maximum heap size of JVM : " + maxMemory);
        // }

        // We read each line of th reducedFileB and put them into fileBList array of Strings
        TreeMap<String, ArrayList<String>> fileBTreeMap = new TreeMap<>();

        try (GZIPInputStream reduceGz = new GZIPInputStream(new FileInputStream(reduceFileB));
                Reader decoder = new InputStreamReader(reduceGz);
                BufferedReader br = new BufferedReader(decoder)) {

            int indexChr = 0;
            int indexPosition = 0;
            int indexAlleleA = 0;
            int indexAlleleB = 0;

            // First: read the header and avoid it
            String line = null;
            String line2 = br.readLine();
            if (line != null && !line.isEmpty()) {
                reduceFileBHashTableIndex = createHashWithHeader(line, "\t");
                reduceFileBHashTableIndexReversed = createHashWithHeaderReversed(line, "\t");

                indexChr = reduceFileAHashTableIndex.get("chr");
                indexPosition = reduceFileAHashTableIndex.get("position");
                indexAlleleA = reduceFileAHashTableIndex.get("alleleA");
                indexAlleleB = reduceFileAHashTableIndex.get("alleleB");
            }

            while ((line = br.readLine()) != null) {
                ArrayList<String> fileBList = new ArrayList<>();
                String[] splitted = line.split("\t");// delimiter I assume single space.

                fileBList.add(line);
                String variantKey = splitted[indexChr] + "_" + splitted[indexPosition] + "_" + splitted[indexAlleleA] + "_"
                        + splitted[indexAlleleB];

                // We only store the ones that are not repeated.
                // Question for Siliva: what is the criteria?
                if (!fileBTreeMap.containsKey(variantKey)) {
                    // Now, we put this casesList into the treemap with the key position
                    fileBTreeMap.put(variantKey, fileBList);
                } else {
                    fileATreeMap.remove(variantKey);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // A place to store the results of this merge
        TreeMap<String, ArrayList<String>> fileCTreeMap = new TreeMap<>();
        Set<Entry<String, ArrayList<String>>> mySet = fileATreeMap.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();

        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, ArrayList<String>> m = iter.next();
            // getKey is used to get key of Map
            // position=(Integer)m.getKey();
            String variantKey = (String) m.getKey();

            // getValue is used to get value of key in Map
            ArrayList<String> fileTmp = new ArrayList<>();
            fileTmp = (ArrayList<String>) m.getValue();
            // We look for the casesPosition key in the controlsTreeMap.
            // If found, we get the value, otherwise we get null

            fileCTreeMap.put(variantKey, fileTmp);
        }
        mySet = fileBTreeMap.entrySet();
        // Move next key and value of Map by iterator
        iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, ArrayList<String>> m = iter.next();
            // getKey is used to get key of Map
            // position=(Integer)m.getKey();
            String variantKey = (String) m.getKey();

            // getValue is used to get value of key in Map
            ArrayList<String> fileTmp = new ArrayList<String>();
            fileTmp = (ArrayList<String>) m.getValue();
            // We look for the casesPosition key in the controlsTreeMap.
            // If found, we get the value, otherwise we get null
            fileCTreeMap.put(variantKey, fileTmp);
        }

        // Finally we put the fileCTreeMap into the outputFile
        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File fCombination = new File(reduceFileC);
        try {
            fCombination.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fCombination))) {

            String valueReversed = null;
            int index;
            if (reduceFileAHashTableIndexReversed.size() >= reduceFileBHashTableIndexReversed.size()) {
                for (index = 0; index < reduceFileAHashTableIndexReversed.size() - 1; index++) {
                    valueReversed = reduceFileAHashTableIndexReversed.get(index);
                    writer.write(valueReversed + "\t");
                }
                valueReversed = reduceFileAHashTableIndexReversed.get(index);
                writer.write(valueReversed);
                writer.newLine();
            } else {
                for (index = 0; index < reduceFileBHashTableIndexReversed.size() - 1; index++) {
                    valueReversed = reduceFileBHashTableIndexReversed.get(index);
                    writer.write(valueReversed + "\t");
                }
                valueReversed = reduceFileBHashTableIndexReversed.get(index);
                writer.write(valueReversed);
                writer.newLine();
            }

            mySet = fileCTreeMap.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                // getKey is used to get key of Map
                // int pos=(Integer)m.getKey();

                ArrayList<String> lineTmp = new ArrayList<>();
                lineTmp = (ArrayList<String>) m.getValue();
                writer.write(lineTmp.get(0));
                for (int j = 1; j < lineTmp.size(); j++) {
                    writer.write("\t" + lineTmp.get(j));
                }
                writer.newLine();
            }

            // if (debug) {
            // long freeMemory = Runtime.getRuntime().freeMemory()/1048576;
            // long totalMemory = Runtime.getRuntime().totalMemory()/1048576;
            // long maxMemory = Runtime.getRuntime().maxMemory()/1048576;

            // System.out.println("JVM freeMemory: " + freeMemory);
            // System.out.println("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
            // System.out.println("JVM maxMemory also equals to maximum heap size of JVM : " + maxMemory);
            // }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(reduceFileC, reduceFileC + ".gz");
        File fc = new File(reduceFileC);
        File fGz = new File(reduceFileC + ".gz");
        try {
            copyFile(fGz, fc);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] mergeTwoChunks startTime: " + startTime);
            System.out.println("\n[DEBUG] mergeTwoChunks endTime: " + stopTime);
            System.out.println("\n[DEBUG] mergeTwoChunks elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of mergeTwoChunks.");
        }
    }

    /**
     * Method to collect the summary
     * 
     * @param chr
     * @param firstImputeFileInfo
     * @param snptestOutFile
     * @param reduceFile
     * @param mafThresholdS
     * @param infoThresholdS
     * @param hweCohortThresholdS
     * @param hweCasesThresholdS
     * @param hweControlsThresholdS
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void collectSummary(String chr, String firstImputeFileInfo, String snptestOutFile, String reduceFile,
            String mafThresholdS, String infoThresholdS, String hweCohortThresholdS, String hweCasesThresholdS,
            String hweControlsThresholdS, String cmdToStore) throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running collectSummary with parameters:");
            System.out.println("[DEBUG] \t- Input chromosome             : " + chr);
            System.out.println("[DEBUG] \t- Input casesImputeFileInfo    : " + firstImputeFileInfo);
            System.out.println("[DEBUG] \t- Input snptestOutFile         : " + snptestOutFile);
            System.out.println("[DEBUG] \t- Output reduceFile            : " + reduceFile);
            System.out.println("[DEBUG] \t- Input mafThresholdS          : " + mafThresholdS);
            System.out.println("[DEBUG] \t- Input infoThresholdS         : " + infoThresholdS);
            System.out.println("[DEBUG] \t- Input hweCohortThresholdS    : " + hweCohortThresholdS);
            System.out.println("[DEBUG] \t- Input hweCasesThresholdS     : " + hweCasesThresholdS);
            System.out.println("[DEBUG] \t- Input hweControlsThresholdS  : " + hweControlsThresholdS);
            System.out.println("\n");
            System.out.println("[DEBUG] \t- COMMAND         : " + cmdToStore);

            // long freeMemory = Runtime.getRuntime().freeMemory()/1048576;
            // long totalMemory = Runtime.getRuntime().totalMemory()/1048576;
            // long maxMemory = Runtime.getRuntime().maxMemory()/1048576;

            // System.out.println("JVM freeMemory: " + freeMemory);
            // System.out.println("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
            // System.out.println("JVM maxMemory also equals to maximum heap size of JVM : " + maxMemory);

            // Map<String, String> env = System.getenv();
            // System.out.println("--------------------------------------");
            // System.out.println("Environmental Variables in Master:");
            // for (String envName : env.keySet()) {
            // System.out.format("%s=%s%n",envName,env.get(envName));
            // }
        }

        long startTime = System.currentTimeMillis();
        // Indexes for impute 2.3.2
        /*
         * int indexRsId = 1; int indexPosition = 2; int indexInfo = 6; int indexCertainty= 7;
         */

        int indexRsId = 0;
        int indexPosition = 0;
        int indexInfo = 0;
        int indexCertainty = 0;

        int length_entry_assoc_list = 0;

        // genOrBimFile is the mixedGenFile and
        // is not necessary to process it.

        // A place to store the results of the merge
        TreeMap<String, ArrayList<String>> summaryTotal = new TreeMap<>();
        Hashtable<String, Integer> snptestHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> snptestHashTableIndexReversed = new Hashtable<>();

        // We read each line of the firstImputeFileInfo and put them into firstList array of Strings
        TreeMap<String, ArrayList<String>> firstTreeMap = new TreeMap<>();
        String positionAndRsId = null;

        try (FileReader fr = new FileReader(firstImputeFileInfo); BufferedReader br = new BufferedReader(fr)) {
            String line = "";
            // Read the header and avoid the header
            line = br.readLine();
            // ---->
            Hashtable<String, Integer> imputeHashTableIndex = new Hashtable<>();
            if (line != null && !line.isEmpty()) {
                // System.err.println("[collectSummary] IMPRIMO line: " + line);
                // If we are here, the file is not empty.
                imputeHashTableIndex = createHashWithHeader(line, " ");

                indexPosition = imputeHashTableIndex.get("position");
                indexRsId = imputeHashTableIndex.get("rs_id");
                indexInfo = imputeHashTableIndex.get("info");
                indexCertainty = imputeHashTableIndex.get("certainty");

                // System.out.println("indexPosition: "+indexPosition);
                // System.out.println("indexRsId: " + indexRsId);
                // System.out.println("indexInfo: " + indexInfo);
                // System.out.println("indexCertainty: " + indexCertainty);
            }

            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter I assume single space.
                String[] splitted = line.split(" ");

                // Store Position:Store rsIDCases:Store infoCases:Store certCases
                firstList.add(splitted[indexPosition]);
                firstList.add(splitted[indexRsId]);
                firstList.add(splitted[indexInfo]);
                firstList.add(splitted[indexCertainty]);
                positionAndRsId = splitted[indexPosition] + "_" + splitted[indexRsId];

                // If there is not a previous snp with this combination of position and rsID, we store it.
                if (!firstTreeMap.containsKey(positionAndRsId)) {
                    // We, put this in the firstTreeMap
                    firstTreeMap.put(positionAndRsId, firstList);
                } else { // If there is a snp with this combination we should remove it.
                    firstTreeMap.remove(positionAndRsId);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // We read each line of the snptestOutFile and put them into assocTreeMap array of Strings
        TreeMap<String, ArrayList<String>> assocTreeMap = new TreeMap<>();
        try (GZIPInputStream snptestOutGz = new GZIPInputStream(new FileInputStream(snptestOutFile));
                Reader decoder = new InputStreamReader(snptestOutGz);
                BufferedReader br = new BufferedReader(decoder)) {

            String line = null;
            while ((line = br.readLine()) != null) {
                // We have to avoid reading all the comment lines that start with "#" character, and one more that is
                // the
                // header and start with "alternate".
                // System.err.println("[DEBUGING1]: " + line);
                String[] splitted = line.split(" ");
                String firstHeader = splitted[0];
                char firstChar = line.charAt(0);

                // TODO: the next if is ugly.
                if (firstHeader.equals("alternate_ids")) {
                    snptestHashTableIndex = createHashWithHeader(line, " ");
                    snptestHashTableIndexReversed = createHashWithHeaderReversed(line, " ");
                }

                if ((firstChar != '#') && (firstChar != 'a')) {
                    ArrayList<String> assocList = new ArrayList<String>();
                    // delimiter I assume single space.
                    // System.err.println("[DEBUGING1]: " + line);
                    // String[] splitted = line.split(" ");

                    // REVIEW: Should we store position?
                    // assocList.add(splitted[3]);

                    // We store everything, from the line.
                    int index_field = 0;
                    for (index_field = 0; index_field < splitted.length; index_field++) {
                        assocList.add(splitted[index_field]);
                    }

                    // Now, store the array of string assocList in the assocTreeMap
                    positionAndRsId = splitted[snptestHashTableIndex.get("position")] + "_" + splitted[snptestHashTableIndex.get("rsid")];
                    if (!assocTreeMap.containsKey(positionAndRsId)) {
                        assocTreeMap.put(positionAndRsId, assocList);
                        length_entry_assoc_list = assocList.size();
                    } else {
                        assocTreeMap.remove(positionAndRsId);
                    }
                } // The line does not start with "#" niether "alternate"
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Now we process all the SNP information that are common
        // Get TreeMap in Set interface to get key and value
        Set<Entry<String, ArrayList<String>>> mySet = firstTreeMap.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, ArrayList<String>> m = iter.next();
            // getKey is used to get key of Map
            String firstPositionAndRsId = (String) m.getKey();

            // getValue is used to get value of key in Map
            ArrayList<String> firstTmp = m.getValue();
            // The same for assocTreeMap. If found, we get the value, otherwise we get null
            ArrayList<String> assocTmp = new ArrayList<>();
            assocTmp = assocTreeMap.get(firstPositionAndRsId);
            ArrayList<String> summaryTmp = new ArrayList<>();
            try {
                summaryTmp = mergeArrays(chr, firstTmp, assocTmp, length_entry_assoc_list, mafThresholdS, infoThresholdS,
                        hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
            summaryTotal.put(firstPositionAndRsId, summaryTmp);
            assocTreeMap.remove(firstPositionAndRsId);
        }

        // ---------------
        // Finally we put the summaryTotal into the outputFile
        // We have to create the outputFile for this combination:
        // We verify that a file with the same name does not exist.
        File fCombination = new File(reduceFile);
        try {
            fCombination.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(fCombination))) {
            writer.write("chr\tposition\trs_id_all\tinfo_all\tcertainty_all\t");
            // We do not store the first 4 field because they are not necessary or are repeated:
            // These four fields are:
            // alternative_ids, rsid, chromosome, position
            for (int index = 4; index < snptestHashTableIndexReversed.size(); index++) {
                String valueReversed = snptestHashTableIndexReversed.get(index);
                writer.write(valueReversed + "\t");
            }
            writer.write("\n");

            mySet = summaryTotal.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                ArrayList<String> lineTmp = m.getValue();
                // System.out.println("VALUE = " + lineTmp);

                writer.write(lineTmp.get(0));
                for (int j = 1; j < lineTmp.size(); j++) {
                    writer.write("\t" + lineTmp.get(j));
                }
                writer.newLine();
            }

            // if (debug) {
            // long freeMemory = Runtime.getRuntime().freeMemory()/1048576;
            // long totalMemory = Runtime.getRuntime().totalMemory()/1048576;
            // long maxMemory = Runtime.getRuntime().maxMemory()/1048576;

            // System.out.println("JVM freeMemory: " + freeMemory);
            // System.out.println("JVM totalMemory also equals to initial heap size of JVM : " + totalMemory);
            // System.out.println("JVM maxMemory also equals to maximum heap size of JVM : " + maxMemory);
            // }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(reduceFile, reduceFile + ".gz");
        File fb = new File(reduceFile);
        File fGz = new File(reduceFile + ".gz");
        try {
            copyFile(fGz, fb);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] collectSummary startTime: " + startTime);
            System.out.println("\n[DEBUG] collectSummary endTime: " + stopTime);
            System.out.println("\n[DEBUG] collectSummary elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of collectSummary.");
        }
    }

    /**
     * Method to initialize the pheno matrix
     * 
     * @param topHitsFile
     * @param ttName
     * @param rpName
     * @param phenomeFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void initPhenoMatrix(String topHitsFile, String ttName, String rpName, String phenomeFile, String cmdToStore)
            throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running initPhenoMatrix with parameters:");
            System.out.println("[DEBUG] \t- Input topHitsFile       : " + topHitsFile);
            System.out.println("[DEBUG] \t- Input ttName            : " + ttName);
            System.out.println("[DEBUG] \t- Input rpName            : " + rpName);
            System.out.println("[DEBUG] \t- Output phenomeFile      : " + phenomeFile);
            System.out.println("\n");

            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // First, let's create the header for this file

        String headerPhenomeFile = "chr\tposition";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "rs_id_all";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "alleleA";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "alleleB";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "all_maf";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_pvalue";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_beta_1";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
        // "frequentist_add_beta_1:genotype/sex=1";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
        // "frequentist_add_beta_2:genotype/sex=2";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_se_1";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
        // "frequentist_add_se_1:genotype/sex=1";
        // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
        // "frequentist_add_se_2:genotype/sex=2";

        // Then, read the input file

        String chrAndPosition = null;
        // TODO HERE
        // replace commas in the header variable by a TAB
        String newHeader = headerPhenomeFile;
        // String newHeader = headerPhenomeFile.replace(',', '\t');
        if (DEBUG) {
            System.out.println("[DEBUG] \t- The new header will be : [" + newHeader + "]");
        }

        // First, we load the whole topHitsFile into a TreeMap

        // A place to store the results of the merge
        TreeMap<String, ArrayList<String>> phenomeTreeMap = new TreeMap<>();

        // Now we create the header for the phenomeFile. The label variable should contain all the information
        // to do it, as follows (//Everything in one line)
        // chr:position:rsId:pval_phe_[0]_pan_[0]:pval_phe_[0]_pa_[1]:...:pval_phe_[0]_pa_[K-1]:
        // pval_phe_[1]_pan_[0]:pval_phe_[1]_pa_[1]:...:pval_phe_[1]_pa_[K-1]:
        // pval_phe_[2]_pan_[0]:pval_phe_[2]_pa_[1]:...:pval_phe_[2]_pa_[K-1]:
        // ...
        // pval_phe_[T-1]_pan_[0]:pval_phe_[T-1]_pa_[1]:...:pval_phe_(T-1)_pa_[K-1]:
        // Where:
        // T : number of phenotypes.
        // phe_[i]: name of phenotype i.
        // K : number of panels.
        // pa_[j]: name of panel j.
        Hashtable<String, Integer> phenomeHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> phenomeHashTableIndexReversed = new Hashtable<>();

        phenomeHashTableIndex = createHashWithHeader(newHeader, "\t");

        // We start reading the topHits File
        try (GZIPInputStream topHitsFileGz = new GZIPInputStream(new FileInputStream(topHitsFile));
                Reader decoder = new InputStreamReader(topHitsFileGz);
                BufferedReader br = new BufferedReader(decoder);) {
            String line = br.readLine();
            Hashtable<String, Integer> topHitsHashTableIndex = new Hashtable<>();
            topHitsHashTableIndex = createHashWithHeader(line, "\t");

            int indexChrInTopHitsFile = topHitsHashTableIndex.get("chr");
            int indexPositionInTopHitsFile = topHitsHashTableIndex.get("position");

            chrAndPosition = null;
            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInTopHitsFile] + "_" + splitted[indexPositionInTopHitsFile];

                // Store chr:position:rsId:pvalues for the all combination of phenotypes and panels
                // It means
                firstList.add(splitted[indexChrInTopHitsFile]);
                firstList.add(splitted[indexPositionInTopHitsFile]);

                // System.out.println("\n[DEBUG] phenomeHashTableIndex.size() " + phenomeHashTableIndex.size());

                // Finally, we put this data into the firstTreeMap, using chrPosition as key and the firstList as value.
                phenomeTreeMap.put(chrAndPosition, firstList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Finally, we print the phenomeThreeMap into the output file.
        // We start with the header:
        File outPhenomeFile = new File(phenomeFile);

        // Try to create the file
        boolean bool = false;
        try {
            bool = outPhenomeFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about the existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outPhenomeFile + " was succesfuly created? " + bool);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outPhenomeFile))) {
            writer.write(newHeader);
            writer.newLine();

            Set<Entry<String, ArrayList<String>>> mySet = phenomeTreeMap.entrySet();
            // Move next key and value of Map by iterator
            Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                // getKey is used to get key of Map
                ArrayList<String> firstTmp = m.getValue();

                writer.write(firstTmp.get(0));
                for (int j = 1; j < firstTmp.size(); j++) {
                    writer.write("\t" + firstTmp.get(j));
                }
                // writer.write(myLine);
                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(phenomeFile, phenomeFile + ".gz");
        File fb = new File(phenomeFile);
        File fGz = new File(phenomeFile + ".gz");
        try {
            copyFile(fGz, fb);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] initPhenoMatrix startTime: " + startTime);
            System.out.println("\n[DEBUG] initPhenoMatrix endTime: " + stopTime);
            System.out.println("\n[DEBUG] initPhenoMatrix elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of initPhenoMatrix.");
        }
    }

    /**
     * Method to add to the pheno matrix
     * 
     * @param phenomeAFile
     * @param topHitsFile
     * @param ttName
     * @param rpName
     * @param phenomeBFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void addToPhenoMatrix(String phenomeAFile, String topHitsFile, String ttName, String rpName, String phenomeBFile,
            String cmdToStore) throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running addToPhenoMatrix with parameters:");
            System.out.println("[DEBUG] \t- Input phenomeFileA      : " + phenomeAFile);
            System.out.println("[DEBUG] \t- Input topHitsFile       : " + topHitsFile);
            System.out.println("[DEBUG] \t- Input ttName            : " + ttName);
            System.out.println("[DEBUG] \t- Input rpName            : " + rpName);
            System.out.println("[DEBUG] \t- Output phenomeFileB     : " + phenomeBFile);
            System.out.println("\n");

            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }

        long startTime = System.currentTimeMillis();

        // A place to store the results
        TreeMap<String, ArrayList<String>> phenomeATreeMap = new TreeMap<>();

        // Now we create the header for the phenomeFile. The label variable should contain all the information
        // to do it, as follows (//Everything in one line)
        // chr:position:rsId:pval_phe_[0]_pan_[0]:pval_phe_[0]_pa_[1]:...:pval_phe_[0]_pa_[K-1]:
        // pval_phe_[1]_pan_[0]:pval_phe_[1]_pa_[1]:...:pval_phe_[1]_pa_[K-1]:
        // pval_phe_[2]_pan_[0]:pval_phe_[2]_pa_[1]:...:pval_phe_[2]_pa_[K-1]:
        // ...
        // pval_phe_[T-1]_pan_[0]:pval_phe_[T-1]_pa_[1]:...:pval_phe_(T-1)_pa_[K-1]:
        // Where:
        // T : number of phenotypes.
        // phe_[i]: name of phenotype i.
        // K : number of panels.
        // pa_[j]: name of panel j.
        Hashtable<String, Integer> phenomeAHashTableIndex = new Hashtable<>();

        // We start reading the phenomeFileA
        String phenomeAHeader = null;
        try (GZIPInputStream phenomeAFileGz = new GZIPInputStream(new FileInputStream(phenomeAFile));
                Reader decoder = new InputStreamReader(phenomeAFileGz);
                BufferedReader br = new BufferedReader(decoder)) {
            // First of all, the header
            phenomeAHeader = br.readLine();

            phenomeAHeader = "chr\tposition";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "rs_id_all";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "alleleA";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "alleleB";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "all_maf";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_pvalue";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_beta_1";
            // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
            // "frequentist_add_beta_1:genotype/sex=1";
            // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
            // "frequentist_add_beta_2:genotype/sex=2";
            // phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_se_1";
            // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
            // "frequentist_add_se_1:genotype/sex=1";
            // headerPhenomeFile = headerPhenomeFile + "\t" + ttName + ":" + rpName + ":" +
            // "frequentist_add_se_2:genotype/sex=2";

            phenomeAHashTableIndex = createHashWithHeader(phenomeAHeader, "\t");

            int indexChrInPhenomeAFile = phenomeAHashTableIndex.get("chr");
            int indexPositionInPhenomeAFile = phenomeAHashTableIndex.get("position");

            // Then, we read the rest of the file and put the information into the phenomeATreeMap
            String chrAndPosition = null;
            String line = null;

            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter,I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInPhenomeAFile] + "_" + splitted[indexPositionInPhenomeAFile];

                // Store chr:position
                firstList.add(splitted[indexChrInPhenomeAFile]);
                firstList.add(splitted[indexPositionInPhenomeAFile]);

                // Finally, we put this data into the phenomeATreeMap, using chrPosition as key and the firstList as
                // value.
                phenomeATreeMap.put(chrAndPosition, firstList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Now we read the topHitsFile
        try (GZIPInputStream topHitsFileGz = new GZIPInputStream(new FileInputStream(topHitsFile));
                Reader decoder = new InputStreamReader(topHitsFileGz);
                BufferedReader br = new BufferedReader(decoder);) {

            // We start reading the topHits File
            String line = br.readLine();
            Hashtable<String, Integer> topHitsHashTableIndex = new Hashtable<>();
            topHitsHashTableIndex = createHashWithHeader(line, "\t");

            int indexChrInTopHitsFile = topHitsHashTableIndex.get("chr");
            int indexPositionInTopHitsFile = topHitsHashTableIndex.get("position");

            String chrAndPosition = null;
            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInTopHitsFile] + "_" + splitted[indexPositionInTopHitsFile];

                if (!phenomeATreeMap.containsKey(chrAndPosition)) {
                    firstList.add(splitted[indexChrInTopHitsFile]);
                    firstList.add(splitted[indexPositionInTopHitsFile]);

                    // Finally, we put this data into the phenomeATreeMap, using chrPosition as key and the firstList as
                    // value.
                    phenomeATreeMap.put(chrAndPosition, firstList);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Finally, we print the phenomeThreeMap into the output file.
        // We start with the header:
        File outPhenomeFile = new File(phenomeBFile);

        // Try to create the file
        boolean bool = false;
        try {
            bool = outPhenomeFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about the existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outPhenomeFile + " was succesfuly created? " + bool);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outPhenomeFile))) {
            writer.write(phenomeAHeader);
            writer.newLine();

            Set<Entry<String, ArrayList<String>>> mySet = phenomeATreeMap.entrySet();
            // Move next key and value of Map by iterator
            Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                // getKey is used to get key of Map
                ArrayList<String> firstTmp = m.getValue();

                writer.write(firstTmp.get(0));
                for (int j = 1; j < firstTmp.size(); j++) {
                    writer.write("\t" + firstTmp.get(j));
                }
                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(phenomeBFile, phenomeBFile + ".gz");
        File fb = new File(phenomeBFile);
        File fGz = new File(phenomeBFile + ".gz");
        try {
            copyFile(fGz, fb);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] addToPhenoMatrix startTime: " + startTime);
            System.out.println("\n[DEBUG] addToPhenoMatrix endTime: " + stopTime);
            System.out.println("\n[DEBUG] addToPhenoMatrix elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of addToPhenoMatrix.");
        }
    }

    /**
     * Method to add X chromosome to the pheno matrix
     * 
     * @param phenomeAFile
     * @param filteredByAllFile
     * @param ttName
     * @param rpName
     * @param phenomeBFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void addToPhenoMatrixX(String phenomeAFile, String filteredByAllFile, String ttName, String rpName, String phenomeBFile,
            String cmdToStore) throws GuidanceTaskException {
        if (DEBUG) {
            System.out.println("\n[DEBUG] Running addToPhenoMatrixX with parameters:");
            System.out.println("[DEBUG] \t- Input phenomeFileA      : " + phenomeAFile);
            System.out.println("[DEBUG] \t- Input filteredByAllFile : " + filteredByAllFile);
            System.out.println("[DEBUG] \t- Input ttName            : " + ttName);
            System.out.println("[DEBUG] \t- Input rpName            : " + rpName);
            System.out.println("[DEBUG] \t- Output phenomeFileB     : " + phenomeBFile);
            System.out.println("\n");

            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // A place to store the results
        TreeMap<String, ArrayList<String>> phenomeATreeMap = new TreeMap<>();

        // Now we create the header for the phenomeFile. The label variable should contain all the information
        // to do it, as follows (//Everything in one line)
        // chr:position:rsId:pval_phe_[0]_pan_[0]:pval_phe_[0]_pa_[1]:...:pval_phe_[0]_pa_[K-1]:
        // pval_phe_[1]_pan_[0]:pval_phe_[1]_pa_[1]:...:pval_phe_[1]_pa_[K-1]:
        // pval_phe_[2]_pan_[0]:pval_phe_[2]_pa_[1]:...:pval_phe_[2]_pa_[K-1]:
        // ...
        // pval_phe_[T-1]_pan_[0]:pval_phe_[T-1]_pa_[1]:...:pval_phe_(T-1)_pa_[K-1]:
        // Where:
        // T : number of phenotypes.
        // phe_[i]: name of phenotype i.
        // K : number of panels.
        // pa_[j]: name of panel j.
        Hashtable<String, Integer> phenomeAHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> phenomeAHashTableIndexReversed = new Hashtable<>();

        // We start reading the phenomeFileA
        // First of all, the header
        String phenomeAHeader = null;
        try (GZIPInputStream phenomeAFileGz = new GZIPInputStream(new FileInputStream(phenomeAFile));
                Reader decoder = new InputStreamReader(phenomeAFileGz);
                BufferedReader br = new BufferedReader(decoder);) {
            phenomeAHeader = br.readLine();
            phenomeAHashTableIndex = createHashWithHeader(phenomeAHeader, "\t");

            int indexChrInPhenomeAFile = phenomeAHashTableIndex.get("chr");
            int indexPositionInPhenomeAFile = phenomeAHashTableIndex.get("position");

            // Then, we read the rest of the file and put the information into the phenomeATreeMap
            String chrAndPosition = null;
            String line = null;
            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter,I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInPhenomeAFile] + "_" + splitted[indexPositionInPhenomeAFile];

                // Store chr:position
                firstList.add(splitted[indexChrInPhenomeAFile]);
                firstList.add(splitted[indexPositionInPhenomeAFile]);

                // Finally, we put this data into the phenomeATreeMap, using chrPosition as key and the firstList as
                // value.
                phenomeATreeMap.put(chrAndPosition, firstList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // ------------->>>

        // Now we load the whole filteredByAllFile into a TreeMap
        TreeMap<String, ArrayList<String>> filteredTreeMap = new TreeMap<>();
        try (GZIPInputStream filteredByAllGz = new GZIPInputStream(new FileInputStream(filteredByAllFile));
                Reader decoder = new InputStreamReader(filteredByAllGz);
                BufferedReader br = new BufferedReader(decoder)) {

            Hashtable<String, Integer> filteredByAllHashTableIndex = new Hashtable<>();
            Hashtable<Integer, String> filteredByAllHashTableIndexReversed = new Hashtable<>();

            // We start reading the filteredByAllFile
            // First of all, the header
            String filteredHeader = br.readLine();
            filteredByAllHashTableIndex = createHashWithHeader(filteredHeader, "\t");

            int indexChrInFiltered = filteredByAllHashTableIndex.get("chr");
            int indexPositionInFiltered = filteredByAllHashTableIndex.get("position");
            int indexRsIdInFiltered = filteredByAllHashTableIndex.get("rs_id_all");
            int indexAlleleAInFiltered = filteredByAllHashTableIndex.get("alleleA");
            int indexAlleleBInFiltered = filteredByAllHashTableIndex.get("alleleB");
            int indexAllMafInFiltered = filteredByAllHashTableIndex.get("all_maf");
            int indexFreqAddPvalueInFiltered = filteredByAllHashTableIndex.get("frequentist_add_pvalue");

            int indexFreqAddBetaInFiltered = filteredByAllHashTableIndex.get("frequentist_add_beta_1");
            // int indexFreqAddBeta1sex1InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_beta_1:genotype/sex=1");
            // int indexFreqAddBeta2sex2InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_beta_2:genotype/sex=2");

            int indexFreqAddSeInFiltered = filteredByAllHashTableIndex.get("frequentist_add_se_1");
            // int indexFreqAddSe1sex1InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_se_1:genotype/sex=1");
            // int indexFreqAddSe2sex2InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_se_2:genotype/sex=2");

            String line = null;
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split("\t");
                String chrAndPosition = splitted[indexChrInFiltered] + "_" + splitted[indexPositionInFiltered];

                ArrayList<String> reducedList = new ArrayList<>();
                reducedList.add(splitted[indexRsIdInFiltered]);
                reducedList.add(splitted[indexAlleleAInFiltered]);
                reducedList.add(splitted[indexAlleleBInFiltered]);
                reducedList.add(splitted[indexAllMafInFiltered]);
                reducedList.add(splitted[indexFreqAddPvalueInFiltered]);

                reducedList.add(splitted[indexFreqAddBetaInFiltered]);
                // reducedList.add(splitted[indexFreqAddBeta1sex1InFiltered]);
                // reducedList.add(splitted[indexFreqAddBeta2sex2InFiltered]);

                reducedList.add(splitted[indexFreqAddSeInFiltered]);
                // reducedList.add(splitted[indexFreqAddSe1sex1InFiltered]);
                // reducedList.add(splitted[indexFreqAddSe2sex2InFiltered]);

                filteredTreeMap.put(chrAndPosition, reducedList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Finally, we print the phenomeThreeMap into the output file.
        // We start with the header:

        Set<Entry<String, ArrayList<String>>> mySet = phenomeATreeMap.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, ArrayList<String>> m = iter.next();
            // getKey is used to get key of Map
            String chrAndPosition = (String) m.getKey();

            ArrayList<String> currentList = (ArrayList<String>) m.getValue();
            ArrayList<String> reducedList = null;

            if (filteredTreeMap.containsKey(chrAndPosition)) {
                reducedList = filteredTreeMap.get(chrAndPosition);

                for (int i = 0; i < reducedList.size(); i++) {
                    currentList.add(reducedList.get(i));
                }
            } else {
                for (int i = 0; i < 11; i++) {
                    currentList.add("NA");
                }
            }

            filteredTreeMap.put(chrAndPosition, currentList);
        }

        // Finally, we print the phenomeThreeMap into the output file.
        // We start with the header:
        File outPhenomeFile = new File(phenomeBFile);

        // Try to create the file
        boolean bool = false;
        try {
            bool = outPhenomeFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about the existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outPhenomeFile + " was succesfuly created? " + bool);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outPhenomeFile))) {
            writer.write(phenomeAHeader);
            writer.newLine();

            mySet = phenomeATreeMap.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                // getKey is used to get key of Map
                ArrayList<String> firstTmp = m.getValue();

                writer.write(firstTmp.get(0));
                for (int j = 1; j < firstTmp.size(); j++) {
                    writer.write("\t" + firstTmp.get(j));
                }
                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(phenomeBFile, phenomeBFile + ".gz");
        File fb = new File(phenomeBFile);
        File fGz = new File(phenomeBFile + ".gz");
        try {
            copyFile(fGz, fb);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] addToPhenoMatrix startTime: " + startTime);
            System.out.println("\n[DEBUG] addToPhenoMatrix endTime: " + stopTime);
            System.out.println("\n[DEBUG] addToPhenoMatrix elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of addToPhenoMatrixX");
        }
    }

    /**
     * Method to fill out the pheno matrix
     * 
     * @param phenomeAFile
     * @param filteredByAllFile
     * @param filteredByAllXFile
     * @param endChrS
     * @param ttName
     * @param rpName
     * @param phenomeBFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void filloutPhenoMatrix(String phenomeAFile, String filteredByAllFile, String filteredByAllXFile, String endChrS,
            String ttName, String rpName, String phenomeBFile, String cmdToStore) throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running filloutPhenoMatrix with parameters:");
            System.out.println("[DEBUG] \t- Input phenomeFileA       : " + phenomeAFile);
            System.out.println("[DEBUG] \t- Input filteredByAllFile  : " + filteredByAllFile);
            System.out.println("[DEBUG] \t- Input filteredByAllXFile : " + filteredByAllXFile);
            System.out.println("[DEBUG] \t- Input endChrS            : " + endChrS);
            System.out.println("[DEBUG] \t- Input ttName             : " + ttName);
            System.out.println("[DEBUG] \t- Input rpName             : " + rpName);
            System.out.println("[DEBUG] \t- Output phenomeFileB      : " + phenomeBFile);
            System.out.println("\n");

            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // A place to store the data
        TreeMap<String, ArrayList<String>> phenomeATreeMap = new TreeMap<>();

        Hashtable<String, Integer> phenomeAHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> phenomeAHashTableIndexReversed = new Hashtable<>();

        // We start reading the phenomeFileA
        String phenomeAHeader = null;
        try (GZIPInputStream phenomeAFileGz = new GZIPInputStream(new FileInputStream(phenomeAFile));
                Reader decoder = new InputStreamReader(phenomeAFileGz);
                BufferedReader br = new BufferedReader(decoder)) {
            // First of all, the header
            phenomeAHeader = br.readLine();

            // phenomeAHeader = "chr\tposition";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "rs_id_all";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "alleleA";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "alleleB";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "all_maf";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_pvalue";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_beta_1";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_se_1";

            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_beta_1:genotype/sex=1";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_beta_2:genotype/sex=2";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_se_1:genotype/sex=1";
            phenomeAHeader = phenomeAHeader + "\t" + ttName + ":" + rpName + ":" + "frequentist_add_se_2:genotype/sex=2";

            phenomeAHashTableIndex = createHashWithHeader(phenomeAHeader, "\t");

            int indexChrInPhenomeAFile = phenomeAHashTableIndex.get("chr");
            int indexPositionInPhenomeAFile = phenomeAHashTableIndex.get("position");

            int indexRsIdInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "rs_id_all");
            int indexAlleleAInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "alleleA");
            int indexAlleleBInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "alleleB");
            int indexallMAFInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "all_maf");
            int indexFreqAddPvalueInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "frequentist_add_pvalue");
            int indexFreqAddBetaInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "frequentist_add_beta_1");
            int indexFreqAddSeInPhenomeAFile = phenomeAHashTableIndex.get(ttName + ":" + rpName + ":" + "frequentist_add_se_1");

            int indexFreqAddBeta1Sex1InPhenomeAFile = phenomeAHashTableIndex
                    .get(ttName + ":" + rpName + ":" + "frequentist_add_beta_1:genotype/sex=1");
            int indexFreqAddBeta2Sex2InPhenomeAFile = phenomeAHashTableIndex
                    .get(ttName + ":" + rpName + ":" + "frequentist_add_beta_2:genotype/sex=2");
            int indexFreqAddSe1Sex1InPhenomeAFile = phenomeAHashTableIndex
                    .get(ttName + ":" + rpName + ":" + "frequentist_add_se_1:genotype/sex=1");
            int indexFreqAddSe2Sex2InPhenomeAFile = phenomeAHashTableIndex
                    .get(ttName + ":" + rpName + ":" + "frequentist_add_se_2:genotype/sex=2");

            // --->
            // Then, we read the rest of the file and put the information into the phenomeATreeMap
            String chrAndPosition = null;
            String line = null;
            while ((line = br.readLine()) != null) {
                ArrayList<String> currentList = new ArrayList<>();
                // delimiter,I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInPhenomeAFile] + "_" + splitted[indexPositionInPhenomeAFile];

                // currentList.add(splitted[indexChrInPhenomeAFile]);
                // currentList.add(splitted[indexPositionInPhenomeAFile]);
                for (int i = 0; i < splitted.length; i++) {
                    currentList.add(splitted[i]);
                }

                // We update the phenomeATreeMap with the currentList
                phenomeATreeMap.put(chrAndPosition, currentList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        ////

        // Then, we need to extract the information of each snp from the filteredByAllFile
        // Now we load the whole filteredByAllFile into a TreeMap
        TreeMap<String, ArrayList<String>> filteredTreeMap = new TreeMap<>();
        try (GZIPInputStream filteredByAllGz = new GZIPInputStream(new FileInputStream(filteredByAllFile));
                Reader decoder = new InputStreamReader(filteredByAllGz);
                BufferedReader br = new BufferedReader(decoder)) {

            Hashtable<String, Integer> filteredByAllHashTableIndex = new Hashtable<>();
            Hashtable<Integer, String> filteredByAllHashTableIndexReversed = new Hashtable<>();

            // We start reading the filteredByAllFile
            // First of all, the header
            String filteredHeader = br.readLine();
            filteredByAllHashTableIndex = createHashWithHeader(filteredHeader, "\t");

            int indexChrInFiltered = filteredByAllHashTableIndex.get("chr");
            int indexPositionInFiltered = filteredByAllHashTableIndex.get("position");
            int indexRsIdInFiltered = filteredByAllHashTableIndex.get("rs_id_all");
            int indexAlleleAInFiltered = filteredByAllHashTableIndex.get("alleleA");
            int indexAlleleBInFiltered = filteredByAllHashTableIndex.get("alleleB");
            int indexAllMafInFiltered = filteredByAllHashTableIndex.get("all_maf");
            int indexFreqAddPvalueInFiltered = filteredByAllHashTableIndex.get("frequentist_add_pvalue");
            int indexFreqAddBetaInFiltered = filteredByAllHashTableIndex.get("frequentist_add_beta_1");
            // int indexFreqAddBeta1sex1InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_beta_1:genotype/sex=1");
            // int indexFreqAddBeta2sex2InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_beta_2:genotype/sex=2");
            int indexFreqAddSeInFiltered = filteredByAllHashTableIndex.get("frequentist_add_se_1");
            // int indexFreqAddSe1sex1InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_se_1:genotype/sex=1");
            // int indexFreqAddSe2sex2InFiltered =
            // filteredByAllHashTableIndex.get("frequentist_add_se_2:genotype/sex=2");

            String line = null;
            while ((line = br.readLine()) != null) {
                String[] splitted = line.split("\t");
                String chrAndPosition = splitted[indexChrInFiltered] + "_" + splitted[indexPositionInFiltered];

                ArrayList<String> reducedList = new ArrayList<>();
                reducedList.add(splitted[indexRsIdInFiltered]);
                reducedList.add(splitted[indexAlleleAInFiltered]);
                reducedList.add(splitted[indexAlleleBInFiltered]);
                reducedList.add(splitted[indexAllMafInFiltered]);
                reducedList.add(splitted[indexFreqAddPvalueInFiltered]);
                reducedList.add(splitted[indexFreqAddBetaInFiltered]);
                reducedList.add(splitted[indexFreqAddSeInFiltered]);

                // Now we put 4 values more that are the ones for chrX
                reducedList.add("NA");
                reducedList.add("NA");
                reducedList.add("NA");
                reducedList.add("NA");

                filteredTreeMap.put(chrAndPosition, reducedList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Here we have to do some similar with filteredByAllXFile (the results for chr23)
        TreeMap<String, ArrayList<String>> filteredXTreeMap = new TreeMap<>();
        Hashtable<String, Integer> filteredByAllXHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> filteredByAllXHashTableIndexReversed = new Hashtable<>();

        if (endChrS.equals("23")) {
            // Then, we need to extract the information of each snp from the filteredByAllFile
            // Now we load the whole filteredByAllFile into a TreeMap
            try (GZIPInputStream filteredByAllXGz = new GZIPInputStream(new FileInputStream(filteredByAllXFile));
                    Reader decoder = new InputStreamReader(filteredByAllXGz);
                    BufferedReader br = new BufferedReader(decoder)) {

                // We start reading the filteredByAllXFile
                // First of all, the header
                String filteredXHeader = br.readLine();
                filteredByAllXHashTableIndex = createHashWithHeader(filteredXHeader, "\t");

                int indexChrInFilteredX = filteredByAllXHashTableIndex.get("chr");
                int indexPositionInFilteredX = filteredByAllXHashTableIndex.get("position");
                int indexRsIdInFilteredX = filteredByAllXHashTableIndex.get("rs_id_all");
                int indexAlleleAInFilteredX = filteredByAllXHashTableIndex.get("alleleA");
                int indexAlleleBInFilteredX = filteredByAllXHashTableIndex.get("alleleB");
                int indexAllMafInFilteredX = filteredByAllXHashTableIndex.get("all_maf");
                int indexFreqAddPvalueInFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_pvalue");
                // int indexFreqAddBetaInFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_beta_1");
                // int indexFreqAddSeInFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_se_1");

                int indexFreqAddBeta1sex1InFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_beta_1:genotype/sex=1");
                int indexFreqAddBeta2sex2InFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_beta_2:genotype/sex=2");
                int indexFreqAddSe1sex1InFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_se_1:genotype/sex=1");
                int indexFreqAddSe2sex2InFilteredX = filteredByAllXHashTableIndex.get("frequentist_add_se_2:genotype/sex=2");

                String line = null;
                while ((line = br.readLine()) != null) {
                    String[] splitted = line.split("\t");
                    String chrAndPosition = splitted[indexChrInFilteredX] + "_" + splitted[indexPositionInFilteredX];

                    ArrayList<String> reducedList = new ArrayList<>();
                    reducedList.add(splitted[indexRsIdInFilteredX]);
                    reducedList.add(splitted[indexAlleleAInFilteredX]);
                    reducedList.add(splitted[indexAlleleBInFilteredX]);
                    reducedList.add(splitted[indexAllMafInFilteredX]);
                    reducedList.add(splitted[indexFreqAddPvalueInFilteredX]);

                    // This to values for chr23
                    reducedList.add("NA");
                    reducedList.add("NA");

                    // Now we put 4 values more that are the ones for chrX
                    reducedList.add(splitted[indexFreqAddBeta1sex1InFilteredX]);
                    reducedList.add(splitted[indexFreqAddBeta2sex2InFilteredX]);
                    reducedList.add(splitted[indexFreqAddSe1sex1InFilteredX]);
                    reducedList.add(splitted[indexFreqAddSe2sex2InFilteredX]);

                    filteredXTreeMap.put(chrAndPosition, reducedList);
                }
            } catch (IOException ioe) {
                throw new GuidanceTaskException(ioe);
            }
        }

        // Now, we print the information of each snp from filteredByAllFile into the phenomeATreeMap

        Set<Entry<String, ArrayList<String>>> mySet = phenomeATreeMap.entrySet();
        // Move next key and value of Map by iterator
        Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();
        while (iter.hasNext()) {
            // key=value separator this by Map.Entry to get key and value
            Entry<String, ArrayList<String>> m = iter.next();

            ArrayList<String> currentList = m.getValue();
            ArrayList<String> reducedList = null;

            // getKey is used to get key of Map
            String chrAndPosition = (String) m.getKey();

            String[] divideKey = chrAndPosition.split("_");
            String Chr = divideKey[0];

            if (!Chr.equals("23")) {
                if (filteredTreeMap.containsKey(chrAndPosition)) {
                    reducedList = filteredTreeMap.get(chrAndPosition);

                    for (int i = 0; i < reducedList.size(); i++) {
                        currentList.add(reducedList.get(i));
                    }
                } else {
                    for (int i = 0; i < 11; i++) {
                        currentList.add("NA");
                    }
                }
            } else {
                if (filteredXTreeMap.containsKey(chrAndPosition)) {
                    reducedList = filteredXTreeMap.get(chrAndPosition);

                    for (int i = 0; i < reducedList.size(); i++) {
                        currentList.add(reducedList.get(i));
                    }
                } else {
                    for (int i = 0; i < 11; i++) {
                        currentList.add("NA");
                    }
                }
            }

            phenomeATreeMap.put(chrAndPosition, currentList);
        }

        // Finally, we print the phenomeAThreeMap into the output file.
        // We start with the header:
        File outPhenomeFile = new File(phenomeBFile);

        // Try to create the file
        boolean bool = false;
        try {
            bool = outPhenomeFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about the existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outPhenomeFile + " was succesfuly created? " + bool);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outPhenomeFile))) {
            writer.write(phenomeAHeader);
            writer.newLine();

            mySet = phenomeATreeMap.entrySet();
            // Move next key and value of Map by iterator
            iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                // getKey is used to get key of Map
                ArrayList<String> firstTmp = m.getValue();

                writer.write(firstTmp.get(0));
                for (int j = 1; j < firstTmp.size(); j++) {
                    writer.write("\t" + firstTmp.get(j));
                }

                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(phenomeBFile, phenomeBFile + ".gz");
        File fb = new File(phenomeBFile);
        File fGz = new File(phenomeBFile + ".gz");
        try {
            copyFile(fGz, fb);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] filloutPhenoMatrix startTime: " + startTime);
            System.out.println("\n[DEBUG] filloutPhenoMatrix endTime: " + stopTime);
            System.out.println("\n[DEBUG] filloutPhenoMatrix elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of filloutPhenoMatrix.");
        }

    }

    /**
     * Method to finalize the pheno matrix
     * 
     * @param phenomeAFile
     * @param phenomeBFile
     * @param ttName
     * @param rpName
     * @param phenomeCFile
     * @param cmdToStore
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void finalizePhenoMatrix(String phenomeAFile, String phenomeBFile, String ttName, String rpName, String phenomeCFile,
            String cmdToStore) throws GuidanceTaskException {

        if (DEBUG) {
            System.out.println("\n[DEBUG] Running finalizePhenoMatrix with parameters:");
            System.out.println("[DEBUG] \t- Input phenomeAFile  : " + phenomeAFile);
            System.out.println("[DEBUG] \t- Input phenomeBFile  : " + phenomeBFile);
            System.out.println("[DEBUG] \t- Input ttName        : " + ttName);
            System.out.println("[DEBUG] \t- Input rpName        : " + rpName);
            System.out.println("[DEBUG] \t- Output phenomeCFile : " + phenomeCFile);
            System.out.println("\n");

            System.out.println("\n");
            System.out.println("[DEBUG] \t- Command: " + cmdToStore);
        }
        long startTime = System.currentTimeMillis();

        // A place to store the data
        TreeMap<String, ArrayList<String>> phenomeATreeMap = new TreeMap<>();
        Hashtable<String, Integer> phenomeAHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> phenomeAHashTableIndexReversed = new Hashtable<>();

        // We start reading the phenomeFileA
        String phenomeAHeader = null;
        try (GZIPInputStream phenomeAFileGz = new GZIPInputStream(new FileInputStream(phenomeAFile));
                Reader decoder = new InputStreamReader(phenomeAFileGz);
                BufferedReader br = new BufferedReader(decoder)) {
            // First of all, the header
            phenomeAHeader = br.readLine();
            phenomeAHashTableIndex = createHashWithHeader(phenomeAHeader, "\t");

            int indexChrInPhenomeAFile = phenomeAHashTableIndex.get("chr");
            int indexPositionInPhenomeAFile = phenomeAHashTableIndex.get("position");

            // Then, we read the rest of the file and put the information into the phenomeATreeMap
            String chrAndPosition = null;
            String line = null;
            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter,I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInPhenomeAFile] + "_" + splitted[indexPositionInPhenomeAFile];

                // Store chr:position
                for (int i = 0; i < splitted.length; i++) {
                    firstList.add(splitted[i]);
                }
                // Finally, we put this data into the phenomeATreeMap, using chrPosition as key and the firstList as
                // value.
                phenomeATreeMap.put(chrAndPosition, firstList);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Second we read the phenomeBFile and load the information into the phenomeATreeMap
        String phenomeBHeader = null;
        Hashtable<String, Integer> phenomeBHashTableIndex = new Hashtable<>();
        Hashtable<Integer, String> phenomeBHashTableIndexReversed = new Hashtable<>();
        try (GZIPInputStream phenomeBFileGz = new GZIPInputStream(new FileInputStream(phenomeBFile));
                Reader decoder = new InputStreamReader(phenomeBFileGz);
                BufferedReader br = new BufferedReader(decoder)) {

            // We start reading the phenomeFileA
            // First of all, the header
            phenomeBHeader = br.readLine();
            phenomeBHashTableIndex = createHashWithHeader(phenomeBHeader, "\t");
            phenomeBHashTableIndexReversed = createHashWithHeaderReversed(phenomeBHeader, "\t");

            int indexChrInPhenomeBFile = phenomeBHashTableIndex.get("chr");
            int indexPositionInPhenomeBFile = phenomeBHashTableIndex.get("position");

            // Then, we read the rest of the file and put the information into the phenomeATreeMap
            String chrAndPosition = null;
            String line = null;
            while ((line = br.readLine()) != null) {
                ArrayList<String> firstList = new ArrayList<>();
                // delimiter,I assume TAP space.
                String[] splitted = line.split("\t");

                chrAndPosition = splitted[indexChrInPhenomeBFile] + "_" + splitted[indexPositionInPhenomeBFile];

                if (phenomeATreeMap.containsKey(chrAndPosition)) {
                    firstList = phenomeATreeMap.get(chrAndPosition);

                    for (int i = 2; i < splitted.length; i++) {
                        firstList.add(splitted[i]);
                    }

                    // Finally, we put this data into the phenomeATreeMap, using chrPosition as key and the firstList as
                    // value.
                    phenomeATreeMap.put(chrAndPosition, firstList);
                }
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        String finalHeader = phenomeAHeader;
        for (int i = 2; i < phenomeBHashTableIndex.size(); i++) {
            finalHeader = finalHeader + "\t" + phenomeBHashTableIndexReversed.get(i);
        }

        // Finally, we print the phenomeAThreeMap into the output file.
        // We start with the header:
        File outPhenomeFile = new File(phenomeCFile);

        // Try to create the file
        boolean bool = false;
        try {
            bool = outPhenomeFile.createNewFile();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
        // Print information about the existence of the file
        System.out.println("\n[DEBUG] \t- Output file " + outPhenomeFile + " was succesfuly created? " + bool);

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outPhenomeFile))) {
            writer.write(finalHeader);
            writer.newLine();

            Set<Entry<String, ArrayList<String>>> mySet = phenomeATreeMap.entrySet();
            // Move next key and value of Map by iterator
            Iterator<Entry<String, ArrayList<String>>> iter = mySet.iterator();
            while (iter.hasNext()) {
                // key=value separator this by Map.Entry to get key and value
                Entry<String, ArrayList<String>> m = iter.next();
                // getKey is used to get key of Map
                ArrayList<String> firstTmp = m.getValue();

                writer.write(firstTmp.get(0));
                for (int j = 1; j < firstTmp.size(); j++) {
                    writer.write("\t" + firstTmp.get(j));
                }

                writer.newLine();
            }

            writer.flush();
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // Then, we create the gz file and rename it
        gzipFile(phenomeCFile, phenomeCFile + ".gz");
        File fb = new File(phenomeCFile);
        File fGz = new File(phenomeCFile + ".gz");
        try {
            copyFile(fGz, fb);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }

        // <----------
        long stopTime = System.currentTimeMillis();
        long elapsedTime = (stopTime - startTime) / 1000;
        if (DEBUG) {
            System.out.println("\n[DEBUG] finalizePhenoMatrix startTime: " + startTime);
            System.out.println("\n[DEBUG] finalizePhenoMatrix endTime: " + stopTime);
            System.out.println("\n[DEBUG] finalizePhenoMatrix elapsedTime: " + elapsedTime + " seconds");
            System.out.println("\n[DEBUG] Finished execution of finalizePhenoMatrix.");
        }

    }

    /**
     * Method to merge arrays
     * 
     * @param chr
     * @param caseArray
     * @param assocArray
     * @param length_entry_assoc_list
     * @param mafThresholdS
     * @param infoThresholdS
     * @param hweCohortThresholdS
     * @param hweCasesThresholdS
     * @param hweControlsThresholdS
     * @return
     * @throws IOException
     */
    private static ArrayList<String> mergeArrays(String chr, ArrayList<String> caseArray, ArrayList<String> assocArray,
            int length_entry_assoc_list, String mafThresholdS, String infoThresholdS, String hweCohortThresholdS, String hweCasesThresholdS,
            String hweControlsThresholdS) throws IOException {

        Double mafThreshold = Double.parseDouble(mafThresholdS);
        Double infoThreshold = Double.parseDouble(infoThresholdS);
        Double hweCohortThreshold = Double.parseDouble(hweCohortThresholdS);
        Double hweCasesThreshold = Double.parseDouble(hweCasesThresholdS);
        Double hweControlsThreshold = Double.parseDouble(hweCohortThresholdS);
        int real_length_assoc = 0;
        if (chr.equals("23")) {
            real_length_assoc = 67;
        } else {
            // real_length_assoc = 69;
            real_length_assoc = 67;
        }

        ArrayList<String> summaryTmp = new ArrayList<>();
        // We just need to put the information of the differente arraysList into the returned arrayList
        // First the position (as string):

        // We put the chromosome number
        summaryTmp.add(chr);

        if (caseArray != null) {
            for (int i = 0; i < caseArray.size(); i++) {
                summaryTmp.add(caseArray.get(i));
            }
        } else {
            System.out.println("\n[DEBUG] Extranisimo, este caso no deberia darse.");
        }

        if (assocArray != null) {
            int assoc_length = assocArray.size();
            // We do not store the first 4 field because they are not necessary or are repeated:
            // These four fields are:
            // alternative_ids, rsid, chromosome, position
            for (int i = 4; i < assoc_length; i++) {
                summaryTmp.add(assocArray.get(i));
            }
        } else {
            // 67 is the number of field in a mixed assocArray
            for (int i = 4; i < real_length_assoc; i++) {
                summaryTmp.add("NA");
            }
        }

        return summaryTmp;
    }

    /**
     * Method to get an allele
     * 
     * @param allele1
     * @param allele2
     * @param typeAllele
     * @return
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static String getAllele(String allele1, String allele2, String typeAllele) throws GuidanceTaskException {

        // Lets compute the reverse of allele1
        // String a1Reverse = new StringBuffer().reverse(allele1).toString();

        // Lets compute the reverse of allele2
        // String a2Reverse = new StringBuffer().reverse(allele2).toString();

        char aChar;
        // Lets compute the complement of allele1
        String a1Complement = "";
        for (int i = 0; i < allele1.length(); i++) {
            aChar = allele1.charAt(i);
            switch (allele1.charAt(i)) {
                case 'A':
                    aChar = 'T';
                    break;
                case 'T':
                    aChar = 'A';
                    break;
                case 'C':
                    aChar = 'G';
                    break;
                case 'G':
                    aChar = 'C';
                    break;
                default:
                    aChar = 'X';
            }
            a1Complement = a1Complement + aChar;
        }

        // Lets compute the complement of allele2
        String a2Complement = "";
        for (int i = 0; i < allele2.length(); i++) {
            aChar = allele2.charAt(i);
            switch (allele2.charAt(i)) {
                case 'A':
                    aChar = 'T';
                    break;
                case 'T':
                    aChar = 'A';
                    break;
                case 'C':
                    aChar = 'G';
                    break;
                case 'G':
                    aChar = 'C';
                    break;
                default:
                    aChar = 'X';
            }
            a2Complement = a2Complement + aChar;
        }

        // Lets compute the reverse of a1Complement
        // String a1ComplementAndReverse = new StringBuffer().reverse(a1Complement).toString();
        // Lets compute the reverse of a2Complement
        // String a2ComplementAndReverse = new StringBuffer().reverse(a2Complement).toString();

        if (typeAllele.equals("reverse")) {
            return allele2 + "_" + allele1;
        } else if (typeAllele.equals("complement")) {
            return a1Complement + "_" + a2Complement;
        } else if (typeAllele.equals("complementAndReverse")) {
            return a2Complement + "_" + a1Complement;
        } else {
            throw new GuidanceTaskException("Error, the option (" + typeAllele + ") for creating alternative alleles is now correct!!");
        }

    }

    /**
     * Method to create a hashtable from the header of particular files
     * 
     * @param line
     * @param separator
     * @return
     */
    private static Hashtable<String, Integer> createHashWithHeader(String line, String separator) {
        Hashtable<String, Integer> myHashLine = new Hashtable<>();

        String[] splitted = line.split(separator);
        for (int i = 0; i < splitted.length; i++) {
            myHashLine.put(splitted[i], i);
        }
        return myHashLine;
    }

    /**
     * Method to create a hashtable from the header of particular files
     * 
     * @param line
     * @param separator
     * @return
     */
    private static Hashtable<Integer, String> createHashWithHeaderReversed(String line, String separator) {
        Hashtable<Integer, String> myHashLine = new Hashtable<>();

        String[] splitted = line.split(separator);
        for (int i = 0; i < splitted.length; i++) {
            myHashLine.put(i, splitted[i]);
        }
        return myHashLine;
    }

    /**
     * Method to delete file sourceFile
     * 
     * @param sourceFile
     * @throws IOException
     * @throws InterruptedException
     * @throws Exception
     */
    public static void deleteFile(String sourceFile) throws GuidanceTaskException {

        // We cannot delete the file
        // File tmpFile = new File(sourceFile);
        // tmpFile.delete();

        try (FileOutputStream fis = new FileOutputStream(new File(sourceFile))) {
            FileChannel destination = fis.getChannel();
            String newData = "This file has been compressed. See the .gz version";
            ByteBuffer buf = ByteBuffer.allocate(64);
            buf.clear();
            buf.put(newData.getBytes());
            buf.flip();

            while (buf.hasRemaining()) {
                destination.write(buf);
            }
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
    }

    /**
     * Method to copy file source to dest
     * 
     * @param source
     * @param dest
     * @throws IOException
     */
    private static void copyFile(File source, File dest) throws IOException {
        Files.copy(source.toPath(), dest.toPath(), REPLACE_EXISTING);
    }

    /**
     * Method to zip a file
     * 
     * @param sourceFilePath
     * @param destZipFilePath
     */
    private static void gzipFile(String sourceFilePath, String destZipFilePath) {
        byte[] buffer = new byte[1024];
        try {
            FileOutputStream fileOutputStream = new FileOutputStream(destZipFilePath);
            GZIPOutputStream gzipOuputStream = new GZIPOutputStream(fileOutputStream);

            FileInputStream fileInput = new FileInputStream(sourceFilePath);

            int bytes_read;

            while ((bytes_read = fileInput.read(buffer)) > 0) {
                gzipOuputStream.write(buffer, 0, bytes_read);
            }

            fileInput.close();

            gzipOuputStream.finish();
            gzipOuputStream.close();

            // System.out.println("The file was compressed successfully!");
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
            ex.printStackTrace();
        }
    }

    /**
     * Method to unzip a file
     * 
     * @param compressedFile
     * @param decompressedFile
     */
    private static void gunzipFile(String compressedFile, String decompressedFile) {
        byte[] buffer = new byte[1024];

        try {
            FileInputStream fileIn = new FileInputStream(compressedFile);
            GZIPInputStream gZIPInputStream = new GZIPInputStream(fileIn);

            FileOutputStream fileOutputStream = new FileOutputStream(decompressedFile);

            int bytes_read;

            while ((bytes_read = gZIPInputStream.read(buffer)) > 0) {
                fileOutputStream.write(buffer, 0, bytes_read);
            }

            gZIPInputStream.close();
            fileOutputStream.close();

            // System.out.println("The file was decompressed successfully!");

        } catch (IOException ex) {
            System.err.println(ex.getMessage());
            ex.printStackTrace();
        }
    }

    /**
     * Reads the given output and error of a process and stores to the given file paths
     * 
     * @param in
     * @param outputPath
     * @param error
     * @param errorPath
     * @throws IOException
     */
    private static void readOutputAndError(InputStream in, String outputPath, InputStream error, String errorPath) throws IOException {
        try (BufferedInputStream bisInp = new BufferedInputStream(in);
                BufferedOutputStream bosInp = new BufferedOutputStream(new FileOutputStream(outputPath))) {

            byte[] b = new byte[1024];
            int read;
            while ((read = bisInp.read(b)) >= 0) {
                bosInp.write(b, 0, read);
            }
        }

        try (BufferedInputStream bisErr = new BufferedInputStream(error);
                BufferedOutputStream bosErr = new BufferedOutputStream(new FileOutputStream(errorPath))) {

            byte[] b = new byte[1024];
            int read;
            while ((read = bisErr.read(b)) >= 0) {
                bosErr.write(b, 0, read);
            }
        }
    }

}
