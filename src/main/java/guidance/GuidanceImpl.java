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

import guidance.exceptions.GuidanceTaskException;
import guidance.files.FileUtils;
import guidance.processes.ProcessUtils;
import guidance.utils.Headers;

import java.io.File;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.FileReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.TreeMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class GuidanceImpl {

	// Debug flag
	private static final boolean DEBUG = true;

	// Environment variable names
	private static final String PLINKBINARY = "PLINKBINARY";
	private static final String GTOOLBINARY = "GTOOLBINARY";
	private static final String QCTOOLBINARY = "QCTOOLBINARY";
	private static final String SHAPEITBINARY = "SHAPEITBINARY";
	private static final String EAGLEBINARY = "EAGLEBINARY";
	private static final String SAMTOOLSBINARY = "SAMTOOLSBINARY";
	private static final String IMPUTE2BINARY = "IMPUTE2BINARY";
	private static final String MINIMACBINARY = "MINIMACBINARY";
	private static final String RSCRIPTBINDIR = "RSCRIPTBINDIR";
	private static final String RSCRIPTDIR = "RSCRIPTDIR";
	private static final String SNPTESTBINARY = "SNPTESTBINARY";
	private static final String BCFTOOLSBINARY = "BCFTOOLSBINARY";

	// Method headers
	private static final String HEADER_CONVERT_FROM_BED_TO_BED = "[convertFromBedToBed]";
	private static final String HEADER_CONVERT_FROM_BED_TO_PED = "[convertFromBedToPed]";
	private static final String HEADER_CONVERT_FROM_PED_TO_GEN = "[convertFromPedToGen]";
	private static final String HEADER_CREATE_RSID_LIST = "[CreateRsIdList]";
	private static final String HEADER_WRITE_OUT_PAIRS_FILE = "[writeOutPairsFile]";
	private static final String HEADER_GTOOLS = "[gtoolS]";
	private static final String HEADER_QCTOOLS = "[qctools]";
	private static final String HEADER_PHASING = "[phasing]";
	private static final String HEADER_SAMTOOLS = "[samtoolsBgzip]";
	private static final String HEADER_FILTER_HAPLOTYPES = "[filterHaplotypes]";
	private static final String HEADER_IMPUTE = "[impute]";
	private static final String HEADER_MINIMAC = "[minimac]";
	private static final String HEADER_GENERATE_QQ_MANHATTAN_PLOTS = "[generateQQManhattanPlots]";
	private static final String HEADER_SNPTEST = "[snptest]";
	private static final String HEADER_PHENO = "[phenoAnalysis]";

	private static final String SEX1 = "males";
	private static final String SEX2 = "females";
	private static final String NO_SEX = "";

	// Commonly used characters
	private static final String NEW_LINE = "\n";
	private static final String TAB = "\t";
	private static final String SPACE = "\\s+"; // "\\s+";
	private static final String SPACE_WRITE = " "; // "\u0020";
	private static final String CHR_23 = "23";

	// File extensions
	private static final String TEMP_EXTENSION = ".temp";
	private static final String STDOUT_EXTENSION = ".stdout";
	private static final String STDERR_EXTENSION = ".stderr";

	// Debug messages
	private static final String MSG_CMD = "Command: ";
	private static final String FILE_SUFFIX = " file";

	// Error messages
	private static final String ERROR_ON_FILE = "Error: the file ";
	private static final String ERROR_SUFFIX_RENAMED_FILE = " was not succesfully renamed";
	private static final String ERROR_BINARY_EXEC = "Error executing binary job, exit value is: ";
	private static final String ERROR_FILE_CREATION = "Error: Cannot create ";

	/**
	 * Private constructor to avoid instantiation
	 */
	private GuidanceImpl() {
		// Private constructor to avoid instantiation
	}

	public static String getSex1() {
		return SEX1;
	}

	public static String getSex2() {
		return SEX2;
	}

	public static String getNoSex() {
		return NO_SEX;
	}

	private static String loadFromEnvironment(String envVarName, String methodHeader) throws GuidanceTaskException {
		String envVar = System.getenv(envVarName);
		if (envVar == null) {
			throw new GuidanceTaskException(
					methodHeader + "Error, " + envVarName + " environment variable is not defined in .bashrc!!!");
		}

		return envVar;
	}

	/**
	 * Method to split the gmapFile in case of eagle/impute combination
	 * 
	 * @param gmapFile
	 * @param gmapFileChr
	 * @param chromo
	 * @param cmdToStore
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws Exception
	 */
	public static void createSplitedFiles(String gmapFile, String gmapFileChr, String chromo, String cmdToStore)
			throws IOException, InterruptedException, Exception {

		// String gmapFileChrGz = gmapFileChr + ".gz";

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running createSplitedFiles with parameters:");
			System.out.println("[DEBUG] \t- Input gmapFile      : " + gmapFile);
			System.out.println("[DEBUG] \t- Output gmapFileChr  : " + gmapFileChr);
			System.out.println("[DEBUG] \t- Chromosome          : " + chromo);
			System.out.println("\n");
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// Input file
		FileInputStream fr = new FileInputStream(gmapFile);
		GZIPInputStream gzInput = new GZIPInputStream(fr);
		InputStreamReader unzipFile = new InputStreamReader(gzInput);
		BufferedReader br = new BufferedReader(unzipFile);

		// Output file
		FileWriter fw = new FileWriter(gmapFileChr);
		BufferedWriter writerInfo = new BufferedWriter(fw);

		// Header handling
		String line = br.readLine();

		// Header without chr column for impute2
		String[] headerSplited = line.split("\\s");
		String headerLine = "";
		String finalLine = "";
		for (int elem = 1; elem < headerSplited.length - 1; elem++)
			headerLine += headerSplited[elem] + "\t";
		headerLine += headerSplited[headerSplited.length - 1] + "\n";
		System.out.println("[DEBUG] \t- Header: " + headerLine);

		writerInfo.write(headerLine);

		while ((line = br.readLine()) != null) {

			String[] lineSplited = line.split("\\s");
			String chr = lineSplited[0];
			if (Integer.parseInt(chr) == Integer.parseInt(chromo)) {
				for (int elem = 1; elem < lineSplited.length - 1; elem++) {
					finalLine += lineSplited[elem] + "\t";
				}
				finalLine += lineSplited[lineSplited.length - 1] + "\n";
				// writerInfo.write(finalLine);
			}
		}
		writerInfo.write(finalLine);
		writerInfo.close();
		br.close();

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] createSplitedFiles startTime: " + startTime);
			System.out.println("\n[DEBUG] createSplitedFiles endTime: " + stopTime);
			System.out.println("\n[DEBUG] createSplitedFiles elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of createSplitedFiles.");
		}
	}

	/**
	 * Method to perform the conversion from Bed to Bed Format file
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
	public static void convertFromBedToBed(String bedFile, String bimFile, String famFile, String newBedFile,
			String newBimFile, String newFamFile, String logFile, String chromo, String cmdToStore)
			throws GuidanceTaskException {

		String plinkBinary = loadFromEnvironment(PLINKBINARY, HEADER_CONVERT_FROM_BED_TO_BED);

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
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		String basePath = newBedFile.substring(0, newBedFile.length() - 4);

		String cmd = plinkBinary + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile + " --chr " + chromo
				+ " --out " + basePath + " --make-bed";

		if (DEBUG) {
			System.out.println(HEADER_CONVERT_FROM_BED_TO_BED + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, basePath + STDOUT_EXTENSION, basePath + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}
		System.out.println("The exit value is " + exitValue);

		// Check process exit value
		if (exitValue != 0) {
			throw new GuidanceTaskException(HEADER_CONVERT_FROM_BED_TO_BED + ERROR_BINARY_EXEC + exitValue);
		}

		// If there is not output in the convertFromBedToBed process. Then we have to
		// create some empty outputs
		try {
			FileUtils.createEmptyFile(newBedFile, HEADER_CONVERT_FROM_BED_TO_BED);
			FileUtils.createEmptyFile(newBimFile, HEADER_CONVERT_FROM_BED_TO_BED);
			FileUtils.createEmptyFile(newFamFile, HEADER_CONVERT_FROM_BED_TO_BED);
			FileUtils.createEmptyFile(logFile, HEADER_CONVERT_FROM_BED_TO_BED);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] convertFromBedToBed startTime: " + startTime);
			System.out.println("\n[DEBUG] convertFromBedToBed endTime: " + stopTime);
			System.out.println("\n[DEBUG] convertFromBedToBed elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of convertFromBedToBed.");
		}
	}

	/**
	 * Method to perform split chr 23 in males and females
	 *
	 * @param myPrefix
	 * @param bedFile
	 * @param bimFile
	 * @param famFile
	 * @param bedChr23File
	 * @param bimChr23File
	 * @param famChr23File
	 * @param logFile
	 * @param sex
	 * @param theChromo
	 * @param cmdToStore
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws Exception
	 */
	public static void splitChr23(String bedFile, String bimFile, String famFile, String bedChr23File,
			String bimChr23File, String famChr23File, String logFile, String sex, String theChromo, String cmdToStore)
			throws IOException, InterruptedException, Exception {

		String myPrefix = bedChr23File.substring(0, bedChr23File.length() - 4);

		String plinkBinary = System.getenv("PLINKBINARY");
		if (plinkBinary == null) {
			throw new Exception("[splitChr23] Error, PLINKBINARY environment variable is not defined in .bashrc!!!");
		}

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running splitChr23 with parameters:");
			System.out.println("[DEBUG] \t- Input prefix    : " + myPrefix);
			System.out.println("[DEBUG] \t- Input bedFile   : " + bedFile);
			System.out.println("[DEBUG] \t- Input bimFile   : " + bimFile);
			System.out.println("[DEBUG] \t- Input famFile   : " + famFile);
			System.out.println("[DEBUG] \t- Output newBedFile  : " + bedChr23File);
			System.out.println("[DEBUG] \t- Output newBimFile  : " + bimChr23File);
			System.out.println("[DEBUG] \t- Output newFamFile  : " + famChr23File);
			System.out.println("[DEBUG] \t- Output logFile  : " + logFile);
			System.out.println("[DEBUG] \t- Gender          : " + sex);
			System.out.println("[DEBUG] \t- Chromosome      : " + theChromo);
			System.out.println("\n");
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();
		String cmd = null;

		if (sex.equals("males")) {
			cmd = plinkBinary + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
					+ " --filter-males --out " + myPrefix + " --make-bed";

		} else if (sex.equals("females")) {
			cmd = plinkBinary + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile
					+ " --filter-females --out " + myPrefix + " --make-bed";

		}

		if (DEBUG) {
			System.out.println("[splitChr23] Command: " + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, bedFile + STDOUT_EXTENSION, bedFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		// Check process exit value
		if (exitValue == 12) {
			System.out
					.println("Look at the log file. Error: All people removed due to gender filter (--filter-males).");
			try {
				FileUtils.createEmptyFile(bedChr23File, HEADER_CONVERT_FROM_BED_TO_BED);
				FileUtils.createEmptyFile(bimChr23File, HEADER_CONVERT_FROM_BED_TO_BED);
				FileUtils.createEmptyFile(famChr23File, HEADER_CONVERT_FROM_BED_TO_BED);
				FileUtils.createEmptyFile(logFile, HEADER_CONVERT_FROM_BED_TO_BED);
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}
		} else if (exitValue != 0) {
			throw new GuidanceTaskException(HEADER_CONVERT_FROM_BED_TO_BED + ERROR_BINARY_EXEC + exitValue);
		}

		/*
		 * FileUtils.move(myPrefix + ".bed", bedChr23File); FileUtils.move(myPrefix +
		 * ".bim", bimChr23File); FileUtils.move(myPrefix + ".fam", famChr23File);
		 * FileUtils.move(myPrefix + ".log", logFile);
		 */

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] splitChr23 startTime: " + startTime);
			System.out.println("\n[DEBUG] splitChr23 endTime: " + stopTime);
			System.out.println("\n[DEBUG] splitChr23 elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of splitChr23.");
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
	public static void convertFromBedToPed(String bedPrefix, String bedFile, String bimFile, String famFile,
			String pedFile, String mapFile, String logFile, String chromo, String cmdToStore)
			throws GuidanceTaskException {

		String plinkBinary = loadFromEnvironment(PLINKBINARY, HEADER_CONVERT_FROM_BED_TO_PED);

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
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		String cmd = plinkBinary + " --noweb --bfile " + bedPrefix + " --chr " + chromo + " --recode --out " + pedFile;

		if (DEBUG) {
			System.out.println(HEADER_CONVERT_FROM_BED_TO_PED + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, pedFile + STDOUT_EXTENSION, pedFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			throw new GuidanceTaskException(HEADER_CONVERT_FROM_BED_TO_PED + ERROR_BINARY_EXEC + exitValue);
		}

		// Rename file (or directory)
		boolean success = FileUtils.move(pedFile + ".ped", pedFile);
		if (!success) {
			throw new GuidanceTaskException(
					HEADER_CONVERT_FROM_BED_TO_PED + ERROR_ON_FILE + pedFile + ERROR_SUFFIX_RENAMED_FILE);
		}
		// Rename file (or directory)
		success = FileUtils.move(pedFile + ".map", mapFile);
		if (!success) {
			throw new GuidanceTaskException(
					HEADER_CONVERT_FROM_BED_TO_PED + ERROR_ON_FILE + mapFile + ERROR_SUFFIX_RENAMED_FILE);
		}
		// Rename file (or directory)
		success = FileUtils.move(pedFile + ".log", logFile);
		if (!success) {
			throw new GuidanceTaskException(
					HEADER_CONVERT_FROM_BED_TO_PED + ERROR_ON_FILE + logFile + ERROR_SUFFIX_RENAMED_FILE);
		}

		// If there is not output in the convertFromBedToPed process. Then we have to
		// create some empty outputs.
		try {
			FileUtils.createEmptyFile(mapFile, HEADER_CONVERT_FROM_BED_TO_PED);
			FileUtils.createEmptyFile(pedFile, HEADER_CONVERT_FROM_BED_TO_PED);
			FileUtils.createEmptyFile(logFile, HEADER_CONVERT_FROM_BED_TO_PED);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void convertFromPedToGen(String pedFile, String mapFile, String genFile, String sampleFile,
			String logFile, String cmdToStore) throws GuidanceTaskException {

		String gtoolBinary = loadFromEnvironment(GTOOLBINARY, HEADER_CONVERT_FROM_PED_TO_GEN);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running convertFromPedToGen with parameters:");
			System.out.println("[DEBUG] \t- Input pedFile      : " + pedFile);
			System.out.println("[DEBUG] \t- Input mapFile      : " + mapFile);
			System.out.println("[DEBUG] \t- Output genFile     : " + genFile);
			System.out.println("[DEBUG] \t- Output sampleFile  : " + sampleFile);
			System.out.println("[DEBUG] \t- Output logFile     : " + logFile);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		String cmd = gtoolBinary + " -P --ped " + pedFile + " --map " + mapFile + " --og " + genFile + " --os "
				+ sampleFile + " --binary_phenotype --order --log " + logFile;

		if (DEBUG) {
			System.out.println(HEADER_CONVERT_FROM_PED_TO_GEN + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, genFile + STDOUT_EXTENSION, genFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			throw new GuidanceTaskException(HEADER_CONVERT_FROM_PED_TO_GEN + ERROR_BINARY_EXEC + exitValue);
		}

		// If there is not output in the convertFromPedToGen process. Then we have to
		// create some empty outputs.
		try {
			FileUtils.createEmptyFile(genFile, HEADER_CONVERT_FROM_PED_TO_GEN);
			FileUtils.createEmptyFile(sampleFile, HEADER_CONVERT_FROM_PED_TO_GEN);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Now we have to change the missing/unaffection/affection coding in the sample
		// file from -9/1/2 to NA/0/1
		// Then, we go to the sample file and change the last column.
		String sampleFileChangedName = sampleFile + ".changed";
		File changedSampleFile = new File(sampleFileChangedName);
		try {
			if (!changedSampleFile.createNewFile()) {
				throw new IOException(ERROR_FILE_CREATION + changedSampleFile + FILE_SUFFIX);
			}
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// We read each line of the genFile and look for the rsID into the newPosMap
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
					System.out.println(HEADER_CONVERT_FROM_PED_TO_GEN
							+ "Error changing the sample file. Invalid affection coding in line " + myNewLine);
					// throw new Exception("Error changing the sample file. Invalid affection coding
					// in line " +
					// counter);
				}
				writerPos.write(myNewLine);
				writerPos.newLine();
			}

			writerPos.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Rename file (or directory)
		boolean success = FileUtils.move(sampleFileChangedName, sampleFile);
		if (!success) {
			throw new GuidanceTaskException(
					HEADER_CONVERT_FROM_PED_TO_GEN + ERROR_ON_FILE + sampleFile + ERROR_SUFFIX_RENAMED_FILE);
		}

		try {
			FileUtils.createEmptyFile(logFile, HEADER_CONVERT_FROM_PED_TO_GEN);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void createRsIdList(String genOrBimFile, String exclCgatFlag, String pairsFile, String inputFormat,
			String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running :q with parameters:");
			System.out.println("[DEBUG] \t- Input genOrBimFile : " + genOrBimFile);
			System.out.println("[DEBUG] \t- Input exclCgatFlag : " + exclCgatFlag);
			System.out.println("[DEBUG] \t- Output pairsFile   : " + pairsFile);
			System.out.println("[DEBUG] \t- InputFormat        : " + inputFormat);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// Create output file
		File outPairsFile = new File(pairsFile);
		try {
			if (!outPairsFile.createNewFile()) {
				throw new IOException(HEADER_CREATE_RSID_LIST + ERROR_FILE_CREATION + outPairsFile + FILE_SUFFIX);
			}
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check if the file is gzip or not
		boolean thisIsGz = genOrBimFile.endsWith(".gz");
		// Process the file depending on its format
		if (thisIsGz) {
			System.out.println(HEADER_CREATE_RSID_LIST + "It seems the file " + genOrBimFile + " is a gzip file");

			try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(genOrBimFile));
					InputStreamReader decoder = new InputStreamReader(inputGz);
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
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] createRsIdList StartTime: " + startTime);
			System.out.println("\n[DEBUG] createRsIdList endTime: " + stopTime);
			System.out.println("\n[DEBUG] createRsIdList elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of createRsIdList.");
		}

	}

	private static void writeOutPairsFile(BufferedReader br, File outPairsFile, String inputFormat, String exclCgatFlag,
			String pairsFile) throws GuidanceTaskException {

		try (BufferedWriter writerPairs = new BufferedWriter(new FileWriter(outPairsFile))) {
			if (inputFormat.equals("BED")) {
				String line = null;
				while ((line = br.readLine()) != null) {
					String[] splittedLine = line.split(SPACE);// delimiter I assume single space.
					String allele = splittedLine[4] + splittedLine[5]; // store Allele (AC,GC,AG, GT,..., etc.)

					// Store rsID of the SNP which its allele is AT or TA or GC or CG into the
					// .pairs file
					if (exclCgatFlag.equals("YES")) {
						// Then we have to see if allele is AT TA GC CG to put the rsID into the .pairs
						// file.
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

					// Store rsID of the SNP which its allele is AT or TA or GC or CG into the
					// .pairs file
					if (exclCgatFlag.equals("YES")) {
						// Then we have to see if allele is AT TA GC CG to put the rsID into the .pairs
						// file.
						if (allele.equals("AT") || allele.equals("TA") || allele.equals("GC") || allele.equals("CG")) {
							writerPairs.write(splittedLine[1]);
							writerPairs.newLine();
						}
					}
				}
			} else {
				throw new GuidanceTaskException(HEADER_WRITE_OUT_PAIRS_FILE + "Error, It was not possible to generate "
						+ pairsFile + ". The " + inputFormat + " is not valid!!");
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
	public static void gtoolS(String newGenFile, String modSampleFile, String gtoolGenFile, String gtoolSampleFile,
			String sampleExclFile, String snpWtccFile, String gtoolLogFile, String cmdToStore)
			throws GuidanceTaskException {

		String gtoolBinary = loadFromEnvironment(GTOOLBINARY, HEADER_GTOOLS);

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
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		String cmd = gtoolBinary + " -S --g " + newGenFile + " --s " + modSampleFile + " --og " + gtoolGenFile
				+ " --os " + gtoolSampleFile + " --sample_excl " + sampleExclFile + " --exclusion " + snpWtccFile
				+ " --log " + gtoolLogFile;

		if (DEBUG) {
			System.out.println(HEADER_GTOOLS + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, gtoolGenFile + STDOUT_EXTENSION, gtoolGenFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			throw new GuidanceTaskException(HEADER_GTOOLS + ERROR_BINARY_EXEC + exitValue);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] gtoolS startTime: " + startTime);
			System.out.println("\n[DEBUG] gtoolS endTime: " + stopTime);
			System.out.println("\n[DEBUG] gtoolS elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of gtoolS.");
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
	public static void qctoolS(String imputationTool, String imputeFile, String inclusionRsIdFile, String mafThresholdS,
			String filteredFile, String filteredLogFile, String cmdToStore) throws GuidanceTaskException {

		String qctoolBinary = null;
		if (imputationTool.equals("impute")) {
			qctoolBinary = loadFromEnvironment(QCTOOLBINARY, HEADER_QCTOOLS);
		} else if (imputationTool.equals("minimac")) {
			qctoolBinary = loadFromEnvironment(QCTOOLBINARY, HEADER_QCTOOLS);
			// qctoolBinary = loadFromEnvironment(QCTOOLBINARY, HEADER_QCTOOLS);
		}

		if (DEBUG) {
			System.out.println("\nRunning qctoolS for generation a subset of rsids with parameters:");
			System.out.println("\t- qctoolBinary               : " + qctoolBinary);
			System.out.println("\t- Input imputeFile           : " + imputeFile);
			System.out.println("\t- Input inclusionRsIdFile    : " + inclusionRsIdFile);
			System.out.println("\t- Input mafThreshold         : " + mafThresholdS);
			System.out.println("\t- Output filteredFile        : " + filteredFile);
			System.out.println("\t- Output filteredLogFile     : " + filteredLogFile);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}
		long startTime = System.currentTimeMillis();

		// We have to make sure whether we are using renamed files of the original gz
		// files.
		// We detect this situation by scanning the last three characters:
		String extension = imputeFile.substring(Math.max(0, imputeFile.length() - 3));

		String imputeFileGz = null;
		if (extension.equals(".gz")) {
			imputeFileGz = imputeFile;
		} else {
			// If imputeFile exists, then imputeFileGz exists also.
			// We reused the imputFileGz
			imputeFileGz = imputeFile + ".gz";
			FileUtils.gzipFile(imputeFile, imputeFileGz);

		}

		String cmd = null;

		boolean execution = false;

		if (imputationTool.equals("impute")) {
			cmd = qctoolBinary + " -g " + imputeFileGz + " -og " + filteredFile + " -incl-rsids " + inclusionRsIdFile
					+ " -omit-chromosome -force -log " + filteredLogFile + " -maf " + mafThresholdS + " 1";

			execution = true;

		} else if (imputationTool.equals("minimac")) {
			cmd = qctoolBinary + " -g " + imputeFileGz + " -og " + filteredFile + " -incl-rsids " + inclusionRsIdFile
					+ " -omit-chromosome -force -log " + filteredLogFile + " -vcf-genotype-field GP";

			File file = new File(imputeFileGz);
			if (file.length() <= 28) {
				System.out.println("\n[DEBUG] Execution: FALSE");
				execution = false;
			} else {
				System.out.println("\n[DEBUG] Execution: TRUE");
				execution = true;
			}

		} else {
			System.out.println("[qctoolS]: The Imputation Tool is not valid. Only IMPUTE2 and MINIMAC3 are available");
		}

		if (DEBUG) {
			System.out.println("\n[DEBUG] Command: " + cmd);
		}

		if (execution) {
			// Execute the command retrieving its exitValue, output and error
			int exitValue = -1;
			try {
				String outputBase = filteredFile.substring(0, filteredFile.length() - 3);
				// exitValue = ProcessUtils.execute(cmd, outputBase + STDOUT_EXTENSION,
				// outputBase + STDERR_EXTENSION,
				// "LD_LIBRARY_PATH");
				exitValue = ProcessUtils.execute(cmd, outputBase + STDOUT_EXTENSION, outputBase + STDERR_EXTENSION);
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}

			if (exitValue != 0) {
				throw new GuidanceTaskException(HEADER_QCTOOLS + ERROR_BINARY_EXEC + exitValue);
			}
		} else {
			System.out.println("[qctoolS]: Empty chunk. Output files being created...");
			try {
				FileUtils.delete(filteredFile);
			} catch (GuidanceTaskException gte) {
				System.out.println(gte);
			}
			try {
				FileUtils.delete(filteredLogFile);
			} catch (GuidanceTaskException gte) {
				System.out.println(gte);
			}
			try {
				String plainOutputFilteredFile = filteredFile.substring(0, filteredFile.length() - 3);
				FileUtils.createEmptyFile(plainOutputFilteredFile, HEADER_GTOOLS);
				FileUtils.gzipFile(plainOutputFilteredFile, filteredFile);
				FileUtils.delete(plainOutputFilteredFile);
				FileUtils.createEmptyFile(filteredLogFile, HEADER_GTOOLS);
			} catch (IOException gte) {
				throw new GuidanceTaskException(HEADER_QCTOOLS + ERROR_BINARY_EXEC + gte);
			}

		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void createListOfExcludedSnps(String phasingHapsFile, String excludedSnpsFile, String exclCgatFlag,
			String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running createListOfExcludedSnps method:");
			System.out.println("[DEBUG] \t- Input shapeitHapsFile   : " + phasingHapsFile);
			System.out.println("[DEBUG] \t- Output excludedSnpsFile : " + excludedSnpsFile);
			System.out.println("[DEBUG] \t- Input exclCgatFlag      : " + exclCgatFlag);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// We have to check the haplotypes file and extract tha allele information. For
		// the format file used by shapeit,
		// the indices for alleles are:

		// int rsIdIndex = 1;
		int posIndex = 2;
		int a1Index = 3;
		int a2Index = 4;
		// And the column separator is:
		String separator = " ";

		// We have to make sure whether we are using renamed files of the original gz
		// files.
		// We detect this situation by scanning the last three characters of the file
		// name:
		String extension = phasingHapsFile.substring(Math.max(0, phasingHapsFile.length() - 3));

		String phasingHapsFileGz = null;
		if (extension.equals(".gz")) {
			phasingHapsFileGz = phasingHapsFile;
		} else {
			// If the renamed shapeitHapsFile exists, then an extended .gz version exists
			// also.
			phasingHapsFileGz = phasingHapsFile + ".gz";
		}

		// Array of string to store positions of SNPs to exclude
		ArrayList<String> excludeList = new ArrayList<>();

		// Then, we read the gz File line by line
		try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(phasingHapsFileGz));
				InputStreamReader decoder = new InputStreamReader(inputGz);
				BufferedReader br = new BufferedReader(decoder)) {

			String line = "";
			while ((line = br.readLine()) != null) {
				String[] splitted = line.split(separator);// delimiter defined previously.
				// String allele1 = splitted[a1Index];
				// String allele2 = splitted[a2Index];
				String positionS = splitted[posIndex];

				String allele = splitted[a1Index] + splitted[a2Index];
				if (exclCgatFlag.equals("YES")) {
					// Then we have to see if allele is AT TA GC CG to exclude it.
					if (allele.equals("AT") || allele.equals("TA") || allele.equals("GC") || allele.equals("CG")) {
						excludeList.add(positionS);
					}
				}
			}
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Finally we put the excludedList into the outputFile
		// We have to create the outputFile:
		// We verify that a file with the same name does not exist. (It should not
		// exist!!)
		File outputFile = new File(excludedSnpsFile);
		try {
			if (!outputFile.createNewFile()) {
				throw new IOException(ERROR_FILE_CREATION + outputFile + FILE_SUFFIX);
			}
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile))) {
			for (int i = 0; i < excludeList.size(); i++) {
				writer.write(excludeList.get(i) + NEW_LINE);
			}
			writer.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void phasingBed(String chromo, String bedFile, String bimFile, String famFile, String gmapFile,
			String phasingHapsFile, String phasingSampleFile, String phasingLogFile, String sex, String phasingTool,
			String cmdToStore) throws GuidanceTaskException, IOException {

		String phasingBinary = "";

		if (phasingTool.equals("shapeit")) {
			phasingBinary = loadFromEnvironment(SHAPEITBINARY, HEADER_PHASING);
		} else if (phasingTool.equals("eagle")) {
			phasingBinary = loadFromEnvironment(EAGLEBINARY, HEADER_PHASING);
		} else {
			System.out.println("\n[DEBUG] Only Shapeit and Eagle are available for phasing");
		}

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running phasing with parameters:");
			System.out.println("[DEBUG] \t- phasingBinary            : " + phasingBinary);
			System.out.println("[DEBUG] \t- Input bedFile            : " + bedFile);
			System.out.println("[DEBUG] \t- Input bimFile            : " + bimFile);
			System.out.println("[DEBUG] \t- Input famFile            : " + famFile);
			System.out.println("[DEBUG] \t- Input gmapFile           : " + gmapFile);
			System.out.println("[DEBUG] \t- Output phasingHapsFile   : " + phasingHapsFile);
			System.out.println("[DEBUG] \t- Output phasingSampleFile : " + phasingSampleFile);
			System.out.println("[DEBUG] \t- Output phasingLogFile    : " + phasingLogFile);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		String cmd = null;

		boolean execute = true;

		String myPrefix = phasingSampleFile.substring(0, phasingSampleFile.length() - 7);
		System.out.println("myPrefix on phasingBed is: " + myPrefix);

		if (phasingTool.equals("shapeit")) {
			if (chromo.equals("23")) {
				cmd = phasingBinary + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map "
						+ gmapFile + " --chrX --output-max " + phasingHapsFile + " " + phasingSampleFile
						+ " --thread 47 --effective-size 20000 --output-log " + phasingLogFile;
			} else {
				cmd = phasingBinary + " --input-bed " + bedFile + " " + bimFile + " " + famFile + " --input-map "
						+ gmapFile + " --output-max " + phasingHapsFile + " " + phasingSampleFile
						+ " --thread 47 --effective-size 20000 --output-log " + phasingLogFile;
			}
		} else if (phasingTool.equals("eagle")) {

			if (chromo.equals("23")) {

				if (sex.equals(SEX1)) {

					String plinkBinary = loadFromEnvironment(PLINKBINARY, HEADER_PHASING);
					String baseDirOrigin = bedFile.substring(0, bedFile.length() - 4);
					String baseDirDest = phasingSampleFile.substring(0, phasingSampleFile.length() - 7);
					cmd = plinkBinary + " --bfile " + baseDirOrigin + " --recode vcf --out " + baseDirDest;

					int exitValue = -1;
					try {
						exitValue = ProcessUtils.execute(cmd, baseDirOrigin + STDOUT_EXTENSION,
								baseDirOrigin + STDERR_EXTENSION);
					} catch (IOException ioe) {
						throw new GuidanceTaskException(ioe);
					}

					// Check process exit value
					if (exitValue != 0) {
						System.err.println(
								HEADER_PHASING + "Warning executing phasingProc job, exit value is: " + exitValue);
						System.err.println(HEADER_PHASING + "                         (This warning is not fatal).");
					}

					if (DEBUG) {
						System.out.println(HEADER_PHASING + MSG_CMD + cmd);
					}

					// https://www.cog-genomics.org/plink/1.9/data
					String bcfBinary = loadFromEnvironment(BCFTOOLSBINARY, HEADER_PHASING);
					cmd = bcfBinary + " convert " + baseDirDest + ".vcf" + " --hapsample " + baseDirDest + " --vcf-ids";
					System.out.println("bcf call for males: " + cmd);

					exitValue = -1;
					try {
						exitValue = ProcessUtils.execute(cmd, baseDirOrigin + STDOUT_EXTENSION,
								baseDirOrigin + STDERR_EXTENSION);
					} catch (IOException ioe) {
						throw new GuidanceTaskException(ioe);
					}

					// Check process exit value
					if (exitValue != 0) {
						System.err.println(
								HEADER_PHASING + "Warning executing phasingProc job, exit value is: " + exitValue);
						System.err.println(HEADER_PHASING + "                         (This warning is not fatal).");
					}

					if (DEBUG) {
						System.out.println(HEADER_PHASING + MSG_CMD + cmd);
					}

					String generatedSample = phasingSampleFile + "s";

					FileUtils.move(phasingSampleFile, generatedSample);

					cmd = "echo \"ID_1 ID_2 missing\" > " + phasingSampleFile + "; echo \"0 0 0\" >> "
							+ phasingSampleFile + "; tail -n +3 " + generatedSample
							+ " | tr \"_\" \" \" | awk '{ print $1\"_\"$2\" \"$3\" 0\" }' >> " + phasingSampleFile;

					exitValue = -1;

					try {
						exitValue = ProcessUtils.executeWithoutOutputs(cmd);
					} catch (IOException ioe) {
						throw new GuidanceTaskException(ioe);
					}

					// Check process exit value
					if (exitValue != 0) {
						System.err.println(
								HEADER_PHASING + "Warning executing phasingProc job, exit value is: " + exitValue);
						System.err.println(HEADER_PHASING + "                         (This warning is not fatal).");
					}

					if (DEBUG) {
						System.out.println(HEADER_PHASING + MSG_CMD + cmd);
					}

					execute = false;

				} else if (sex.equals(SEX2)) {
					cmd = phasingBinary + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile + " --chrom "
							+ chromo + " --geneticMapFile " + gmapFile + " --numThreads 47 --outPrefix " + myPrefix;
				}

			} else {
				cmd = phasingBinary + " --bed " + bedFile + " --bim " + bimFile + " --fam " + famFile + " --chrom "
						+ chromo + " --geneticMapFile " + gmapFile + " --numThreads 47 --outPrefix " + myPrefix;
			}
		}

		if (DEBUG) {
			System.out.println(HEADER_PHASING + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		// Ugly issue: If we run shapeit_v1, all the stdXXX is done stderr, and there is
		// not stdout
		// Ugly issue: If we run shapeit_v2, all the stdXXX is done stdout, and there is
		// not stderr
		if (execute) {
			int exitValue = -1;
			try {
				exitValue = ProcessUtils.execute(cmd, myPrefix + STDOUT_EXTENSION, myPrefix + STDERR_EXTENSION);
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}

			// Check process exit value
			if (exitValue != 0) {
				System.err.println(HEADER_PHASING + "Warning executing phasingProc job, exit value is: " + exitValue);
				System.err.println(HEADER_PHASING + "                         (This warning is not fatal).");
			}
		}

		// Now we rename shapeitHapsFileGz to shapeitHapsFile
		/*
		 * boolean success = FileUtils.move(phasingHapsFileGz, phasingHapsFile); if
		 * (!success) { throw new GuidanceTaskException(HEADER_PHASING_BED +
		 * ERROR_ON_FILE + phasingHapsFileGz); }
		 */

		if (phasingTool.equals("eagle")) {
			if (new File(myPrefix + ".hap.gz").exists()) {
				FileUtils.move(myPrefix + ".hap.gz", phasingHapsFile);
			} else if (new File(myPrefix + ".haps.gz").exists()) {
				FileUtils.move(myPrefix + ".haps.gz", phasingHapsFile);
			} else {
				throw new GuidanceTaskException("File " + phasingHapsFile
						+ " has not been successfully generated by phasingBed with command " + cmd);
			}
			if (new File(myPrefix + ".sample").exists()) {
				FileUtils.move(myPrefix + ".sample", phasingSampleFile);
			} else if (new File(myPrefix + ".samples").exists()) {
				FileUtils.move(myPrefix + ".samples", phasingSampleFile);
			} else {
				FileUtils.createEmptyFile(phasingSampleFile, HEADER_PHASING);
			}

			if (new File(myPrefix + ".log").exists()) {
				FileUtils.move(myPrefix + ".log", phasingLogFile);
			} else {
				FileUtils.createEmptyFile(phasingLogFile, HEADER_PHASING);
			}
		} else if (phasingTool.equals("shapeit")) {
			// Ugly, because shapeit_v2 automatically puts the .log to the file.
			// If there is not output in the impute process. Then we have to create some
			// empty outputs.
			String tmpFile = phasingLogFile + ".log";
			if (new File(tmpFile).exists()) {
				boolean success = FileUtils.move(tmpFile, phasingLogFile);
				if (!success) {
					throw new GuidanceTaskException(
							HEADER_PHASING + ERROR_ON_FILE + tmpFile + ERROR_SUFFIX_RENAMED_FILE);
				}
			} else {
				FileUtils.createEmptyFile(phasingLogFile, HEADER_PHASING);
			}
			if (!new File(phasingLogFile).exists()) {
				FileUtils.createEmptyFile(phasingLogFile, HEADER_PHASING);
			}
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] phasing startTime: " + startTime);
			System.out.println("\n[DEBUG] phasing endTime: " + stopTime);
			System.out.println("\n[DEBUG] phasing elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of phasing.");
		}
	}

	/**
	 * Method to execute phasing where input files are in GEN format
	 */
	public static void phasing(String chromo, String inputGenFile, String inputSampleFile, String gmapFile,
			String phasingHapsFile, String phasingSampleFile, String phasingLogFile, String phasingTool,
			String cmdToStore) throws IOException, InterruptedException, GuidanceTaskException {

		String phasingBinary = "";
		if (phasingTool.equals("shapeit")) {
			phasingBinary = loadFromEnvironment(SHAPEITBINARY, HEADER_PHASING);
		} else if (phasingTool.equals("eagle")) {
			phasingBinary = loadFromEnvironment(EAGLEBINARY, HEADER_PHASING);
		} else {
			throw new GuidanceTaskException(
					HEADER_PHASING + ERROR_ON_FILE + " Only Shapeit and Eagle are available for phasing");
		}

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running phasing with parameters:");
			System.out.println("[DEBUG] \t- phasingBinary            : " + phasingBinary);
			System.out.println("[DEBUG] \t- Input inputGenFile       : " + inputGenFile);
			System.out.println("[DEBUG] \t- Input inputSampleFile    : " + inputSampleFile);
			System.out.println("[DEBUG] \t- Input gmapFile           : " + gmapFile);
			System.out.println("[DEBUG] \t- Output phasingHapsFile   : " + phasingHapsFile);
			System.out.println("[DEBUG] \t- Output phasingSampleFile : " + phasingSampleFile);
			System.out.println("[DEBUG] \t- Output phasingLogFile    : " + phasingLogFile);
			System.out.println("\n");
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);

		}
		long startTime = System.currentTimeMillis();

		String cmd = null;

		if (phasingTool.equals("shapeit")) {
			if (chromo.equals("23")) {
				// If we process chromoso X (23) then we change the cmdToStore
				cmd = phasingBinary + " --input-gen " + inputGenFile + " " + inputSampleFile + " --input-map "
						+ gmapFile + " --chrX --output-max " + phasingHapsFile + " " + phasingSampleFile
						+ " --thread 48 --effective-size 20000 --output-log " + phasingLogFile;
			} else {
				cmd = phasingBinary + " --input-gen " + inputGenFile + " " + inputSampleFile + " --input-map "
						+ gmapFile + " --output-max " + phasingHapsFile + " " + phasingSampleFile
						+ " --thread 48 --effective-size 20000 --output-log " + phasingLogFile;
			}
		}

		if (DEBUG) {
			System.out.println(HEADER_PHASING + MSG_CMD + cmd);
		}

		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, phasingHapsFile + STDOUT_EXTENSION,
					phasingHapsFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}
		if (exitValue != 0) {
			System.err.println("[phasing] Warning executing phasingProc job, exit value is: " + exitValue);
			System.err.println("[phasing]                         (This warning is not fatal).");
		}

		// Ugly, because phasing_v2 automatically puts the .log to the file.
		// If there is not output in the impute process. Then we have to create some
		// empty outputs.
		String tmpFile = phasingLogFile + ".log";
		File tmpPhasingLogFile = new File(tmpFile);
		if (tmpPhasingLogFile.exists()) {
			tmpPhasingLogFile.renameTo(new File(phasingLogFile));
		}

		/*
		 * // Now we rename phasingHapsFileGz to phasingHapsFile File source = new
		 * File(phasingHapsFile); File dest = new File(phasingHapsFile);
		 * copyFile(source, dest);
		 */

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] phasing startTime: " + startTime);
			System.out.println("\n[DEBUG] phasing endTime: " + stopTime);
			System.out.println("\n[DEBUG] phasing elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of phasing.");
		}
	}

	public static void newSample(String sampleFile, String phasingSampleFile, String phasingNewSampleFile,
			String responseVar, String covariables, String cmdToStore) throws IOException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running newSample with parameters:");
			System.out.println("[DEBUG] \t- Input sampleFile            : " + sampleFile);
			System.out.println("[DEBUG] \t- Input phasingSampleFile     : " + phasingSampleFile);
			System.out.println("[DEBUG] \t- Output phasingNewSampleFile : " + phasingNewSampleFile);
			System.out.println("[DEBUG] \t- Input responseVar           : " + responseVar);
			System.out.println("[DEBUG] \t- Input covariables           : " + covariables);
			System.out.println("\n");
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
		}

		long startTime = System.currentTimeMillis();

		Hashtable<String, ArrayList<String>> loadSampleFile = new Hashtable<String, ArrayList<String>>();
		ArrayList<Integer> columnsHeaderResponse = new ArrayList<Integer>();
		ArrayList<Integer> columnsHeaderCovar = new ArrayList<Integer>();

		ArrayList<String> valuesJoined = new ArrayList<String>();

		// Read File SampleFile
		BufferedReader br = new BufferedReader(new FileReader(sampleFile));
		String line = br.readLine(); // reading the header

		String[] splitHeaderSF = line.split(SPACE);
		String[] namesHeaderCovar = covariables.split(","); // Get the pheno columns
		String[] namesHeaderResponse = responseVar.split(","); // Get the pheno columns

		for (int i = 0; i < splitHeaderSF.length; i++) {

			for (int j = 0; j < namesHeaderCovar.length; j++) {
				if (splitHeaderSF[i].equals(namesHeaderCovar[j])) {
					System.out.println("Adding column " + i + " corresponding to " + splitHeaderSF[i] + " to covars");
					columnsHeaderCovar.add(i);
				}
			}

			for (int j = 0; j < namesHeaderResponse.length; j++) {
				if (splitHeaderSF[i].equals(namesHeaderResponse[j])) {
					System.out.println("Adding column " + i + " corresponding to " + splitHeaderSF[i] + " to reponse");
					columnsHeaderResponse.add(i);
				}
			}

		}

		while ((line = br.readLine()) != null) {

			String[] lineSplited = line.split(SPACE);
			String key = lineSplited[1];

			String valueCovar = lineSplited[columnsHeaderCovar.get(0)];

			for (int i = 1; i < columnsHeaderCovar.size(); i++) {
				valueCovar += SPACE_WRITE + lineSplited[columnsHeaderCovar.get(i)];
			}

			String valueResponse = lineSplited[columnsHeaderResponse.get(0)];

			for (int i = 1; i < columnsHeaderResponse.size(); i++) {
				valueResponse += SPACE_WRITE + lineSplited[columnsHeaderResponse.get(i)];
			}

			valuesJoined = new ArrayList<String>();
			valuesJoined.add(valueCovar); // Pos 1
			valuesJoined.add(valueResponse); // Pos 2

			loadSampleFile.put(key, valuesJoined);
		}
		br.close();

		// Prepearing header output file
		String outputFile = "";
		for (int i = 0; i < 3; i++) {
			outputFile += splitHeaderSF[i] + SPACE_WRITE;
		}

		outputFile += splitHeaderSF[columnsHeaderCovar.get(0)];
		for (int i = 1; i < columnsHeaderCovar.size(); i++) {
			outputFile += SPACE_WRITE + splitHeaderSF[columnsHeaderCovar.get(i)];
		}

		outputFile += SPACE_WRITE + splitHeaderSF[columnsHeaderResponse.get(0)];
		for (int i = 1; i < columnsHeaderResponse.size(); i++) {
			outputFile += SPACE_WRITE + splitHeaderSF[columnsHeaderResponse.get(i)];
		}

		outputFile += "\n";

		// Read File phasingSampleFil
		br = new BufferedReader(new FileReader(phasingSampleFile));
		line = br.readLine(); // reading the header

		while ((line = br.readLine()) != null) {

			String[] lineSplited = line.split(SPACE);
			String key = lineSplited[1];

			valuesJoined = new ArrayList<String>();
			valuesJoined = loadSampleFile.get(key);

			// Add "key" Pos0 Pos1, Pos2
			outputFile += lineSplited[0] + SPACE_WRITE + lineSplited[1] + SPACE_WRITE + lineSplited[2];

			// Add the Covar and Response columns
			for (int i = 0; i < valuesJoined.size(); i++) {

				String values = valuesJoined.get(i);
				outputFile += SPACE_WRITE + values;
			}
			outputFile += "\n";
		}
		br.close();

		FileWriter fw = new FileWriter(phasingNewSampleFile);
		BufferedWriter writerInfo = new BufferedWriter(fw);
		writerInfo.write(outputFile);
		writerInfo.close();
		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] newSample startTime: " + startTime);
			System.out.println("\n[DEBUG] newSample endTime: " + stopTime);
			System.out.println("\n[DEBUG] newSample elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of newSample.");

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
	public static void filterHaplotypes(String phasingHapsFile, String phasingSampleFile, String excludedSnpsFile,
			String filteredLogFile, String filteredHapsVcfFile, String cmdToStore) throws GuidanceTaskException {

		String shapeitBinary = loadFromEnvironment(SHAPEITBINARY, HEADER_FILTER_HAPLOTYPES);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running filterHaplotypes with parameters:");
			System.out.println("[DEBUG] \t- shapeitBinary             : " + SHAPEITBINARY);
			System.out.println("[DEBUG] \t- Input phasingHapsFile     : " + phasingHapsFile);
			System.out.println("[DEBUG] \t- Input phasingSampleFile   : " + phasingSampleFile);
			System.out.println("[DEBUG] \t- Input excludedSnpsFile    : " + excludedSnpsFile);
			System.out.println("[DEBUG] \t- Output filteredLogFile     : " + filteredLogFile);
			System.out.println("[DEBUG] \t- Output filteredHapsVcfFile : " + filteredHapsVcfFile);

			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}
		long startTime = System.currentTimeMillis();

		String cmd = null;
		cmd = shapeitBinary + " -convert --input-haps " + phasingHapsFile + " " + phasingSampleFile + " --exclude-snp "
				+ excludedSnpsFile + " " + " --output-log " + filteredLogFile + " --output-vcf " + filteredHapsVcfFile;

		if (DEBUG) {
			System.out.println(HEADER_FILTER_HAPLOTYPES + MSG_CMD + cmd);
		}

		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, filteredLogFile + STDOUT_EXTENSION,
					filteredLogFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			System.err.println(HEADER_FILTER_HAPLOTYPES
					+ "Warning executing shapeitProc job in mode -convert, exit value is: " + exitValue);
			System.err.println(HEADER_FILTER_HAPLOTYPES + "                         (This warning is not fatal).");
		}

		try {
			File tempFile = new File(filteredHapsVcfFile + ".tmp");
			FileWriter fw = new FileWriter(tempFile);
			BufferedWriter writer = new BufferedWriter(new FileWriter(tempFile));

			FileReader fr = new FileReader(filteredHapsVcfFile);
			BufferedReader br = new BufferedReader(fr);

			String output = "";
			String line = "";

			while (br.ready()) {
				line = br.readLine();
				// String kk = aux.replaceAll("^23", "X");
				output = line.replaceAll("^23", "X");
				writer.write(output + "\n");
			}

			br.close();
			fr.close();
			fw.close();
			writer.close();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Finally replace the original file.
		FileUtils.move(filteredHapsVcfFile + ".tmp", filteredHapsVcfFile);

		FileUtils.move(filteredLogFile + ".log", filteredLogFile);
		FileUtils.move(filteredHapsVcfFile + ".vcf", filteredHapsVcfFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] filterHaplotypes startTime: " + startTime);
			System.out.println("\n[DEBUG] filterHaplotypes endTime: " + stopTime);
			System.out.println("\n[DEBUG] filterHaplotypes elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of filterHaplotypes.");
		}
	}

	public static void samtoolsBgzip(String input, String output, String cmdToStore)
			throws IOException, InterruptedException, Exception {

		String samToolsBinary = loadFromEnvironment(SAMTOOLSBINARY, HEADER_SAMTOOLS);

		long startTime = System.currentTimeMillis();

		String cmd = samToolsBinary + "/bgzip -f " + input + " " + output;

		if (DEBUG) {
			System.out.println("\n[DEBUG] Command: " + cmd);
			System.out.println(" ");
		}

		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, input + STDOUT_EXTENSION, input + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		if (exitValue != 0) {
			System.err.println("[phasing] Warning executing samtoolsBgzip job, exit value is: " + exitValue);
			System.err.println("[phasing]                         (This warning is not fatal).");
		}

		File tmpOutput = new File(input + ".gz");
		if (tmpOutput.exists()) {
			boolean success = FileUtils.move(input + ".gz", output);
			if (!success) {
				throw new GuidanceTaskException(
						"[samtoolsBgzip] Error, the file " + input + ".gz was not succesfully renamed to " + output);
				// File was not successfully renamed
			}
		} else {
			System.out.println("[samtoolsBgzip] No! " + tmpOutput + " does not exist!");
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] samtoolsBgzip startTime: " + startTime);
			System.out.println("\n[DEBUG] samtoolsBgzip endTime: " + stopTime);
			System.out.println("\n[DEBUG] samtoolsBgzip elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of samtoolsBgzip.");
		}
	}

	public static void samtoolsTabix(String inputGz, String outputTbi, String cmdToStore)
			throws IOException, InterruptedException, Exception {

		// Public static void samtoolsTabix(String filteredHaplotypesVcfFileBgzip,
		// String filteredHaplotypesVcfFileBgzipIndexed,
		// String cmdToStore) throws IOException, InterruptedException, Exception {

		// System.out.println("[samtoolsTabix] filteredHaplotypesVcfFileBgzip is: " +
		// filteredHaplotypesVcfFileBgzip);
		// System.out.println("[samtoolsTabix] filteredHaplotypesVcfFileBgzipIndexed is:
		// " + filteredHaplotypesVcfFileBgzipIndexed);

		String samtoolsBinary = loadFromEnvironment(SAMTOOLSBINARY, HEADER_SAMTOOLS);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running samtoolsTabix with parameters:");
			System.out.println("[DEBUG] \t- samtoolsBinary                                : " + samtoolsBinary);
			System.out.println(
					"[DEBUG] \t- Input                                                                         : "
							+ inputGz);
			System.out.println(
					"[DEBUG] \t- Output                                                                            : "
							+ outputTbi);
			System.out.println("\n");
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
		}
		long startTime = System.currentTimeMillis();

		// String filteredHaplotypesVcfFileBgzipIndexed = filteredHaplotypesVcfFileBgzip
		// + ".tbi";

		String cmd = null;

		cmd = samtoolsBinary + "/tabix -p vcf --force " + inputGz;

		if (DEBUG) {
			System.out.println("\n[DEBUG] Command: " + cmd);
			System.out.println(" ");
		}

		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, outputTbi + STDOUT_EXTENSION, outputTbi + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		if (exitValue != 0) {
			System.err.println("[phasing] Warning executing phasingProc job, exit value is: " + exitValue);
			System.err.println("[phasing]                         (This warning is not fatal).");
		}

		String output = inputGz + ".tbi";
		File tmpOutput = new File(output);
		if (tmpOutput.exists()) {
			boolean success = FileUtils.move(inputGz + ".tbi", outputTbi);
			// boolean success = tmpOutput.renameTo(new File(outputTbi));
			if (!success) {
				throw new Exception("[samtoolsTabix] Error, the file " + outputTbi + " was not succesfully renamed");
			}
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] samtoolsTabix startTime: " + startTime);
			System.out.println("\n[DEBUG] samtoolsTabix endTime: " + stopTime);
			System.out.println("\n[DEBUG] samtoolsTabix elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of samtoolsTabix.");
		}
	}

	public static void imputeWithImputeLow(String gmapFile, String knownHapFile, String legendFile,
			String phasingHapsFile, String phasingSampleFile, String lim1S, String lim2S, String pairsFile,
			String imputeFile, String imputeFileInfo, String imputeFileSummary, String imputeFileWarnings,
			String theChromo, String sex, String cmdToStore) throws GuidanceTaskException {
		imputeWithImpute(gmapFile, knownHapFile, legendFile, phasingHapsFile, phasingSampleFile, lim1S, lim2S,
				pairsFile, imputeFile, imputeFileInfo, imputeFileSummary, imputeFileWarnings, theChromo, sex,
				cmdToStore);
	}

	public static void imputeWithImputeMedium(String gmapFile, String knownHapFile, String legendFile,
			String phasingHapsFile, String phasingSampleFile, String lim1S, String lim2S, String pairsFile,
			String imputeFile, String imputeFileInfo, String imputeFileSummary, String imputeFileWarnings,
			String theChromo, String sex, String cmdToStore) throws GuidanceTaskException {
		imputeWithImpute(gmapFile, knownHapFile, legendFile, phasingHapsFile, phasingSampleFile, lim1S, lim2S,
				pairsFile, imputeFile, imputeFileInfo, imputeFileSummary, imputeFileWarnings, theChromo, sex,
				cmdToStore);
	}

	public static void imputeWithImputeHigh(String gmapFile, String knownHapFile, String legendFile,
			String phasingHapsFile, String phasingSampleFile, String lim1S, String lim2S, String pairsFile,
			String imputeFile, String imputeFileInfo, String imputeFileSummary, String imputeFileWarnings,
			String theChromo, String sex, String cmdToStore) throws GuidanceTaskException {
		imputeWithImpute(gmapFile, knownHapFile, legendFile, phasingHapsFile, phasingSampleFile, lim1S, lim2S,
				pairsFile, imputeFile, imputeFileInfo, imputeFileSummary, imputeFileWarnings, theChromo, sex,
				cmdToStore);
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
	public static void imputeWithImpute(String gmapFile, String knownHapFile, String legendFile, String phasingHapsFile,
			String phasingSampleFile, String lim1S, String lim2S, String pairsFile, String imputeFile,
			String imputeFileInfo, String imputeFileSummary, String imputeFileWarnings, String theChromo, String sex,
			String cmdToStore) throws GuidanceTaskException {

		String impute2Binary = loadFromEnvironment(IMPUTE2BINARY, HEADER_IMPUTE);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running impute with parameters:");
			System.out.println("[DEBUG] \t- impute2Binary             : " + impute2Binary);
			System.out.println("[DEBUG] \t- Input gmapFile            : " + gmapFile);
			System.out.println("[DEBUG] \t- Input knownHapFile        : " + knownHapFile);
			System.out.println("[DEBUG] \t- Input legendHapFile       : " + legendFile);
			System.out.println("[DEBUG] \t- Input shapeitHapsFile     : " + phasingHapsFile);
			System.out.println("[DEBUG] \t- Input shapeitSampleFile   : " + phasingSampleFile);
			System.out.println("[DEBUG] \t- Input lim1S               : " + lim1S);
			System.out.println("[DEBUG] \t- Input lim2S               : " + lim2S);
			System.out.println("[DEBUG] \t- Input pairsFile           : " + pairsFile);
			System.out.println("[DEBUG] \t- Output imputeFile         : " + imputeFile);
			System.out.println("[DEBUG] \t- Output imputeFileInfo     : " + imputeFileInfo);
			System.out.println("[DEBUG] \t- Output imputeFileSummary  : " + imputeFileSummary);
			System.out.println("[DEBUG] \t- Output imputeFileWarnings : " + imputeFileWarnings);
			System.out.println("[DEBUG] \t- Input  theChromo          : " + theChromo);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// We have to make sure whether we are using renamed files of the original gz
		// files.
		// We detect this situation by scanning the last three characters:
		String extensionHaps = phasingHapsFile.substring(Math.max(0, phasingHapsFile.length() - 3));

		String phasingHapsFileGz = null;
		if (extensionHaps.equals(".gz")) {
			phasingHapsFileGz = phasingHapsFile;
		} else {
			// If shapeitHapsFile exists, then shapeitHapsFileGz exists also.
			phasingHapsFileGz = phasingHapsFile + ".gz";
		}

		// We erase the extension .gz
		imputeFile = imputeFile.substring(0, imputeFile.length() - 3);

		String cmd = null;
		if (theChromo.equals(CHR_23)) {

			// if (sex.equals(SEX1)) {
			cmd = impute2Binary + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l " + legendFile
					+ " -known_haps_g " + phasingHapsFileGz + " -sample_g " + phasingSampleFile + " -int " + lim1S + " "
					+ lim2S + "  -chrX -exclude_snps_g " + pairsFile + " -impute_excluded -Ne 20000 -o " + imputeFile
					+ " -i " + imputeFileInfo + " -r " + imputeFileSummary + " -w " + imputeFileWarnings
					+ " -no_sample_qc_info -o_gz";
			// } else if (sex.equals(SEX2)) {
			// } else {
			// cmd = impute2Binary + " -use_prephased_g -m " + gmapFile + " -h " +
			// knownHapFile + " -l " + legendFile
			// + " -known_haps_g " + phasingHapsFileGz + " -sample_g " + phasingSampleFile +
			// " -int " + lim1S
			// + " " + lim2S + " -chrX -exclude_snps_g " + pairsFile + " -impute_excluded
			// -Ne 20000 -o "
			// + imputeFile + " -i " + imputeFileInfo + " -r " + imputeFileSummary + " -w "
			// + imputeFileWarnings + " -no_sample_qc_info -o_gz";
			// }

			theChromo = "X";

		} else {
			cmd = impute2Binary + " -use_prephased_g -m " + gmapFile + " -h " + knownHapFile + " -l " + legendFile
					+ " -known_haps_g " + phasingHapsFileGz + " -int " + lim1S + " " + lim2S + " -exclude_snps_g "
					+ pairsFile + " -impute_excluded -Ne 20000 -o " + imputeFile + " -i " + imputeFileInfo + " -r "
					+ imputeFileSummary + " -w " + imputeFileWarnings + " -no_sample_qc_info -o_gz";
		}

		if (DEBUG) {
			System.out.println(HEADER_IMPUTE + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, imputeFile + STDOUT_EXTENSION, imputeFile + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			System.err.println(HEADER_IMPUTE + " Warning executing imputeProc job, exit value is: " + exitValue);
			System.err.println(HEADER_IMPUTE + "                        (This warning is not fatal).");
		}

		// With the -o_gz option in the command, the outputs are imputeFile.gz
		// If there is not output in the impute process. Then we have to create some
		// empty outputs
		String imputeGZFile = imputeFile + ".gz";
		if (!new File(imputeGZFile).exists()) {
			try {
				if (!FileUtils.createEmptyFile(imputeFile, HEADER_IMPUTE)) {
					throw new GuidanceTaskException(HEADER_IMPUTE + ERROR_FILE_CREATION + imputeFile + FILE_SUFFIX);
				}
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}
			FileUtils.gzipFile(imputeFile, imputeFile + ".gz");
		}

		// The result has the gz extension
		imputeFile = imputeFile + ".gz";

		// boolean success = FileUtils.move(imputeGZFile, imputeFile);

		// if (!success)

		// {
		// throw new GuidanceTaskException(HEADER_IMPUTE + ERROR_ON_FILE +
		// imputeGZFile);
		// }

		try {
			FileUtils.createEmptyFile(imputeFileInfo, HEADER_IMPUTE);
			FileUtils.createEmptyFile(imputeFileSummary, HEADER_IMPUTE);
			FileUtils.createEmptyFile(imputeFileWarnings, HEADER_IMPUTE);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		String imputeFileColumnTransformation = "zcat " + imputeFile + " | awk -v chr=" + theChromo
				+ " '{out=$1 \" \" chr \":\" $3 \"_\" $4 \"_\" $5 ; for(i=3;i<=NF;i++){out=out\" \"$i}; print out}' | gzip > "
				+ imputeFile + "_tmp; mv " + imputeFile + "_tmp " + imputeFile;

		try {
			ProcessUtils.executeWithoutOutputs(imputeFileColumnTransformation);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		String imputeInfoFileColumnTransformation = "head -n1 " + imputeFileInfo + " > " + imputeFileInfo
				+ "_tmp;tail -n +2 " + imputeFileInfo + " | awk -v chr=" + theChromo
				+ " '{out=$1 \" \" chr \":\" $3 \"_\" $4 \"_\" $5 ; for(i=3;i<=NF;i++){out=out\" \"$i}; print out}' >> "
				+ imputeFileInfo + "_tmp; mv " + imputeFileInfo + "_tmp " + imputeFileInfo;

		try {
			ProcessUtils.executeWithoutOutputs(imputeInfoFileColumnTransformation);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		System.out.println(imputeFileColumnTransformation);
		System.out.println(imputeInfoFileColumnTransformation);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] imputeWithImpute startTime: " + startTime);
			System.out.println("\n[DEBUG] imputeWithImpute endTime: " + stopTime);
			System.out.println("\n[DEBUG] imputeWithImpute elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of imputeWithImpute with parameters:");
			System.out.println("[DEBUG] \t- impute2Binary             : " + impute2Binary);
			System.out.println("[DEBUG] \t- Input gmapFile            : " + gmapFile);
			System.out.println("[DEBUG] \t- Input knownHapFile        : " + knownHapFile);
			System.out.println("[DEBUG] \t- Input legendHapFile       : " + legendFile);
			System.out.println("[DEBUG] \t- Input shapeitHapsFile     : " + phasingHapsFile);
			System.out.println("[DEBUG] \t- Input shapeitSampleFile   : " + phasingSampleFile);
			System.out.println("[DEBUG] \t- Input lim1S               : " + lim1S);
			System.out.println("[DEBUG] \t- Input lim2S               : " + lim2S);
			System.out.println("[DEBUG] \t- Input pairsFile           : " + pairsFile);
			System.out.println("[DEBUG] \t- Output imputeFile         : " + imputeFile);
			System.out.println("[DEBUG] \t- Output imputeFileInfo     : " + imputeFileInfo);
			System.out.println("[DEBUG] \t- Output imputeFileSummary  : " + imputeFileSummary);
			System.out.println("[DEBUG] \t- Output imputeFileWarnings : " + imputeFileWarnings);
			System.out.println("[DEBUG] \t- Input  theChromo          : " + theChromo);
		}

	}

	public static void imputeWithMinimacLow(String vcfFile, String filteredHapsVcfFileBgzip, String imputeFile,
			String imputeFileInfo, String imputeFileErate, String imputeFileRec, String imputeFileM3vcf,
			String imputeFileLog, String chrS, String lim1S, String lim2S, String myPrefix, String sex,
			String cmdToStore) throws GuidanceTaskException {

		imputeWithMinimac(vcfFile, filteredHapsVcfFileBgzip, imputeFile, imputeFileInfo, imputeFileErate, imputeFileRec,
				imputeFileM3vcf, imputeFileLog, chrS, lim1S, lim2S, myPrefix, sex, cmdToStore);
	}

	public static void imputeWithMinimacMedium(String vcfFile, String filteredHapsVcfFileBgzip, String imputeFile,
			String imputeFileInfo, String imputeFileErate, String imputeFileRec, String imputeFileM3vcf,
			String imputeFileLog, String chrS, String lim1S, String lim2S, String myPrefix, String sex,
			String cmdToStore) throws GuidanceTaskException {

		imputeWithMinimac(vcfFile, filteredHapsVcfFileBgzip, imputeFile, imputeFileInfo, imputeFileErate, imputeFileRec,
				imputeFileM3vcf, imputeFileLog, chrS, lim1S, lim2S, myPrefix, sex, cmdToStore);
	}

	public static void imputeWithMinimacHigh(String vcfFile, String filteredHapsVcfFileBgzip, String imputeFile,
			String imputeFileInfo, String imputeFileErate, String imputeFileRec, String imputeFileM3vcf,
			String imputeFileLog, String chrS, String lim1S, String lim2S, String myPrefix, String sex,
			String cmdToStore) throws GuidanceTaskException {

		imputeWithMinimac(vcfFile, filteredHapsVcfFileBgzip, imputeFile, imputeFileInfo, imputeFileErate, imputeFileRec,
				imputeFileM3vcf, imputeFileLog, chrS, lim1S, lim2S, myPrefix, sex, cmdToStore);
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
	private static void imputeWithMinimac(String vcfFile, String filteredHapsVcfFileBgzip, String imputeFile,
			String imputeFileInfo, String imputeFileErate, String imputeFileRec, String imputeFileM3vcf,
			String imputeFileLog, String chrS, String lim1S, String lim2S, String myPrefix, String sex,
			String cmdToStore) throws GuidanceTaskException {

		String minimacBinary = loadFromEnvironment(MINIMACBINARY, HEADER_MINIMAC);
		String realPrefix = imputeFileErate.substring(0, imputeFileErate.length() - 6);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running imputation with parameters:");
			System.out.println("[DEBUG] \t- minimacBinary                     : " + minimacBinary);
			System.out.println("[DEBUG] \t- Input vcfFile                     : " + vcfFile);
			System.out.println("[DEBUG] \t- Input filteredHapsVcfFile         : " + filteredHapsVcfFileBgzip);
			System.out.println("[DEBUG] \t- Output imputeMMInfoFile           : " + imputeFileInfo);
			System.out.println("[DEBUG] \t- Output imputeMMErateFile          : " + imputeFileErate);
			System.out.println("[DEBUG] \t- Output imputeMMRecFile            : " + imputeFileRec);
			System.out.println("[DEBUG] \t- Output imputeMMDoseVCFFile        : " + imputeFile);
			System.out.println("[DEBUG] \t- Output imputeMMM3VCFFile          : " + imputeFileM3vcf);
			System.out.println("[DEBUG] \t- Output imputeMMLogFile            : " + imputeFileLog);
			System.out.println("[DEBUG] \t- Input chrS                        : " + chrS);
			System.out.println("[DEBUG] \t- Input lim1S                       : " + lim1S);
			System.out.println("[DEBUG] \t- Input lim2S                       : " + lim2S);
			System.out.println("[DEBUG] \t- Input myPrefix                    : " + myPrefix);
			System.out.println("\n");
			System.out.println("[DEBUG] \t- Command                           : " + cmdToStore);
		}

		long startTime = System.currentTimeMillis();

		String cmd = null;
		// Submitting the impute task per chunk
		if (chrS.equals("23")) {
			//if (sex.equals(SEX1)) {
				cmd = minimacBinary + " --refHaps " + vcfFile + " --haps " + filteredHapsVcfFileBgzip + " --start "
						+ lim1S + " --end " + lim2S + " --chr X --window 500000 --prefix " + realPrefix
						+ " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
			//} else if (sex.equals(SEX2)) {
			//} else {
			//	cmd = minimacBinary + " --refHaps " + vcfFile + " --haps " + filteredHapsVcfFileBgzip + " --start "
			//			+ lim1S + " --end " + lim2S + " --chr X --window 500000 --prefix " + realPrefix
			//			+ " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
			//}
		} else {
			cmd = minimacBinary + " --refHaps " + vcfFile + " --haps " + filteredHapsVcfFileBgzip + " --start " + lim1S
					+ " --end " + lim2S + " --chr " + chrS + " --window 500000 --prefix " + realPrefix
					+ " --log --allTypedSites --noPhoneHome --format GT,DS,GP --nobgzip";
		}

		if (DEBUG) {
			System.out.println(HEADER_MINIMAC + MSG_CMD + cmd);
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.execute(cmd, filteredHapsVcfFileBgzip + STDOUT_EXTENSION,
					filteredHapsVcfFileBgzip + STDERR_EXTENSION);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			System.err.println(HEADER_MINIMAC + " Warning executing minimacProc job, exit value is: " + exitValue);
			System.err.println(HEADER_MINIMAC + "                        (This warning is not fatal).");
		}

		// If there is not output in the impute process. Then we have to create some
		// empty outputs.

		File MMDoseVCFFile = new File(realPrefix + ".dose.vcf");
		// If this file does not exists, we assume that all the others neither exist
		if (!MMDoseVCFFile.exists()) {
			// Here we create all the new files
			try {
				FileUtils.createEmptyFile(realPrefix + ".m3vcf", HEADER_MINIMAC);

				FileUtils.createEmptyFile(realPrefix + ".dose.vcf", HEADER_MINIMAC);

				FileUtils.createEmptyFile(realPrefix + ".rec", HEADER_MINIMAC);

				FileUtils.createEmptyFile(realPrefix + ".erate", HEADER_MINIMAC);

				FileUtils.createEmptyFile(realPrefix + ".info", HEADER_MINIMAC);

				FileUtils.createEmptyFile(realPrefix + ".logfile", HEADER_MINIMAC);

			} catch (IOException gte) {
				throw new GuidanceTaskException(gte);
			}

		} else if (MMDoseVCFFile.exists()) {
			// We assume here that all the files have been created
			try {
				FileWriter fw = new FileWriter(MMDoseVCFFile.getAbsoluteFile(), true);
				BufferedWriter bw = new BufferedWriter(fw);
				bw.write("\n");
				bw.close();
				fw.close();
			} catch (IOException gte) {
				throw new GuidanceTaskException(gte);
			}
			FileUtils.move(realPrefix + ".m3vcf", realPrefix + ".m3vcf");

			FileUtils.move(realPrefix + ".dose.vcf", realPrefix + ".dose.vcf");

			FileUtils.move(realPrefix + ".rec", imputeFileRec);

			FileUtils.move(realPrefix + ".erate", imputeFileErate);

			FileUtils.move(realPrefix + ".info", imputeFileInfo);

			FileUtils.move(realPrefix + ".logfile", imputeFileLog);

			/*
			 * try {
			 * 
			 * FileUtils.createEmptyFile(realPrefix + ".rec", HEADER_MINIMAC);
			 * 
			 * FileUtils.createEmptyFile(realPrefix + ".erate", HEADER_MINIMAC);
			 * 
			 * FileUtils.createEmptyFile(realPrefix + ".logfile", HEADER_MINIMAC);
			 * 
			 * FileUtils.createEmptyFile(realPrefix + ".info", HEADER_MINIMAC);
			 * 
			 * } catch (IOException ioe) { throw new GuidanceTaskException(ioe); }
			 */
		}
		FileUtils.bgzipFile(realPrefix + ".dose.vcf", imputeFile);
		FileUtils.gzipFile(realPrefix + ".m3vcf", imputeFileM3vcf);

		// FileUtils.delete(realPrefix + ".dose.vcf");
		// FileUtils.delete(realPrefix + ".m3vcf");

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] imputeWithMinimac startTime: " + startTime);
			System.out.println("\n[DEBUG] imputeWithMinimac endTime: " + stopTime);
			System.out.println("\n[DEBUG] imputeWithMinimac elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of imputeWithMinimac with parameters:");
		}

	}

	// This is not a task!! Calling it will imply a sincronization!!!
	// If someday this can be a task, combinedTopHits -> FILE_IN_ARRAY
	// topHitsAllPheno -> FILE_OUT
	public static void generateTopHitsAllPhenos(List<String> combinedTopHits, String topHitsAllPheno)
			throws GuidanceTaskException {
		String rScriptBinDir = loadFromEnvironment(RSCRIPTBINDIR, HEADER_PHENO);
		String rScriptDir = loadFromEnvironment(RSCRIPTDIR, HEADER_PHENO);

		String combinedTopHitsString = combinedTopHits.get(0);
		for (int i = 1; i < combinedTopHits.size(); ++i) {
			combinedTopHitsString += ("," + combinedTopHits.get(i));
		}

		String command = rScriptBinDir + "Rscript --verbose " + rScriptDir + "/tophits_all_phenotypes.R "
				+ combinedTopHitsString + " " + topHitsAllPheno;

		long startTime = 0;
		if (DEBUG) {
			startTime = System.currentTimeMillis();
			System.out.println("\n[DEBUG] [pheno] Launched command:                 : " + command);
		}

		try {
			ProcessUtils.executeWithoutOutputs(command);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// FileUtils.getFile(topHitsAllPheno);

		if (DEBUG) {
			long stopTime = System.currentTimeMillis();
			long elapsedTime = (stopTime - startTime) / 1_000;
			System.out.println("\n[DEBUG] [pheno] generateTopHitsAllPhenos startTime   : " + startTime);
			System.out.println("\n[DEBUG] [pheno] generateTopHitsAllPhenos endTime     : " + stopTime);
			System.out.println("\n[DEBUG] [pheno] generateTopHitsAllPhenos elapsedTime : " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] [pheno] Finished execution of generateTopHitsAllPhenos.");
		}
	}

	public static void generateMergedPhenoTopHits(String topHitsAllPheno, String condensedFile, String mergedPhenoFile,
			String pheno) throws GuidanceTaskException {

		String rScriptBinDir = loadFromEnvironment(RSCRIPTBINDIR, HEADER_PHENO);
		String rScriptDir = loadFromEnvironment(RSCRIPTDIR, HEADER_PHENO);

		String command = rScriptBinDir + "Rscript --verbose " + rScriptDir + "/merging_tophits_all_pheno.R "
				+ topHitsAllPheno + " " + condensedFile + " " + mergedPhenoFile + " " + pheno;

		long startTime = 0;
		if (DEBUG) {
			startTime = System.currentTimeMillis();
			System.out.println("\n[DEBUG] [pheno] Launched command:                 : " + command);
		}

		try {
			ProcessUtils.executeWithoutOutputs(command);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		if (DEBUG) {
			long stopTime = System.currentTimeMillis();
			long elapsedTime = (stopTime - startTime) / 1_000;
			System.out.println("\n[DEBUG] [pheno] generateMergedPhenoTopHits startTime   : " + startTime);
			System.out.println("\n[DEBUG] [pheno] generateMergedPhenoTopHits endTime     : " + stopTime);
			System.out
					.println("\n[DEBUG] [pheno] generateMergedPhenoTopHits elapsedTime : " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] [pheno] Finished execution of generateMergedPhenoTopHits.");
		}
	}

	// This is not a task!! Calling it will imply a sincronization!!!
	// If someday this can be a task, phenoMergedTopHits -> FILE_IN_ARRAY
	// crossPhenoAll -> FILE_OUT
	// crossPhenoRanges -> FILE_OUT
	// crossPhenoTopVariants -> FILE_OUT
	// pvaThreshold -> STRING_IN
	public static void computeCrossPheno(List<String> phenoMergedTopHits, String crossPhenoAll, String crossPhenoRanges,
			String crossPhenoTopVariants, String pvaThreshold) throws GuidanceTaskException {

		String rScriptBinDir = loadFromEnvironment(RSCRIPTBINDIR, HEADER_PHENO);
		String rScriptDir = loadFromEnvironment(RSCRIPTDIR, HEADER_PHENO);

		String mergedTopHitsString = phenoMergedTopHits.get(0);
		for (int i = 1; i < phenoMergedTopHits.size(); ++i) {
			mergedTopHitsString += ("," + phenoMergedTopHits.get(i));
		}

		String command = rScriptBinDir + "Rscript --verbose " + rScriptDir + "/crossphenotype.R " + mergedTopHitsString
				+ " " + crossPhenoAll + " " + crossPhenoRanges + " " + crossPhenoTopVariants + " " + pvaThreshold;

		long startTime = 0;
		if (DEBUG) {
			startTime = System.currentTimeMillis();
			System.out.println("\n[DEBUG] [pheno] Launched command:                 : " + command);
		}

		try {
			ProcessUtils.executeWithoutOutputs(command);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		if (DEBUG) {
			long stopTime = System.currentTimeMillis();
			long elapsedTime = (stopTime - startTime) / 1_000;
			System.out.println("\n[DEBUG] [pheno] computeCrossPheno startTime   : " + startTime);
			System.out.println("\n[DEBUG] [pheno] computeCrossPheno endTime     : " + stopTime);
			System.out.println("\n[DEBUG] [pheno] computeCrossPheno elapsedTime : " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] [pheno] Finished execution of computeCrossPheno.");
		}
	}

	public static void generateCondensedAndTopHitsFile(String filteredFile, String filteredMalesFile,
			String filteredFemalesFile, String condensedFile, String topHitsFile, String crossRangesFile,
			String pvaThresholdStr, String cmdToStore) throws GuidanceTaskException {

		String command = null;
		long startTime = System.currentTimeMillis();
		String condensedPlain = condensedFile.substring(0, condensedFile.length() - 3);
		String topHitsPlain = topHitsFile.substring(0, topHitsFile.length() - 3);
		String crossRangesPlain = crossRangesFile.substring(0, crossRangesFile.length() - 3);

		String rScriptBinDir = loadFromEnvironment(RSCRIPTBINDIR, HEADER_GENERATE_QQ_MANHATTAN_PLOTS);
		String rScriptDir = loadFromEnvironment(RSCRIPTDIR, HEADER_GENERATE_QQ_MANHATTAN_PLOTS);

		command = rScriptBinDir + "Rscript --verbose " + rScriptDir + "/condensed_tophits_crossmodel.R " + filteredFile
				+ " " + filteredMalesFile + " " + filteredFemalesFile + " " + condensedPlain + " " + topHitsPlain + " "
				+ crossRangesPlain + " " + pvaThresholdStr;

		if (DEBUG) {
			System.out.println("\n[DEBUG] Launched command:                 : " + command);
		}

		try {
			ProcessUtils.executeWithoutOutputs(command);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		FileUtils.gzipFile(condensedPlain, condensedFile);
		FileUtils.gzipFile(topHitsPlain, topHitsFile);
		FileUtils.gzipFile(crossRangesPlain, crossRangesFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] generateCondensedFile startTime   : " + startTime);
			System.out.println("\n[DEBUG] generateCondensedFile endTime     : " + stopTime);
			System.out.println("\n[DEBUG] generateCondensedFile elapsedTime : " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of generateCondensedFile.");
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
	public static void filterByInfo(String imputationTool, String imputeFileInfo, String inclusionRsIdFile,
			String infoThresholdS, String mafThresholdS, String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running filterByInfo with parameters:");
			System.out.println("[DEBUG] \t- Imputation Tool          : " + imputationTool);
			System.out.println("[DEBUG] \t- Input inputeFileInfo     : " + imputeFileInfo);
			System.out.println("[DEBUG] \t- Output inclusionRsIdFile : " + inclusionRsIdFile);
			System.out.println("[DEBUG] \t- Input threshold          : " + infoThresholdS);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}
		long startTime = System.currentTimeMillis();

		// The position of info and rsId values in the imputeFileInfo
		int infoIndex = 6;
		int rsIdIndex = -1;
		int mafIndex = 4;
		int typeIndex = 7;

		if (imputationTool.equals("impute")) {
			rsIdIndex = 1;
		} else if (imputationTool.equals("minimac")) {
			rsIdIndex = 0;
		} else {
			System.out.println(
					"[filterByInfo]: The Imputation Tool is not valid. Only IMPUTE2 and MINIMAC3 are available");
		}

		// Convert threshold string into thresholdDouble
		double thresholdInfoDouble = Double.parseDouble(infoThresholdS); // store info value in Double format
		double thresholdMafDouble = Double.parseDouble(mafThresholdS);

		File outInclusionRsIdFile = new File(inclusionRsIdFile);
		// Tries to create the file
		if (outInclusionRsIdFile.exists()) {
			outInclusionRsIdFile.delete();
		}

		boolean bool = false;
		try {
			bool = outInclusionRsIdFile.createNewFile();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}
		// Print information about the existence of the file
		System.out.println("\n[DEBUG] \t- Output file " + inclusionRsIdFile + " was succesfuly created? " + bool);

		try (FileReader fr = new FileReader(imputeFileInfo);
				BufferedReader br = new BufferedReader(fr);
				BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(outInclusionRsIdFile))) {

			// We read each line of the imputeFileInfo and put them into string.
			// I read the header
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				String[] splittedLine = line.split(SPACE);
				if (imputationTool.equals("minimac")) {
					String type = splittedLine[typeIndex];

					if (type.equals("Typed_Only") || (type.equals("Genotyped")
							&& (Double.parseDouble(splittedLine[mafIndex]) > thresholdMafDouble))) {

						writerFiltered.write(splittedLine[rsIdIndex]);
						writerFiltered.newLine();
					} else {
						double info = Double.parseDouble(splittedLine[infoIndex]); // store info value in Double format

						// int retval = Double.compare(info, thresholdDouble);

						if ((info >= thresholdInfoDouble)
								&& (Double.parseDouble(splittedLine[mafIndex]) > thresholdMafDouble)) {
							// The info value is greater or equal to the threshold, then store the rsID into
							// the output file.
							writerFiltered.write(splittedLine[rsIdIndex]);
							writerFiltered.newLine();
						}
					}

				} else {

					double info = Double.parseDouble(splittedLine[infoIndex]); // store info value in Double format

					// Store rsID into filteredFile if info >= threshold
					// int retval = Double.compare(info, thresholdDouble);

					if (info >= thresholdInfoDouble) {
						// The info value is greater or equal to the threshold, then store the rsID into
						// the output file.
						writerFiltered.write(splittedLine[rsIdIndex]);
						writerFiltered.newLine();
					}
				}
			}

			writerFiltered.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void filterByAll(String imputationTool, String inputFile, String outputFile,
			/* String outputCondensedFile, */ String mafThresholdS, String infoThresholdS, String hweCohortThresholdS,
			String hweCasesThresholdS, String hweControlsThresholdS, String sex, String rpanelName, String cmdToStore)
			throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running filterByAll with parameters:");
			System.out.println("[DEBUG] \t- Input imputationTool          : " + imputationTool);
			System.out.println("[DEBUG] \t- Input summaryFile             : " + inputFile);
			System.out.println("[DEBUG] \t- Output outputFile             : " + outputFile);
			System.out.println("[DEBUG] \t- Input maf threshold           : " + mafThresholdS);
			System.out.println("[DEBUG] \t- Input info threshold          : " + infoThresholdS);
			System.out.println("[DEBUG] \t- Input hwe cohort threshold    : " + hweCohortThresholdS);
			System.out.println("[DEBUG] \t- Input hwe controls threshold  : " + hweCasesThresholdS);
			System.out.println("[DEBUG] \t- Input hwe cases threshold     : " + hweControlsThresholdS);
			System.out.println("[DEBUG] \t- Input sex                     : " + sex);
			System.out.println("[DEBUG] \t- Input rpanelName              : " + rpanelName);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}
		long startTime = System.currentTimeMillis();

		// The output will be first generated in plain format, then compressed
		String plainOutputFile = outputFile.substring(0, outputFile.length() - 3);

		// Convert threshold string into thresholdDouble
		Double mafThreshold = Double.parseDouble(mafThresholdS);
		Double infoThreshold = Double.parseDouble(infoThresholdS);
		Double hweCohortThreshold = Double.parseDouble(hweCohortThresholdS);
		Double hweCasesThreshold = Double.parseDouble(hweCasesThresholdS);
		Double hweControlsThreshold = Double.parseDouble(hweControlsThresholdS);

		HashMap<String, Integer> inputFileHashTableIndex = new HashMap<>();

		try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(inputFile));
				InputStreamReader decoder = new InputStreamReader(inputGz);
				BufferedReader br = new BufferedReader(decoder);
				BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(plainOutputFile));) {

			// I read the header
			String line = br.readLine();
			// Put the header in the output file.
			writerFiltered.write(line + "\trefpanel");
			writerFiltered.newLine();

			inputFileHashTableIndex = Headers.createHashWithHeader(line, TAB);

			while ((line = br.readLine()) != null) {
				String[] splittedLine = line.split(TAB);// delimiter I assume single space.

				String chromo = splittedLine[inputFileHashTableIndex.get("chr")];
				String infoS = null;
				if (imputationTool.equals("impute") || sex.equals(SEX1) || sex.equals(SEX2)) {
					infoS = splittedLine[inputFileHashTableIndex.get("info_all")];
				} else if (imputationTool.equals("minimac")) {
					infoS = splittedLine[inputFileHashTableIndex.get("info_all")];
				}

				// String alleleA = splittedLine[inputFileHashTableIndex.get("alleleA")];
				// String alleleB = splittedLine[inputFileHashTableIndex.get("alleleB")];

				// We start with these values for hwe values just to allows the X chromosome to
				// pass the if statement of
				// the
				// next lines
				// Just remember that hwe filtering when chromo X is being processed does not
				// make sense.
				String hwe_cohortS = "1.0";
				String hwe_casesS = "1.0";
				String hwe_controlsS = "1.0";
				String cases_mafS = "1.0";
				String controls_mafS = "1.0";

				if (!chromo.equals(CHR_23)) {
					hwe_cohortS = splittedLine[inputFileHashTableIndex.get("cohort_1_hwe")];
					try {
						hwe_casesS = splittedLine[inputFileHashTableIndex.get("cases_hwe")];
						hwe_controlsS = splittedLine[inputFileHashTableIndex.get("controls_hwe")];
					} catch (Exception e) {
						if (DEBUG) {
							System.out.println("There are not the columns cases_hwe and controls_hew in this case");
						}
					}
				}

				try {
					cases_mafS = splittedLine[inputFileHashTableIndex.get("cases_maf")];
					controls_mafS = splittedLine[inputFileHashTableIndex.get("controls_maf")];
				} catch (Exception e) {
					if (DEBUG) {
						System.out.println("There are not the columns cases_maf and controls_maf in this case");
					}
					cases_mafS = splittedLine[inputFileHashTableIndex.get("all_maf")];
					controls_mafS = splittedLine[inputFileHashTableIndex.get("all_maf")];
				}
				// String position = splittedLine[inputFileHashTableIndex.get("position")];
				String pva = splittedLine[inputFileHashTableIndex.get("frequentist_add_pvalue")];

				if (!cases_mafS.equals("NA") && !controls_mafS.equals("NA") && !infoS.equals("NA")
						&& !hwe_cohortS.equals("NA") && !hwe_casesS.equals("NA") && !hwe_controlsS.equals("NA")
						&& !pva.equals("NA")) {

					Double cases_maf = Double.parseDouble(cases_mafS);
					Double controls_maf = Double.parseDouble(controls_mafS);

					Double infoD = 0.0;
					String infoStr = " ";
					Double hweCohort = 1.0;
					Double hweCases = 1.0;
					Double hweControls = 1.0;

					if (!chromo.equals(CHR_23)) {
						hweCohort = Double.parseDouble(hwe_cohortS);
						hweCases = Double.parseDouble(hwe_casesS);
						hweControls = Double.parseDouble(hwe_controlsS);
					}

					if (imputationTool.equals("impute")) {
						infoD = Double.parseDouble(infoS);
					} else if (imputationTool.equals("minimac")) {
						if (!infoS.equals("-")) {
							infoD = Double.parseDouble(infoS);
						} else {
							infoD = 1.0;
						}
					}

					if (infoD >= infoThreshold || infoStr.equals("-")) {
						if (cases_maf >= mafThreshold && controls_maf >= mafThreshold && // VERIFICAR LA CONDICION
								hweCohort >= hweCohortThreshold && hweCases >= hweCasesThreshold
								&& hweControls >= hweControlsThreshold) {

							writerFiltered.write(line + TAB + rpanelName);
							writerFiltered.newLine();
						}
					}
				}
			}

			writerFiltered.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Should never be empty since we always write the header but lets create the
		// file if empty
		try {
			FileUtils.createEmptyFile(plainOutputFile, "[filterByAll]");
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Then, we compress the output files and clean them
		FileUtils.gzipFile(plainOutputFile, outputFile);
		FileUtils.delete(plainOutputFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void jointFilteredByAllFiles(String filteredByAllA, String filteredByAllB, String filteredByAllC,
			String rpanelName, String rpanelFlag, String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running jointFilteredByAllFiles with parameters:");
			System.out.println("[DEBUG] \t- Input filteredByAllA          : " + filteredByAllA);
			System.out.println("[DEBUG] \t- Input filteredByAllB          : " + filteredByAllB);
			System.out.println("[DEBUG] \t- Output filteredByAllC         : " + filteredByAllC);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// First we generate the plain file, then we will compress it
		String plainFilteredByAll = filteredByAllC.substring(0, filteredByAllC.length() - 3);
		try (BufferedWriter writerFiltered = new BufferedWriter(new FileWriter(plainFilteredByAll))) {
			try (GZIPInputStream filteredByAllGz = new GZIPInputStream(new FileInputStream(filteredByAllA));
					InputStreamReader decoder = new InputStreamReader(filteredByAllGz);
					BufferedReader br = new BufferedReader(decoder)) {

				boolean putRefpanel = false;

				// I read the header
				String line = br.readLine();
				// I put the refpanel column in the header:
				String[] splittedHeader = line.split(TAB);
				if (!splittedHeader[splittedHeader.length - 1].equals("refpanel")) {
					line = line + "\trefpanel";
					putRefpanel = true;
				}

				// Put the header in the output file.
				writerFiltered.write(line);
				writerFiltered.newLine();

				while ((line = br.readLine()) != null) {
					if (putRefpanel) {
						line = line + TAB + rpanelName;
					}
					writerFiltered.write(line);
					writerFiltered.newLine();
				}
			}

			// Do the same with the filteredByAllB file if this is different to the
			// filteredByAllA file
			// OK, I explain now: The only way filteredByAllB = filteredByAllA is when there
			// is only one chromosome to
			// process.
			// In that case, in the main program, we put the same file as filteredByAllA and
			// filteredByAllB.
			if (!filteredByAllA.equals(filteredByAllB)) {
				try (GZIPInputStream filteredByAllGz = new GZIPInputStream(new FileInputStream(filteredByAllB));
						InputStreamReader decoder = new InputStreamReader(filteredByAllGz);
						BufferedReader br = new BufferedReader(decoder)) {

					boolean putRefpanel = false;
					// I read the header and skip it.
					String line = br.readLine();
					String[] splittedHeader = line.split(TAB);
					if (!splittedHeader[splittedHeader.length - 1].equals("refpanel")) {
						putRefpanel = true;
					}

					while ((line = br.readLine()) != null) {
						if (putRefpanel) {
							line = line + TAB + rpanelName;
						}

						writerFiltered.write(line);
						writerFiltered.newLine();
					}
				}
			}

			writerFiltered.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// We create the file if empty
		try {
			FileUtils.createEmptyFile(plainFilteredByAll, "[jointFilteredByAllFiles]");
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Then, we create the gz file and erase the plain file
		FileUtils.gzipFile(plainFilteredByAll, filteredByAllC);
		FileUtils.delete(plainFilteredByAll);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// We first write the plain text, then we will compress it
		String plainOutCondensed = outputFile.substring(0, outputFile.length() - 3);

		try (BufferedWriter writerCondensed = new BufferedWriter(new FileWriter(plainOutCondensed))) {
			try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(inputAFile));
					InputStreamReader decoder = new InputStreamReader(inputGz);
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
			// OK, I explain now: The only way inputAFile = inputBFile is when there is only
			// one chromosome to process.
			// In that case, in the main program, we put he same file as inputAFile and
			// inputBFile.
			if (!inputAFile.equals(inputBFile)) {
				// Now the next file: inputBFile
				try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(inputBFile));
						InputStreamReader decoder = new InputStreamReader(inputGz);
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

		// We create the file if empty
		try {
			FileUtils.createEmptyFile(plainOutCondensed, "[jointCondensedFiles]");
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Then, we create the gz file and rename it
		FileUtils.gzipFile(plainOutCondensed, outputFile);
		FileUtils.delete(plainOutCondensed);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] jointCondensedFiles startTime: " + startTime);
			System.out.println("\n[DEBUG] jointCondensedFiles endTime: " + stopTime);
			System.out.println("\n[DEBUG] jointCondensedFiles elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of jointCondensedFiles");
		}
	}

	/**
	 * Method to combine panels complex
	 * 
	 * @param resultsPanelA
	 * @param resultsPanelB
	 * @param resultsPanelC
	 * @param lim1f
	 * @param lim2
	 * @param cmdToStore
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws Exception
	 */
	public static void combinePanelsComplex(String resultsPanelA, String resultsPanelB, String resultsPanelC, int lim1,
			int lim2, String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running combinePanelsComplex with parameters:");
			System.out.println("[DEBUG] \t- resultsPanelA             : " + resultsPanelA);
			System.out.println("[DEBUG] \t- resultsPanelB             : " + resultsPanelB);
			System.out.println("[DEBUG] \t- resultsPanelC             : " + resultsPanelC);
			System.out.println("[DEBUG] \t- lim1               : " + lim1);
			System.out.println("[DEBUG] \t- lim2                 : " + lim2);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		int posIdx = 0;
		int a1Idx = 0;
		int a2Idx = 0;
		int chrIdx = 0;
		int infoIdx = 0;

		// First, we uncompress the input files
		String resultsPanelAUnzip = resultsPanelA + ".temp";
		String resultsPanelBUnzip = resultsPanelB + ".temp";
		FileUtils.gunzipFile(resultsPanelA, resultsPanelAUnzip);
		FileUtils.gunzipFile(resultsPanelB, resultsPanelBUnzip);

		HashMap<String, Integer> resultsHashTableIndex = new HashMap<>();

		// Create the first treeMap for the chromo
		TreeMap<String, String> fileTreeMapA = new TreeMap<>();

		final String EMPTY_HEADER = "chr\tposition\trs_id_all\tinfo_all\tcertainty_all\t";
		String finalHeader = EMPTY_HEADER;
		String positionA1A2Chr = null;

		// Treat results Panel A if it is not empty
		if (new File(resultsPanelAUnzip).exists()) {
			try (BufferedReader br = new BufferedReader(new FileReader(resultsPanelAUnzip))) {
				// Get the header
				String header = br.readLine();

				if (!header.equals(EMPTY_HEADER)) {
					// Save header for the end file
					finalHeader = header;
					// Process header
					resultsHashTableIndex = Headers.createHashWithHeader(header, TAB);
					chrIdx = resultsHashTableIndex.get("chr");
					posIdx = resultsHashTableIndex.get("position");
					a1Idx = resultsHashTableIndex.get("alleleA");
					a2Idx = resultsHashTableIndex.get("alleleB");
					infoIdx = resultsHashTableIndex.get("info_all");

					String line = null;
					while ((line = br.readLine()) != null) {
						String[] splitted = line.split(TAB);

						positionA1A2Chr = splitted[posIdx] + "_" + splitted[a1Idx] + "_" + splitted[a2Idx] + "_"
								+ splitted[chrIdx];
						// Now, we put this String into the treemap with the key positionA1A1Chr
						fileTreeMapA.put(positionA1A2Chr, line);
					}
				}
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}
		}

		// Create the second treeMap for the chromo
		TreeMap<String, String> fileTreeMapB = new TreeMap<>();
		// contador=0;

		// Treat results Panel B if it is not empty
		if (new File(resultsPanelBUnzip).exists()) {
			try (BufferedReader br = new BufferedReader(new FileReader(resultsPanelBUnzip))) {
				// Get the header
				String header = br.readLine();

				if (!header.equals(EMPTY_HEADER)) {
					// Store header for end file
					finalHeader = header;
					// Process header
					resultsHashTableIndex = Headers.createHashWithHeader(header, TAB);
					chrIdx = resultsHashTableIndex.get("chr");
					posIdx = resultsHashTableIndex.get("position");
					a1Idx = resultsHashTableIndex.get("alleleA");
					a2Idx = resultsHashTableIndex.get("alleleB");
					infoIdx = resultsHashTableIndex.get("info_all");

					String line = null;
					while ((line = br.readLine()) != null) {
						String[] splitted = line.split(TAB);

						// if(splitted[chrIdx].equals(chromoS)) {
						positionA1A2Chr = splitted[posIdx] + "_" + splitted[a1Idx] + "_" + splitted[a2Idx] + "_"
								+ splitted[chrIdx];
						// Now, we put this String into the treemap with the key positionA1A1Chr
						fileTreeMapB.put(positionA1A2Chr, line);
						// contador++;
						// }
					}
				}
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}
		}
		// System.out.println("\n[DEBUG] We have read the chromo " + chromoS + " from
		// second File. contador = " +
		// contador);

		// A place to store the results of this combining
		Double infoA;
		Double infoB;

		// String posAllelesEqual = null;
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
			positionA1A2Chr = m.getKey();
			String lineA = m.getValue();
			String[] splittedA = lineA.split(TAB);

			String infoAString = splittedA[infoIdx];
			if (infoAString.equals("-")) {
				infoA = 1.0;
			} else {
				infoA = Double.parseDouble(splittedA[infoIdx]);
			}

			// posAllelesEqual = positionA1A2Chr;
			posAllelesReverse = splittedA[posIdx] + "_" + getAllele(splittedA[a1Idx], splittedA[a2Idx], "reverse") + "_"
					+ splittedA[chrIdx];
			posAllelesComplement = splittedA[posIdx] + "_" + getAllele(splittedA[a1Idx], splittedA[a2Idx], "complement")
					+ "_" + splittedA[chrIdx];
			posAllelesComplementAndReverse = splittedA[posIdx] + "_"
					+ getAllele(splittedA[a1Idx], splittedA[a2Idx], "complementAndReverse") + "_" + splittedA[chrIdx];

			// System.out.println("[combinePanelsComplex] " + positionA1A2Chr + " " +
			// posAllelesEqual + " " +
			// posAllelesReverse + " " + posAllelesComplement + " " +
			// posAllelesComplementAndReverse);

			// The same: position, a1 and a2?
			if (fileTreeMapB.containsKey(positionA1A2Chr)) {
				// If the fileTreeMapB contains this positionA1A2Chr combination, then we have
				// to choose
				// the ones that has a better info (that is the ones with greater info).
				String lineB = fileTreeMapB.get(positionA1A2Chr);
				String[] splittedB = lineB.split(TAB);

				String infoBString = splittedB[infoIdx];
				if (infoBString.equals("-")) {
					infoB = 1.0;
				} else {
					infoB = Double.parseDouble(infoBString);
				}

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
				String[] splittedB = lineB.split(TAB);

				String infoBString = splittedB[infoIdx];
				if (infoBString.equals("-")) {
					infoB = 1.0;
				} else {
					infoB = Double.parseDouble(infoBString);
				}

				// Then we have to choose between A and B.
				if (infoA >= infoB) {
					fileTreeMapC.put(positionA1A2Chr, lineA);
				} else {
					fileTreeMapC.put(posAllelesReverse, lineB);
				}
				// Now we remove this value from the fileTreeMapB
				fileTreeMapB.remove(posAllelesReverse);
				// System.out.println("WOW alelos reversos: " + positionA1A2Chr + " " +
				// posAllelesReverse);
			} else if (fileTreeMapB.containsKey(posAllelesComplement)) {
				// If the fileTreeMapB contains this posAllelesComplement, then we have to
				// choose
				// the ones that has a better info (that is the ones with greater info).
				String lineB = fileTreeMapB.get(posAllelesComplement);
				String[] splittedB = lineB.split(TAB);

				String infoBString = splittedB[infoIdx];
				if (infoBString.equals("-")) {
					infoB = 1.0;
				} else {
					infoB = Double.parseDouble(infoBString);
				}

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
				// If the fileTreeMapB contains this posAllelesComplement, then we have to
				// choose
				// the ones that has a better info (that is the ones with greater info).
				String lineB = fileTreeMapB.get(posAllelesComplementAndReverse);
				String[] splittedB = lineB.split(TAB);

				String infoBString = splittedB[infoIdx];
				if (infoBString.equals("-")) {
					infoB = 1.0;
				} else {
					infoB = Double.parseDouble(infoBString);
				}

				// Then we have to choose between A o B.
				if (infoA >= infoB) {
					fileTreeMapC.put(positionA1A2Chr, lineA);
				} else {
					fileTreeMapC.put(posAllelesComplementAndReverse, lineB);
				}
				// Now we remove this value from the fileTreeMapB
				fileTreeMapB.remove(posAllelesComplementAndReverse);
				// System.out.println("WOW alelos complementariosYreversos: " + positionA1A2Chr
				// + " " +
				// posAllelesComplementAndReverse);
			} else {
				// Else means that fileTreeMapB does not contain this SNP or any of its
				// variants.
				// Therefore, we keep the one in fileTreeMapA
				fileTreeMapC.put(positionA1A2Chr, lineA);
				// System.out.println("WOW fileTreeMapB does not contain this SNP: " +
				// positionA1A2Chr);
			}
			// contador++;
		}

		fileTreeMapA.clear();

		// Now we have to put in fileTreeMapC the rest of values that remain in
		// fileTreeMapB.
		// We iterate the fileTreeMapB (the rest of the...)
		mySet = fileTreeMapB.entrySet();
		// Move next key and value of Map by iterator
		iter = mySet.iterator();
		while (iter.hasNext()) {
			Entry<String, String> m = iter.next();
			positionA1A2Chr = m.getKey();
			String lineB = m.getValue();
			// Then we have to store the value in fileTreeMapC
			fileTreeMapC.put(positionA1A2Chr, lineB);
			// contador++;
		}

		fileTreeMapB.clear();

		// Finally we put the fileTreeMapC into the plain output file and then compress
		// it
		String plainResultsPanelC = resultsPanelC.substring(0, resultsPanelC.length() - 3);
		File resultFile = new File(plainResultsPanelC);
		try {
			resultFile.createNewFile();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultFile))) {
			// We print the header which is the same always
			writer.write(finalHeader);
			writer.newLine();

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
		// System.out.println("\n[DEBUG] We have stored snps from chromosome " + chromoS
		// + " in the output file");

		// Then, we create the gz file and rename it
		FileUtils.gzipFile(plainResultsPanelC, resultsPanelC);

		// Erase tmp files
		if (new File(resultsPanelAUnzip).exists()) {
			FileUtils.delete(resultsPanelAUnzip);
		}
		if (new File(resultsPanelBUnzip).exists()) {
			FileUtils.delete(resultsPanelBUnzip);
		}
		FileUtils.delete(plainResultsPanelC);

		System.out.println("\n[DEBUG] Finished all chromosomes");

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void combineCondensedFiles(String filteredA, String filteredX, String combinedCondensedFile,
			String mafThresholdS, String infoThresholdS, String hweCohortThresholdS, String hweCasesThresholdS,
			String hweControlsThresholdS, String cmdToStore) throws GuidanceTaskException {

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
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}
		long startTime = System.currentTimeMillis();

		// Convert threshold string into thresholdDouble
		Double mafThreshold = Double.parseDouble(mafThresholdS);
		Double infoThreshold = Double.parseDouble(infoThresholdS);
		Double hweCohortThreshold = Double.parseDouble(hweCohortThresholdS);
		Double hweCasesThreshold = Double.parseDouble(hweCasesThresholdS);
		Double hweControlsThreshold = Double.parseDouble(hweControlsThresholdS);

		// Plain output
		String plainCombinedCondensedFile = combinedCondensedFile.substring(0, combinedCondensedFile.length() - 3);

		HashMap<String, Integer> inputFileHashTableIndex = new HashMap<>();
		// HashMap<Integer, String> inputFileHashTableIndexReversed = new HashMap<>();

		try (BufferedWriter writerCondensed = new BufferedWriter(new FileWriter(plainCombinedCondensedFile))) {
			try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(filteredA));
					InputStreamReader decoder = new InputStreamReader(inputGz);
					BufferedReader br = new BufferedReader(decoder)) {

				// I read the header
				String line = br.readLine();

				inputFileHashTableIndex = Headers.createHashWithHeader(line, TAB);
				// inputFileHashTableIndexReversed = createHashWithHeaderReversed(line, TAB);

				String headerCondensed = "CHR\tBP\tP";
				writerCondensed.write(headerCondensed);
				writerCondensed.newLine();

				while ((line = br.readLine()) != null) {
					String[] splittedLine = line.split(TAB);// delimiter I assume single space.

					String chromo = splittedLine[inputFileHashTableIndex.get("chr")];

					String infoS = splittedLine[inputFileHashTableIndex.get("info_all")];

					// We start with these values for hwe values just to allows the X chromosome to
					// pass the if
					// statement of the next lines
					// Just remember that hwe filtering when chromo X is being processed does not
					// make sense.
					String hwe_cohortS = "1.0";
					String hwe_casesS = "1.0";
					String hwe_controlsS = "1.0";

					if (!chromo.equals(CHR_23)) {
						hwe_cohortS = splittedLine[inputFileHashTableIndex.get("cohort_1_hwe")];
						hwe_casesS = splittedLine[inputFileHashTableIndex.get("cases_hwe")];
						hwe_controlsS = splittedLine[inputFileHashTableIndex.get("controls_hwe")];
					}

					String cases_mafS = splittedLine[inputFileHashTableIndex.get("cases_maf")];
					String controls_mafS = splittedLine[inputFileHashTableIndex.get("controls_maf")];

					String position = splittedLine[inputFileHashTableIndex.get("position")];
					String pva = splittedLine[inputFileHashTableIndex.get("frequentist_add_pvalue")];

					String chrbpb = chromo + TAB + position + TAB + pva;

					if (!cases_mafS.equals("NA") && !controls_mafS.equals("NA") && !infoS.equals("NA")
							&& !hwe_cohortS.equals("NA") && !hwe_casesS.equals("NA") && !hwe_controlsS.equals("NA")
							&& !pva.equals("NA")) {
						Double cases_maf = Double.parseDouble(cases_mafS);
						Double controls_maf = Double.parseDouble(controls_mafS);
						Double info = Double.parseDouble(infoS);
						Double hweCohort = 1.0;
						Double hweCases = 1.0;
						Double hweControls = 1.0;

						if (!chromo.equals(CHR_23)) {
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

			// Now with crh 23
			// If filteredA != filteredX then there is chr23 file (filteredX, therefore we
			// have to include it in the
			// results. Otherwise, there is nothing to do.
			if (!filteredA.equals(filteredX)) {
				try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(filteredX));
						InputStreamReader decoder = new InputStreamReader(inputGz);
						BufferedReader br = new BufferedReader(decoder);) {

					inputFileHashTableIndex = new HashMap<>();
					// inputFileHashTableIndexReversed = new HashMap<>();

					// I read the header
					String line = br.readLine();

					inputFileHashTableIndex = Headers.createHashWithHeader(line, TAB);
					// inputFileHashTableIndexReversed = createHashWithHeaderReversed(line, TAB);

					while ((line = br.readLine()) != null) {
						String[] splittedLine = line.split(TAB);// delimiter I assume single space.

						String chromo = splittedLine[inputFileHashTableIndex.get("chr")];

						String infoS = splittedLine[inputFileHashTableIndex.get("info_all")];

						// We start with these values for hwe values just to allows the X chromosome to
						// pass the if
						// statement of the next lines
						// Just remember that hwe filtering when chromo X is being processed does not
						// make sense.
						String hwe_cohortS = "1.0";
						String hwe_casesS = "1.0";
						String hwe_controlsS = "1.0";

						if (!chromo.equals(CHR_23)) {
							hwe_cohortS = splittedLine[inputFileHashTableIndex.get("cohort_1_hwe")];
							hwe_casesS = splittedLine[inputFileHashTableIndex.get("cases_hwe")];
							hwe_controlsS = splittedLine[inputFileHashTableIndex.get("controls_hwe")];
						}

						String cases_mafS = splittedLine[inputFileHashTableIndex.get("cases_maf")];
						String controls_mafS = splittedLine[inputFileHashTableIndex.get("controls_maf")];

						String position = splittedLine[inputFileHashTableIndex.get("position")];
						String pva = splittedLine[inputFileHashTableIndex.get("frequentist_add_pvalue")];

						String chrbpb = chromo + TAB + position + TAB + pva;

						if (!cases_mafS.equals("NA") && !controls_mafS.equals("NA") && !infoS.equals("NA")
								&& !hwe_cohortS.equals("NA") && !hwe_casesS.equals("NA") && !hwe_controlsS.equals("NA")
								&& !pva.equals("NA")) {
							Double cases_maf = Double.parseDouble(cases_mafS);
							Double controls_maf = Double.parseDouble(controls_mafS);
							Double info = Double.parseDouble(infoS);
							Double hweCohort = 1.0;
							Double hweCases = 1.0;
							Double hweControls = 1.0;

							if (!chromo.equals(CHR_23)) {
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

		// Create the file if empty
		try {
			FileUtils.createEmptyFile(plainCombinedCondensedFile, "[combineCondensedFiles]");
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Compress the output file to the parameter value
		FileUtils.gzipFile(plainCombinedCondensedFile, combinedCondensedFile);
		FileUtils.delete(plainCombinedCondensedFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] combineCondensedFiles startTime: " + startTime);
			System.out.println("\n[DEBUG] combineCondensedFiles endTime: " + stopTime);
			System.out.println("\n[DEBUG] combineCondensedFiles elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of combinedCondensedFiles");
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
	public static void generateTopHits(String resultsFile, String outputTopHitFile, String pvaThreshold,
			String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running generateTopHits with parameters:");
			System.out.println("[DEBUG] \t- resultsFile                : " + resultsFile);
			System.out.println("[DEBUG] \t- outputTopHitFile           : " + outputTopHitFile);
			System.out.println("[DEBUG] \t- pvaThreshold               : " + pvaThreshold);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();
		double pvaThres = Double.parseDouble(pvaThreshold);
		// double pvaThres = Double.parseDouble(splitted[17]);

		// We read each line of the resultsFile and put them into the String
		TreeMap<String, String> fileTreeMap = new TreeMap<>();

		String header = null;
		String newHeader = null;
		HashMap<String, Integer> resultsFileHashTableIndex = new HashMap<>();
		// First: read the header and avoid it
		try (GZIPInputStream inputGz = new GZIPInputStream(new FileInputStream(resultsFile));
				InputStreamReader decoder = new InputStreamReader(inputGz);
				BufferedReader br = new BufferedReader(decoder)) {

			header = br.readLine();
			resultsFileHashTableIndex = Headers.createHashWithHeader(header, TAB);

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
				String[] splitted = line.split(TAB);
				String positionAndRsId = splitted[indexPosition] + "_" + splitted[indexRsId];
				double myPva = Double.parseDouble(splitted[indexPvalue]);

				if (myPva <= pvaThres && myPva > 0.0) {
					// Now, we put this String into the treemap with the key positionAndRsId
					// reducedLine is chr;position;RSID_ALL;MAF;a1;a2;pval
					String reducedLine = splitted[indexChromo] + TAB + splitted[indexPosition] + TAB
							+ splitted[indexRsId] + TAB + splitted[indexAllMaf] + TAB + splitted[indexAlleleA] + TAB
							+ splitted[indexAlleleB] + TAB + splitted[indexPvalue];
					fileTreeMap.put(positionAndRsId, reducedLine);
				}
			}
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Finally we put the fileTreeMap into the output file
		// We will generate the output first and then compress it
		String plainOutputTopHitFile = outputTopHitFile.substring(0, outputTopHitFile.length() - 3);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(plainOutputTopHitFile))) {
			// We print the header which is the same always!.
			writer.write(newHeader);
			// writer.newLine();

			Set<Entry<String, String>> mySet = fileTreeMap.entrySet();
			// Move next key and value of Map by iterator
			Iterator<Entry<String, String>> iter = mySet.iterator();
			while (iter.hasNext()) {
				Entry<String, String> m = iter.next();
				String myLine = m.getValue();

				writer.newLine();
				writer.write(myLine);
			}

			writer.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Compress the output file to the parameter value
		FileUtils.gzipFile(plainOutputTopHitFile, outputTopHitFile);
		FileUtils.delete(plainOutputTopHitFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
			String qqPlotTiffFile, String manhattanPlotTiffFile, String manhattanOption, String thresh,
			String cmdToStore) throws GuidanceTaskException {

		String rScriptBinDir = loadFromEnvironment(RSCRIPTBINDIR, HEADER_GENERATE_QQ_MANHATTAN_PLOTS);
		String rScriptDir = loadFromEnvironment(RSCRIPTDIR, HEADER_GENERATE_QQ_MANHATTAN_PLOTS);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running generateQQManhattanPlots with parameters:");
			System.out.println("[DEBUG] \t- lastCondensedFile             : " + lastCondensedFile);
			System.out.println("[DEBUG] \t- qqPlotFile                    : " + qqPlotFile);
			System.out.println("[DEBUG] \t- manhattanPlotFile             : " + manhattanPlotFile);
			System.out.println("[DEBUG] \t- qqPlotTiffFile                : " + qqPlotTiffFile);
			System.out.println("[DEBUG] \t- manhattanPlotTiffFile         : " + manhattanPlotTiffFile);
			System.out.println("[DEBUG] \t- manhattanPlotType             : " + manhattanOption);
			System.out.println("[DEBUG] \t- threshold                     : " + thresh);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		// First, we have to uncompress the input file
		// String theInputFile = lastCondensedFile + TEMP_EXTENSION;
		// FileUtils.gunzipFile(lastCondensedFile, theInputFile);

		String cmd = null;
		cmd = rScriptBinDir + "/Rscript " + rScriptDir + "/qqplot_manhattan_all_models.R " + lastCondensedFile + " "
				+ qqPlotFile + " " + manhattanPlotFile + " " + qqPlotTiffFile + " " + manhattanPlotTiffFile + " "
				+ manhattanOption + " " + thresh;

		if (DEBUG) {
			System.out.println("\n[DEBUG] Cmd -> " + cmd);
			System.out.println(" ");
		}

		// Execute the command retrieving its exitValue, output and error
		int exitValue = -1;
		try {
			exitValue = ProcessUtils.executeWithoutOutputs(cmd);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Check process exit value
		if (exitValue != 0) {
			throw new GuidanceTaskException(HEADER_GENERATE_QQ_MANHATTAN_PLOTS
					+ " Error executing generateQQManProc job, exit value is: " + exitValue);
		}

		// Clean unziped file
		// FileUtils.delete(theInputFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void snptest(String mergedGenFile, String mergedSampleFile, String snptestOutFileGz,
			String snptestLogFile, String responseVar, String covariables, String theChromo, String cmdToStore)
			throws GuidanceTaskException {

		String snptestBinary = loadFromEnvironment(SNPTESTBINARY, HEADER_SNPTEST);

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running snptest with parameters:");
			System.out.println("[DEBUG] \t- snptestBinary                    : " + snptestBinary);
			System.out.println("[DEBUG] \t- Input mergedGenFile              : " + mergedGenFile);
			System.out.println("[DEBUG] \t- Input mergedSampleFile           : " + mergedSampleFile);
			System.out.println("[DEBUG] \t- Output snptestOutFile            : " + snptestOutFileGz);
			System.out.println("[DEBUG] \t- Output snptestLogFile            : " + snptestLogFile);
			System.out.println("[DEBUG] \t- Input responseVar               : " + responseVar);
			System.out.println("[DEBUG] \t- Input covariables                : " + covariables);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		// Replace commas in the string covariables
		String newStr = covariables.replace(',', ' ');
		if (DEBUG) {
			System.out.println("[DEBUG] \t- Changing covariable format. New covariables : " + newStr);
		}

		long startTime = System.currentTimeMillis();

		// We have to make sure whether we are using renamed files of the original gz
		// files.
		// We detect this situation by scanning the last three characters:
		String extension = mergedGenFile.substring(Math.max(0, mergedGenFile.length() - 3));

		String mergedGenFileGz = null;
		if (extension.equals(".gz")) {
			mergedGenFileGz = mergedGenFile;
		} else {
			// If mergeGenFile exists, then mergedGenFileGz exists also.
			// We reused the mergedGenFileGz
			mergedGenFileGz = mergedGenFile + ".gz";
		}

		// The SNPTest Out parameter is a GZ file. We retrieve its name for the binary
		// output and we will compress it
		// later
		String snptestOutFile = snptestOutFileGz.substring(0, snptestOutFileGz.length() - 3);

		// Before executing snptest, I have to verify that the input mergedGenFile is
		// not empty
		int nBytes = -1;
		try (FileInputStream fis = new FileInputStream(new File(mergedGenFile))) {
			nBytes = fis.read();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		if (nBytes != -1) {
			String cmd = null;
			if (covariables.equals("none")) {
				cmd = snptestBinary + " -data " + mergedGenFileGz + " " + mergedSampleFile + " -o " + snptestOutFile
						+ " -pheno " + responseVar + " -hwe -log " + snptestLogFile;
			} else {
				cmd = snptestBinary + " -data " + mergedGenFileGz + " " + mergedSampleFile + " -o " + snptestOutFile
						+ " -pheno " + responseVar + " -cov_names " + newStr + " -hwe -log " + snptestLogFile;
			}

			// Different parameters for chromo 23 (X) and the rest.
			if (theChromo.equals(CHR_23)) {
				cmd = cmd + " -method newml -assume_chromosome X -stratify_on sex -frequentist 1 ";
			} else {
				cmd = cmd + " -method em -frequentist 1 2 3 4 5 ";
			}

			if (DEBUG) {
				System.out.println(HEADER_SNPTEST + MSG_CMD + cmd);
			}

			// Execute the command retrieving its exitValue, output and error
			int exitValue = -1;
			try {
				exitValue = ProcessUtils.execute(cmd, snptestLogFile + STDOUT_EXTENSION,
						snptestLogFile + STDERR_EXTENSION);
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}

			// Check process exit value
			if (exitValue != 0) {
				System.err.println(HEADER_SNPTEST + "Warning executing snptestProc job, exit value is: " + exitValue);
				System.err.println(HEADER_SNPTEST + "                         (This error is not fatal).");
			}
		}

		// The SNP Test binary does not create an empty file if there are not outputs.
		// Check it
		try {
			FileUtils.createEmptyFile(snptestOutFile, HEADER_SNPTEST);
			FileUtils.createEmptyFile(snptestLogFile, HEADER_SNPTEST);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Then, we create the gz file and rename it to the output parameter
		FileUtils.gzipFile(snptestOutFile, snptestOutFileGz);
		FileUtils.delete(snptestOutFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] snptest startTime: " + startTime);
			System.out.println("\n[DEBUG] snptest endTime: " + stopTime);
			System.out.println("\n[DEBUG] snptest elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of snptest.");
		}
	}

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
	public static void mergeTwoChunks(String reduceFileA, String reduceFileB, String reduceFileC, String cmdToStore)
			throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running mergeTwoChunks with parameters:");
			System.out.println("[DEBUG] \t- Input reduceFileA            : " + reduceFileA);
			System.out.println("[DEBUG] \t- Input reduceFileB            : " + reduceFileB);
			System.out.println("[DEBUG] \t- Output reduceFileC           : " + reduceFileC);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		long startTime = System.currentTimeMillis();

		String reducePlainFile = reduceFileC.substring(0, reduceFileC.length() - 3);

		/*
		 * String command1 = "zcat " + reduceFileA + " > " + reducePlainFile; String
		 * command2 = "zcat " + reduceFileB + " | tail -n +2 >> " + reducePlainFile;
		 * String command3 = "gzip -f " + reducePlainFile;
		 * 
		 * try { ProcessUtils.executeWithoutOutputs(command1); } catch (IOException ioe)
		 * { throw new GuidanceTaskException(ioe); }
		 * 
		 * try { ProcessUtils.executeWithoutOutputs(command2); } catch (IOException ioe)
		 * { throw new GuidanceTaskException(ioe); }
		 * 
		 * try { ProcessUtils.executeWithoutOutputs(command3); } catch (IOException ioe)
		 * { throw new GuidanceTaskException(ioe); }
		 */

		String command = "zcat " + reduceFileA + " > " + reducePlainFile + ";zcat " + reduceFileB + " | tail -n +2 >> "
				+ reducePlainFile + ";gzip -f " + reducePlainFile;

		try {
			ProcessUtils.executeWithoutOutputs(command);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
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
	public static void collectSummary(String chr, String imputeTool, String firstImputeFileInfo, String snptestOutFile,
			String reduceFile, String mafThresholdS, String hweCohortThresholdS, String hweCasesThresholdS,
			String hweControlsThresholdS, String sex, String cmdToStore) throws GuidanceTaskException {

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running collectSummary with parameters:");
			System.out.println("[DEBUG] \t- Input chromosome             : " + chr);
			System.out.println("[DEBUG] \t- Input impute tool            : " + imputeTool);
			System.out.println("[DEBUG] \t- Input casesImputeFileInfo    : " + firstImputeFileInfo);
			System.out.println("[DEBUG] \t- Input snptestOutFile         : " + snptestOutFile);
			System.out.println("[DEBUG] \t- Output reduceFile            : " + reduceFile);
			System.out.println("[DEBUG] \t- Input mafThresholdS          : " + mafThresholdS);
			System.out.println("[DEBUG] \t- Input hweCohortThresholdS    : " + hweCohortThresholdS);
			System.out.println("[DEBUG] \t- Input hweCasesThresholdS     : " + hweCasesThresholdS);
			System.out.println("[DEBUG] \t- Input hweControlsThresholdS  : " + hweControlsThresholdS);
			System.out.println("[DEBUG] \t- Input sex                    : " + sex);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + cmdToStore);
			System.out.println("--------------------------------------");
		}

		FileUtils.recursiveSearch(firstImputeFileInfo);

		long startTime = System.currentTimeMillis();
		// Indexes for impute 2.3.2
		/*
		 * int indexRsId = 1; int indexPosition = 2; int indexInfo = 6; int
		 * indexCertainty= 7;
		 */

		int indexRsId = 0;
		int indexPosition = 0;
		int indexInfo = 0;
		int indexCertainty = 0;
		int indexAlleleA = 0;
		int indexAlleleB = 0;

		int length_entry_assoc_list = 0;

		// genOrBimFile is the mixedGenFile and is not necessary to process it.

		// A place to store the results of the merge
		TreeMap<String, ArrayList<String>> summaryTotal = new TreeMap<>();
		HashMap<String, Integer> snptestHashTableIndex = new HashMap<>();
		HashMap<Integer, String> snptestHashTableIndexReversed = new HashMap<>();

		// We read each line of the firstImputeFileInfo and put them into firstList
		// array of Strings
		TreeMap<String, ArrayList<String>> firstTreeMap = new TreeMap<>();
		String positionAndRsId = null;

		try (FileReader fr = new FileReader(firstImputeFileInfo); BufferedReader br = new BufferedReader(fr)) {
			String line = "";
			// Read the header and avoid the header
			line = br.readLine();
			// ---->
			HashMap<String, Integer> imputeHashTableIndex = new HashMap<>();
			if (line != null && !line.isEmpty()) {
				// If we are here, the file is not empty.
				if (imputeTool.equals("impute") || chr.equals("23")) {
					imputeHashTableIndex = Headers.createHashWithHeader(line, " ");
					indexPosition = imputeHashTableIndex.get("position");
					indexRsId = imputeHashTableIndex.get("rs_id");
					indexInfo = imputeHashTableIndex.get("info");
					indexCertainty = imputeHashTableIndex.get("certainty");
					indexAlleleA = imputeHashTableIndex.get("a0");
					indexAlleleB = imputeHashTableIndex.get("a1");
				} else {
					imputeHashTableIndex = Headers.createHashWithHeader(line, TAB);
					indexRsId = imputeHashTableIndex.get("SNP");
					indexInfo = imputeHashTableIndex.get("Rsq");
					indexCertainty = imputeHashTableIndex.get("AvgCall");
					indexAlleleA = imputeHashTableIndex.get("REF(0)");
					indexAlleleB = imputeHashTableIndex.get("ALT(1)");
				}
			}

			while ((line = br.readLine()) != null) {
				ArrayList<String> firstList = new ArrayList<>();
				// delimiter I assume single space.
				String[] splitted = null;

				// Store Position:Store rsIDCases:Store infoCases:Store certCases
				String positionStr = null;
				if (imputeTool.equals("impute") || chr.equals("23")) {
					splitted = line.split(" ");
					positionStr = splitted[indexPosition];
				} else {
					splitted = line.split(TAB);
					positionStr = splitted[indexRsId].split(":")[1];
				}
				firstList.add(positionStr);
				firstList.add(splitted[indexRsId]);
				firstList.add(splitted[indexInfo]);
				firstList.add(splitted[indexCertainty]);
				firstList.add(splitted[indexAlleleA]);
				firstList.add(splitted[indexAlleleB]);
				positionAndRsId = positionStr + "_" + splitted[indexRsId] + "_" + splitted[indexAlleleA] + "_"
						+ splitted[indexAlleleB];

				// If there is not a previous snp with this combination of position and rsID, we
				// store it.
				if (!firstTreeMap.containsKey(positionAndRsId)) {
					// We, put this in the firstTreeMap
					firstTreeMap.put(positionAndRsId, firstList);
				}
				/*
				 * else { // If there is a snp with this combination we should remove it.
				 * firstTreeMap.remove(positionAndRsId); }
				 */
			}
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// We read each line of the snptestOutFile and put them into assocTreeMap array
		// of Strings
		TreeMap<String, ArrayList<String>> assocTreeMap = new TreeMap<>();
		try (GZIPInputStream snptestOutGz = new GZIPInputStream(new FileInputStream(snptestOutFile));
				InputStreamReader decoder = new InputStreamReader(snptestOutGz);
				BufferedReader br = new BufferedReader(decoder)) {

			String line = null;
			while ((line = br.readLine()) != null) {
				// We have to avoid reading all the comment lines that start with "#" character,
				// and one more that is
				// the header and start with "alternate".
				String[] splitted = line.split(" ");
				String firstHeader = splitted[0];
				char firstChar = line.charAt(0);

				// TODO: the next if is ugly.
				if (firstHeader.equals("alternate_ids")) {
					snptestHashTableIndex = Headers.createHashWithHeader(line, " ");
					snptestHashTableIndexReversed = Headers.createHashWithHeaderReversed(line, " ");
				}

				if ((firstChar != '#') && (firstChar != 'a')) {
					ArrayList<String> assocList = new ArrayList<String>();
					// delimiter I assume single space.
					// String[] splitted = line.split(" ");

					// REVIEW: Should we store position?
					// assocList.add(splitted[3]);

					// We store everything, from the line.
					int index_field;
					for (index_field = 2; index_field < splitted.length; index_field++) {
						assocList.add(splitted[index_field]);
					}

					// Now, store the array of string assocList in the assocTreeMap
					positionAndRsId = splitted[snptestHashTableIndex.get("position")] + "_"
							+ splitted[snptestHashTableIndex.get("rsid")] + "_"
							+ splitted[snptestHashTableIndex.get("alleleA")] + "_"
							+ splitted[snptestHashTableIndex.get("alleleB")];
					if (!assocTreeMap.containsKey(positionAndRsId)) {
						assocTreeMap.put(positionAndRsId, assocList);
						length_entry_assoc_list = assocList.size();
					}
					/*
					 * else { assocTreeMap.remove(positionAndRsId); }
					 */
				}
				// The line does not start nor with "#" neither "alternate"
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
			Entry<String, ArrayList<String>> m = iter.next();
			String firstPositionAndRsId = m.getKey();
			ArrayList<String> firstTmp = m.getValue();

			// The same for assocTreeMap. If found, we get the value, otherwise we get null
			ArrayList<String> assocTmp = new ArrayList<>();
			assocTmp = assocTreeMap.get(firstPositionAndRsId);
			ArrayList<String> summaryTmp = new ArrayList<>();
			try {
				summaryTmp = mergeArrays(chr, firstTmp, assocTmp, length_entry_assoc_list, mafThresholdS,
						hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS);
			} catch (IOException ioe) {
				throw new GuidanceTaskException(ioe);
			}
			summaryTotal.put(firstPositionAndRsId, summaryTmp);
			assocTreeMap.remove(firstPositionAndRsId);
		}

		// ---------------
		// Finally we put the summaryTotal into the outputFile
		// We have to create the outputFile for this combination

		// The output of the method is a GZ file. First we create the plain file, then
		// we compress it
		String reducePlainFile = reduceFile.substring(0, reduceFile.length() - 3);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(reducePlainFile))) {
			// We do not store the first 4 field because they are not necessary or are
			// repeated:
			// These four fields are:
			// alternative_ids, rsid, chromosome, position
			if (snptestHashTableIndexReversed.size() > 4) {
				writer.write("chr\tposition\trs_id_all\tinfo_all\tcertainty_all\t");
				for (int index = 4; index < snptestHashTableIndexReversed.size(); index++) {
					String valueReversed = snptestHashTableIndexReversed.get(index);
					writer.write(valueReversed + TAB);
				}
			} else {
				if (sex.equals(NO_SEX)) {
					writer.write(Headers.constructHeader());
				} else {
					writer.write(Headers.constructHeaderX(sex));
				}
			}
			writer.write(NEW_LINE);

			mySet = summaryTotal.entrySet();
			// Move next key and value of Map by iterator
			iter = mySet.iterator();
			while (iter.hasNext()) {
				// key=value separator this by Map.Entry to get key and value
				Entry<String, ArrayList<String>> m = iter.next();
				ArrayList<String> lineTmp = m.getValue();

				writer.write(lineTmp.get(0));
				for (int j = 1; j < lineTmp.size(); j++) {
					writer.write(TAB + lineTmp.get(j));
				}
				writer.newLine();
			}

			writer.flush();
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// We create the file if it does not exist
		try {
			FileUtils.createEmptyFile(reducePlainFile, "[collectSummary]");
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		// Then, we create the gz file and rename it
		FileUtils.gzipFile(reducePlainFile, reduceFile);
		FileUtils.delete(reducePlainFile);

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] collectSummary startTime: " + startTime);
			System.out.println("\n[DEBUG] collectSummary endTime: " + stopTime);
			System.out.println("\n[DEBUG] collectSummary elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of collectSummary.");
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
			int length_entry_assoc_list, String mafThresholdS, String hweCohortThresholdS, String hweCasesThresholdS,
			String hweControlsThresholdS) throws IOException {

		int real_length_assoc = 67;
		if (!chr.equals(CHR_23)) {
			// real_length_assoc = 69;
			real_length_assoc = 67;
		}

		ArrayList<String> summaryTmp = new ArrayList<>();
		// We just need to put the information of the different arraysList into the
		// returned arrayList
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
			// We do not store the first 4 field because they are not necessary or are
			// repeated:
			// These four fields are:
			// alternative_ids, rsid, chromosome, position
			for (int i = 4; i < assoc_length; i++) {
				summaryTmp.add(assocArray.get(i));
			}
		} else {
			// 67 is the number of field in a mixed assocArray
			for (int i = 6; i < real_length_assoc; i++) {
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

		if (typeAllele.equals("reverse")) {
			return allele2 + "_" + allele1;
		} else if (typeAllele.equals("complement")) {
			return a1Complement + "_" + a2Complement;
		} else if (typeAllele.equals("complementAndReverse")) {
			return a2Complement + "_" + a1Complement;
		} else {
			throw new GuidanceTaskException(
					"Error, the option (" + typeAllele + ") for creating alternative alleles is now correct!!");
		}
	}

	public static void copyFile(String fileA, String fileB) throws IOException, GuidanceTaskException {

		long startTime = System.currentTimeMillis();

		String command = "cp " + fileA + " " + fileB;

		if (DEBUG) {
			System.out.println("\n[DEBUG] Running copyFile with parameters:");
			System.out.println("[DEBUG] \t- Origin file                 : " + fileA);
			System.out.println("[DEBUG] \t- Destination file            : " + fileB);
			System.out.println(NEW_LINE);
			System.out.println("[DEBUG] \t- Command: " + command);
			System.out.println("--------------------------------------");
		}

		try {
			ProcessUtils.executeWithoutOutputs(command);
		} catch (IOException ioe) {
			throw new GuidanceTaskException(ioe);
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = (stopTime - startTime) / 1_000;
		if (DEBUG) {
			System.out.println("\n[DEBUG] copyFile startTime: " + startTime);
			System.out.println("\n[DEBUG] copyFile endTime: " + stopTime);
			System.out.println("\n[DEBUG] copyFile elapsedTime: " + elapsedTime + " seconds");
			System.out.println("\n[DEBUG] Finished execution of mergeTwoChunks.");
		}
	}

	public static void getFile(String runtimeFilename, String realFilename) {
		try {
			FileUtils.copy(runtimeFilename, realFilename);
		} catch (IOException e) {
			System.err.println("[DEBUG] Error when bringing back " + realFilename);
		}
	}

}
