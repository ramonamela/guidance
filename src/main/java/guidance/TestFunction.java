package guidance;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import guidance.GuidanceImpl;
import guidance.exceptions.GuidanceTaskException;
import guidance.files.AssocFiles;
import guidance.files.CombinedPanelsFiles;
import guidance.files.MergeFiles;
import guidance.utils.ParseCmdLine;

public class TestFunction {

	private static final String SEX1 = GuidanceImpl.getSex1();
	private static final String SEX2 = GuidanceImpl.getSex2();
	private static final String NO_SEX = GuidanceImpl.getNoSex();

	private static void createSplitedFilesWraper(String[] args) {
		String chr = args[0];
		String gmapFile = args[1];
		String gmapFileChr = args[2];

		try {
			GuidanceImpl.createSplitedFiles(gmapFile, gmapFileChr, chr, "");
		} catch (Exception e) {
			System.err.println("[Guidance] Exception trying the execution of createSplitedFiles task");
			e.printStackTrace();
		}
	}

	private static void convertFromBedToBedWrapper(String[] args) {
		String chr = args[0];
		String bedFile = args[1];
		String bimFile = args[2];
		String famFile = args[3];
		String mixedBedFile = args[4];
		String mixedBimFile = args[5];
		String mixedFamFile = args[6];
		String mixedBedToBedLogFile = args[7];

		try {
			GuidanceImpl.convertFromBedToBed(bedFile, bimFile, famFile, mixedBedFile, mixedBimFile, mixedFamFile,
					mixedBedToBedLogFile, chr, "");
		} catch (GuidanceTaskException e) {
			System.err.println("[Guidance] Exception trying the execution of convertFromBedToBed task");
			e.printStackTrace();
		}
	}

	private static void splitChr23Wrapper(String[] args) {
		String chr = args[1];
		String inputPrefix = args[2];
		String outputPrefix = args[3];
		String logFile = args[4];
		String sex = args[5];

		String bedFile = inputPrefix + ".bed";
		String bimFile = inputPrefix + ".bim";
		String famFile = inputPrefix + ".fam";
		String bedChr23File = outputPrefix + ".bed";
		String bimChr23File = outputPrefix + ".bim";
		String famChr23File = outputPrefix + ".fam";

		try {
			GuidanceImpl.splitChr23(bedFile, bimFile, famFile, bedChr23File, bimChr23File, famChr23File, logFile, sex,
					chr, "");
		} catch (Exception e) {
			System.err.println("[Guidance] Exception trying the execution of splitChr23 task");
			System.err.println(e.getMessage());
		}

	}

	private static void createRsIdListWrapper(String[] args) {

		String mixedBimFile = args[0];
		String exclCgatFlag = args[1];
		String mixedPairsFile = args[2];
		String inputFormat = args[3];

		try {
			GuidanceImpl.createRsIdList(mixedBimFile, exclCgatFlag, mixedPairsFile, inputFormat, "");
		} catch (GuidanceTaskException gte) {
			System.err.println("[Guidance] Exception trying the execution of createRsIdList task");
			gte.printStackTrace();
		}

	}

	private static void collectSummaryWrapper(String[] args) {
		String chr = args[0];
		String imputeTool = args[1];
		String firstImputeFileInfo = args[2];
		String snptestOutFile = args[3];
		String reduceFile = args[4];
		String mafThresholdS = args[5];
		String hweCohortThresholdS = args[6];
		String hweCasesThresholdS = args[7];
		String hweControlsThresholdS = args[8];
		String sex = args[9];

		try {
			GuidanceImpl.collectSummary(chr, imputeTool, firstImputeFileInfo, snptestOutFile, reduceFile, mafThresholdS,
					hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, sex, "");
		} catch (GuidanceTaskException e) { // TODO Auto-generated catch block
			System.err.println("[Guidance] Exception trying the execution of collectSummary task");
			e.printStackTrace();
		}

	}

	private static void filterByAllWrapper(String[] args) {
		String inputFile = args[0];
		String outputFile = args[1];
		String panelName = args[2];
		String mafThresholdS = args[3];
		String infoThresholdS = args[4];
		String hweCohortThresholdS = args[5];
		String hweCasesThresholdS = args[6];
		String hweControlsThresholdS = args[7];
		String cmdToStore = "";
		try {
			GuidanceImpl.filterByAll(inputFile, outputFile, panelName, mafThresholdS, infoThresholdS,
					hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, "", "", cmdToStore);
		} catch (GuidanceTaskException e) {
			System.err.println("[Guidance] Exception trying the execution of filterByAll task");
			e.printStackTrace();
		}
	}

	private static void mergeOfChunksWrapper(String[] args) {
		String reduceA = args[0];
		String reduceB = args[1];
		String reduceC = args[2];
		String cmdToStore = "";
		try {
			GuidanceImpl.mergeTwoChunks(reduceA, reduceB, reduceC, cmdToStore);
		} catch (GuidanceTaskException e) {
			System.err.println("[Guidance] Exception trying the execution of mergeTwoChunks task");
			e.printStackTrace();
		}
	}

	private static void printPaths(String[] args) throws IOException, GuidanceTaskException, InterruptedException {

		ArrayList<String> listOfCommands = new ArrayList<>();

		ParseCmdLine parsingArgs = new ParseCmdLine(args);

		String outDir = parsingArgs.getOutDir();

		ArrayList<String> rpanelTypes = new ArrayList<>(parsingArgs.getRpanelTypes());

		// Create the names for Association files
		AssocFiles assocFilesInfo = new AssocFiles(parsingArgs, outDir, rpanelTypes);

		// Create the names for Merge files
		MergeFiles mergeFilesInfo = new MergeFiles(parsingArgs, outDir, rpanelTypes);

		CombinedPanelsFiles combinedPanelsFilesInfo = new CombinedPanelsFiles(parsingArgs, outDir, rpanelTypes);

		Guidance.makeCombinePanels(parsingArgs, assocFilesInfo, mergeFilesInfo, combinedPanelsFilesInfo, rpanelTypes,
				0);

		File listOfStages = new File("/gpfs/scratch/pr1ees00/pr1ees14/GCAT/SHAPEIT_IMPUTE/list.txt");
		try {
			listOfStages.createNewFile();
		} catch (IOException e) {
			e.printStackTrace();
		}
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(listOfStages))) {
			for (String str : listOfCommands) {
				writer.write(str);
				writer.newLine();
				writer.newLine();
			}
			// Close the file with the list of commands...
			writer.flush();
		}
	}

	private static void newSample(String[] args) throws IOException, InterruptedException, Exception {
		String sampleFile = args[0];
		String phasingSampleFile = args[1];
		String phasingNewSampleFile = args[2];
		String responseVar = args[3];
		String covariables = args[4];

		GuidanceImpl.newSample(sampleFile, phasingSampleFile, phasingNewSampleFile, responseVar, covariables, "");
	}

	private static void imputeWithImpute(String[] args) {
		String gmapFile = args[0];
		String knownHapFile = args[1];
		String legendFile = args[2];
		String phasingHapsFile = args[3];
		String phasingSampleFile = args[4];
		String lim1S = args[5];
		String lim2S = args[6];
		String pairsFile = args[7];
		String imputeFile = args[8];
		String imputeFileInfo = args[9];
		String imputeFileSummary = args[10];
		String imputeFileWarnings = args[11];
		String theChromo = args[12];
		String sex = args[13];
		try {
			System.err.println("[Guidance] Exception trying the execution of imputeWithImpute task");
			GuidanceImpl.imputeWithImpute(gmapFile, knownHapFile, legendFile, phasingHapsFile, phasingSampleFile, lim1S,
					lim2S, pairsFile, imputeFile, imputeFileInfo, imputeFileSummary, imputeFileWarnings, theChromo, sex,
					"");
		} catch (GuidanceTaskException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void filterByInfo(String[] args) {
		String imputationTool = args[0];
		String imputeFileInfo = args[1];
		String inclusionRsIdFile = args[2];
		String infoThresholdS = args[3];
		String mafThresholdS = args[4];
		try {
			GuidanceImpl.filterByInfo(imputationTool, imputeFileInfo, inclusionRsIdFile, infoThresholdS, mafThresholdS,
					"");
		} catch (GuidanceTaskException e) {
			// TODO Auto-generated catch block
			System.err.println("[Guidance] Exception trying the execution of filterByInfo task");
			e.printStackTrace();
		}
	}

	public static void main(String[] args) throws Exception {
		System.out.println(Arrays.toString(args));
		int option = 5;
		switch (option) {
		case 0:
			createSplitedFilesWraper(args);
			break;
		case 1:
			convertFromBedToBedWrapper(args);
		case 2:
			splitChr23Wrapper(args);
		case 3:
			createRsIdListWrapper(args);
		case 8:
			collectSummaryWrapper(args);
			break;
		case 9:
			filterByAllWrapper(args);
			break;
		case 10:
			mergeOfChunksWrapper(args);
			break;
		case 11:
			printPaths(args);
			break;
		case 12:
			newSample(args);
			break;
		case 13:
			imputeWithImpute(args);
			break;
		case 14:
			filterByInfo(args);
			break;
		}
	}

}
