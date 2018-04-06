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

    private static void collectSummaryWrapper(String[] args) {
        String chr = args[0];
        String firstImputeFileInfo = args[1];
        String snptestOutFile = args[2];
        String reduceFile = args[3];
        String mafThresholdS = args[4];
        String infoThresholdS = args[5];
        String hweCohortThresholdS = args[6];
        String hweCasesThresholdS = args[7];
        String hweControlsThresholdS = args[8];
        String cmdToStore = args[9];
        try {
            GuidanceImpl.collectSummary(chr, firstImputeFileInfo, snptestOutFile, reduceFile, mafThresholdS, infoThresholdS,
                    hweCohortThresholdS, hweCasesThresholdS, hweControlsThresholdS, cmdToStore);
        } catch (GuidanceTaskException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private static void filterByAllWrapper(String[] args) {
        String inputFile = args[0];
        String outputFile = args[1];
        String outputCondesedFile = args[2];
        String mafThresholdS = args[3];
        String infoThresholdS = args[4];
        String hweCohortThresholdS = args[5];
        String hweCasesThresholdS = args[6];
        String hweControlsThresholdS = args[7];
        String cmdToStore = "";
        try {
            GuidanceImpl.filterByAll(inputFile, outputFile, outputCondesedFile, mafThresholdS, infoThresholdS, hweCohortThresholdS,
                    hweCasesThresholdS, hweControlsThresholdS, cmdToStore);
        } catch (GuidanceTaskException e) {
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
            e.printStackTrace();
        }
    }
    
    private static void generateTopHits(String[] args) {
        String filteredFile = args[0];
        String filteredXFile = args[1];
        String topHitsResults = args[2];
        String pvaThrS = args[3];
        try {
            GuidanceImpl.generateTopHitsAll(filteredFile, filteredFile, topHitsResults, pvaThrS, "");
        } catch (GuidanceTaskException e) {
            e.printStackTrace();
        }
    }
    
    private static void printPaths(String[] args) throws IOException {
        
        ArrayList<String> listOfCommands = new ArrayList<>();
        
        ParseCmdLine parsingArgs = new ParseCmdLine(args);
        
        String outDir = parsingArgs.getOutDir();
        
        ArrayList<String> rpanelTypes = new ArrayList<>(parsingArgs.getRpanelTypes());
        
        // Create the names for Association files
        AssocFiles assocFilesInfo = new AssocFiles(parsingArgs, outDir, rpanelTypes);

        // Create the names for Merge files
        MergeFiles mergeFilesInfo = new MergeFiles(parsingArgs, outDir, rpanelTypes);
        
        CombinedPanelsFiles combinedPanelsFilesInfo = new CombinedPanelsFiles(parsingArgs, outDir, rpanelTypes);
        
        Guidance.makeCombinePanels(parsingArgs, assocFilesInfo, mergeFilesInfo, combinedPanelsFilesInfo, rpanelTypes, 0, listOfCommands);
        
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

    public static void main(String[] args) throws Exception {
        System.out.println(Arrays.toString(args));
        int option = 4;
        switch (option) {
        case 0:
        	collectSummaryWrapper(args);
        case 1:
        	filterByAllWrapper(args);
        case 2:
        	mergeOfChunksWrapper(args);
        case 3:
        	generateTopHits(args);
        case 4:
        	printPaths(args);
        }
    }
}
