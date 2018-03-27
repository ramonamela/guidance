package guidance;

import java.util.Arrays;

import guidance.GuidanceImpl;
import guidance.exceptions.GuidanceTaskException;


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
        String theChromo = args[3];
        String type = args[4];
        String cmdToStore = "";
        try {
            GuidanceImpl.mergeTwoChunks(reduceA, reduceB, reduceC, theChromo, type, cmdToStore);
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

    public static void main(String[] args) throws Exception {
        System.out.println(Arrays.toString(args));
        // collectSummaryWrapper(args);
        //filterByAllWrapper(args);
        //mergeOfChunksWrapper(args);
        generateTopHits(args);
    }
}
