package guidance.files;

import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Files;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import guidance.exceptions.GuidanceTaskException;
import guidance.utils.ParseCmdLine;


public class FileUtils {

    // Messages
    private static final String GUIDANCE_MSG_HEADER = "[Guidance]";
    private static final String GUIDANCE_MSG_CREATE_FILE = GUIDANCE_MSG_HEADER + " Creating ";

    // File paths
    private static final String COMMON_REL_PATH = "common";
    private static final String MIXED_REL_PATH = "mixed";
    private static final String CHR_REL_PATH = "Chr_";
    private static final String ASSOCIATIONS_REL_PATH = "associations";


    private FileUtils() {
        // Private constructor to avoid instantiation
    }

    /**
     * Method that creates all the output directory structure for the results.
     * 
     * @param parsingArgs
     * @param myOutDir
     * @param refPanels
     * @param startChr
     * @param endChr
     */
    public static void createDirStructure(ParseCmdLine parsingArgs, String myOutDir, List<String> refPanels, int startChr, int endChr) {
        String mixedCohort = null;
        // We create the first directory: the cohort directory.
        mixedCohort = parsingArgs.getCohort();

        String tmpOutDir = "";
        tmpOutDir = myOutDir + File.separator + mixedCohort;

        System.out.println("[Guidance] Creating " + tmpOutDir);
        FileUtils.createDir(tmpOutDir);

        // We create the second directory: the REFPANEL directory.
        // It can be => common: for files from the beginning to phasing
        // => PANEL : for files for imputation for each panel

        // Now I create the directories for common and combined
        String tmpPanelDir = tmpOutDir + File.separator + COMMON_REL_PATH;
        FileUtils.createDir(tmpPanelDir);

        // Next level: Create mixed directories.
        String mixedOutDir = tmpPanelDir + File.separator + MIXED_REL_PATH;
        System.out.println(GUIDANCE_MSG_CREATE_FILE + mixedOutDir);
        FileUtils.createDir(mixedOutDir);
        for (int i = startChr; i <= endChr; i++) {
            String tmpChrDir = mixedOutDir + File.separator + CHR_REL_PATH + i;
            FileUtils.createDir(tmpChrDir);
        }

        // Now I create the directories for imputeOutDir
        for (int j = 0; j < refPanels.size(); j++) {
            String rPanel = refPanels.get(j);
            tmpPanelDir = tmpOutDir + File.separator + rPanel;
            FileUtils.createDir(tmpPanelDir);

            // Next level: Create mixed directories.
            mixedOutDir = tmpPanelDir + File.separator + MIXED_REL_PATH;
            System.out.println(GUIDANCE_MSG_CREATE_FILE + mixedOutDir);
            FileUtils.createDir(mixedOutDir);
            for (int i = startChr; i <= endChr; i++) {
                String tmpChrDir = mixedOutDir + File.separator + CHR_REL_PATH + i;
                FileUtils.createDir(tmpChrDir);
            }
        }

        // Then we create the directories for the Association files
        // This directory is created by default in myOutDir/association.
        String tmpMyOutDir1 = myOutDir + File.separator + ASSOCIATIONS_REL_PATH;
        FileUtils.createDir(tmpMyOutDir1);
        String testTypeName = null;
        int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
        for (int kk = 0; kk < numberOfTestTypes; kk++) {
            testTypeName = parsingArgs.getTestTypeName(kk);
            String tmpMyOutDir = myOutDir + File.separator + ASSOCIATIONS_REL_PATH + File.separator + testTypeName;
            FileUtils.createDir(tmpMyOutDir);

            for (int j = 0; j < refPanels.size(); j++) {
                String rPanel = refPanels.get(j);
                String assocDir = null;
                assocDir = tmpMyOutDir + File.separator + mixedCohort + "_for_" + rPanel;
                System.out.println(GUIDANCE_MSG_CREATE_FILE + assocDir);
                FileUtils.createDir(assocDir);

                for (int i = startChr; i <= endChr; i++) {
                    String tmpChrDir = assocDir + "/Chr_" + i;
                    FileUtils.createDir(tmpChrDir);
                }

                String summaryDir = assocDir + File.separator + "summary";
                System.out.println(GUIDANCE_MSG_CREATE_FILE + summaryDir);
                FileUtils.createDir(summaryDir);
            }

            // Finally, if there are several panel and if combinedRefPanel = 1
            // then we have to create the output directory for the combined ref_panel
            boolean refPanelCombine = parsingArgs.getRefPanelCombine();

            if (refPanelCombine) {
                // We take the first refPanel name
                String rPanel = refPanels.get(0);
                String combinedRefPanel = rPanel;
                for (int j = 1; j < refPanels.size(); j++) {
                    rPanel = refPanels.get(j);
                    combinedRefPanel = combinedRefPanel + "_" + rPanel;
                }

                String combinedAssocOutDir = null;
                combinedAssocOutDir = myOutDir + File.separator + ASSOCIATIONS_REL_PATH + File.separator + testTypeName + File.separator
                        + mixedCohort + "_combined_panels_" + combinedRefPanel;

                System.out.println(GUIDANCE_MSG_CREATE_FILE + combinedAssocOutDir);
                FileUtils.createDir(combinedAssocOutDir);

            } // End if
        } // End for kk < numberOfTestTyes

        // Now create the structure for the phenotype analysis results
        String phenomeAnalysisOutDir = myOutDir + File.separator + ASSOCIATIONS_REL_PATH + File.separator + "pheno_analysis"
                + File.separator;
        System.out.println(GUIDANCE_MSG_CREATE_FILE + phenomeAnalysisOutDir);
        FileUtils.createDir(phenomeAnalysisOutDir);

        phenomeAnalysisOutDir = phenomeAnalysisOutDir + mixedCohort;
        // for(int tt=0; tt< numberOfTestTypes; tt++) {
        // testTypeName = parsingArgs.getTestTypeName(tt);
        // phenomeAnalysisOutDir = phenomeAnalysisOutDir + "_" + testTypeName;
        // }
        System.out.println(GUIDANCE_MSG_CREATE_FILE + phenomeAnalysisOutDir);
        FileUtils.createDir(phenomeAnalysisOutDir);
    }

    /**
     * Method that creates a single directory
     * 
     * @param outputDir
     */
    private static void createDir(String outputDir) {
        File tmpOutDir = new File(outputDir);
        boolean existOutDir = tmpOutDir.exists();
        if (!existOutDir) {
            boolean successExistOutDir = (new File(outputDir)).mkdir();
            if (!successExistOutDir) {
                System.err.println("[createDir] Error, cannot create " + tmpOutDir + " directory");
                System.exit(1);
            }
        }
    }

    /**
     * Method to copy file source to dest
     * 
     * @param source
     * @param dest
     * @throws IOException
     */
    public static void copyFile(File source, File dest) throws GuidanceTaskException {
        try {
            Files.copy(source.toPath(), dest.toPath(), REPLACE_EXISTING);
        } catch (IOException ioe) {
            throw new GuidanceTaskException(ioe);
        }
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
     * Creates an empty file if it does not exist
     * 
     * @param fileName
     * @param moduleName
     * @return true if the file was created or already existant, false otherwise
     * @throws IOException
     */
    public static boolean createEmptyFile(String fileName, String moduleName) throws IOException {
        File fa = new File(fileName);
        if (!fa.exists()) {
            System.out.println(moduleName + " The file " + fileName + " does not exist, then we create an empty file");
            return fa.createNewFile();
        }

        // The file already exists
        return true;
    }

    /**
     * Method to zip a file
     * 
     * @param sourceFilePath
     * @param destZipFilePath
     */
    public static void gzipFile(String sourceFilePath, String destZipFilePath) {
        byte[] buffer = new byte[1024];
        try (FileOutputStream fileOutputStream = new FileOutputStream(destZipFilePath);
                GZIPOutputStream gzipOuputStream = new GZIPOutputStream(fileOutputStream);
                FileInputStream fileInput = new FileInputStream(sourceFilePath)) {

            int bytesRead;
            while ((bytesRead = fileInput.read(buffer)) > 0) {
                gzipOuputStream.write(buffer, 0, bytesRead);
            }

            gzipOuputStream.finish();
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
    public static void gunzipFile(String compressedFile, String decompressedFile) {
        byte[] buffer = new byte[1024];

        try (FileInputStream fileIn = new FileInputStream(compressedFile);
                GZIPInputStream gZIPInputStream = new GZIPInputStream(fileIn);
                FileOutputStream fileOutputStream = new FileOutputStream(decompressedFile)) {

            int bytesRead;
            while ((bytesRead = gZIPInputStream.read(buffer)) > 0) {
                fileOutputStream.write(buffer, 0, bytesRead);
            }

            // System.out.println("The file was decompressed successfully!");
        } catch (IOException ex) {
            System.err.println(ex.getMessage());
            ex.printStackTrace();
        }
    }

}
