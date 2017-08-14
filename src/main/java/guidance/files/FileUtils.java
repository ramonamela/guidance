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
        tmpOutDir = myOutDir + "/" + mixedCohort;

        System.out.println("[Guidance] Creating " + tmpOutDir);
        FileUtils.createDir(tmpOutDir);

        // We create the second directory: the REFPANEL directory.
        // It can be => common: for files from the beginning to phasing
        // => PANEL : for files for imputation for each panel

        // Now I create the directories for common and combined
        String tmpPanelDir = tmpOutDir + "/" + "common";
        FileUtils.createDir(tmpPanelDir);

        // Next level: Create mixed directories.
        String mixedOutDir = tmpPanelDir + "/mixed";
        System.out.println("[Guidance] Creating " + mixedOutDir);
        FileUtils.createDir(mixedOutDir);
        for (int i = startChr; i <= endChr; i++) {
            String tmpChrDir = mixedOutDir + "/Chr_" + i;
            FileUtils.createDir(tmpChrDir);
        }

        // Now I create the directories for imputeOutDir
        for (int j = 0; j < refPanels.size(); j++) {
            String rPanel = refPanels.get(j);
            tmpPanelDir = tmpOutDir + "/" + rPanel;
            FileUtils.createDir(tmpPanelDir);

            // Next level: Create mixed directories.
            mixedOutDir = tmpPanelDir + "/mixed";
            System.out.println("[Guidance] Creating " + mixedOutDir);
            FileUtils.createDir(mixedOutDir);
            for (int i = startChr; i <= endChr; i++) {
                String tmpChrDir = mixedOutDir + "/Chr_" + i;
                FileUtils.createDir(tmpChrDir);
            }
        }

        // Then we create the directories for the Association files
        // This directory is created by default in myOutDir/association.
        String tmpMyOutDir1 = myOutDir + "/associations";
        FileUtils.createDir(tmpMyOutDir1);
        String testTypeName = null;
        int numberOfTestTypes = parsingArgs.getNumberOfTestTypeName();
        for (int kk = 0; kk < numberOfTestTypes; kk++) {
            testTypeName = parsingArgs.getTestTypeName(kk);
            String tmpMyOutDir = myOutDir + "/associations/" + testTypeName;
            FileUtils.createDir(tmpMyOutDir);

            for (int j = 0; j < refPanels.size(); j++) {
                String rPanel = refPanels.get(j);
                String assocDir = null;
                assocDir = tmpMyOutDir + "/" + mixedCohort + "_for_" + rPanel;
                System.out.println("[Guidance] Creating " + assocDir);
                FileUtils.createDir(assocDir);

                for (int i = startChr; i <= endChr; i++) {
                    String tmpChrDir = assocDir + "/Chr_" + i;
                    FileUtils.createDir(tmpChrDir);
                }

                String summaryDir = assocDir + "/summary";
                System.out.println("[Guidance] Creating " + summaryDir);
                FileUtils.createDir(summaryDir);
            }

            // Finally, if there are several panel and if combinedRefPanel = 1
            // then we have to create the output directory for the combined ref_panel
            boolean refPanelCombine = parsingArgs.getRefPanelCombine();

            if (refPanelCombine == true) {
                // We take the first refPanel name
                String rPanel = refPanels.get(0);
                String combinedRefPanel = rPanel;
                for (int j = 1; j < refPanels.size(); j++) {
                    rPanel = refPanels.get(j);
                    combinedRefPanel = combinedRefPanel + "_" + rPanel;
                }

                String combinedAssocOutDir = null;
                combinedAssocOutDir = myOutDir + "/associations/" + testTypeName + "/" + mixedCohort + "_combined_panels_"
                        + combinedRefPanel;

                System.out.println("[Guidance] Creating " + combinedAssocOutDir);
                FileUtils.createDir(combinedAssocOutDir);

            } // End if(refPanelCombine == true)
        } // End for(kk=0;kk<numberOfTestTyes;kk++)

        // Now create the structure for the phenotype analysis results
        String phenomeAnalysisOutDir = myOutDir + "/associations/pheno_analysis/";
        System.out.println("[Guidance] Creating " + phenomeAnalysisOutDir);
        FileUtils.createDir(phenomeAnalysisOutDir);

        phenomeAnalysisOutDir = phenomeAnalysisOutDir + mixedCohort;
        // for(int tt=0; tt< numberOfTestTypes; tt++) {
        // testTypeName = parsingArgs.getTestTypeName(tt);
        // phenomeAnalysisOutDir = phenomeAnalysisOutDir + "_" + testTypeName;
        // }
        System.out.println("[Guidance] Creating " + phenomeAnalysisOutDir);
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
     * @throws IOException
     */
    public static void createEmptyFile(String fileName, String moduleName) throws IOException {
        File fa = new File(fileName);
        if (!fa.exists()) {
            System.out.println(moduleName + " The file " + fileName + " does not exist, then we create an empty file");
            fa.createNewFile();
        }
    }

    /**
     * Method to zip a file
     * 
     * @param sourceFilePath
     * @param destZipFilePath
     */
    public static void gzipFile(String sourceFilePath, String destZipFilePath) {
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
    public static void gunzipFile(String compressedFile, String decompressedFile) {
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

}
