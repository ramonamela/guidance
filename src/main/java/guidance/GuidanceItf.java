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

import es.bsc.compss.types.annotations.Constraints;
import es.bsc.compss.types.annotations.task.Method;
import es.bsc.compss.types.annotations.Parameter;
import es.bsc.compss.types.annotations.parameter.Direction;
import es.bsc.compss.types.annotations.parameter.Type;

public interface GuidanceItf {

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void createSplitedFiles(@Parameter(type = Type.FILE, direction = Direction.IN) String gmapFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String gmapFileChr,
			@Parameter(type = Type.STRING, direction = Direction.IN) String chromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void convertFromBedToBed(@Parameter(type = Type.FILE, direction = Direction.IN) String bedFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bimFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String famFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String newBedFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String newBimFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String newFamFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String logFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String chromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void splitChr23(@Parameter(type = Type.STRING, direction = Direction.IN) String myPrefix,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bedFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bimFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String famFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String bedChr23File,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String bimChr23File,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String famChr23File,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String logFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String sex,
			@Parameter(type = Type.STRING, direction = Direction.IN) String theChromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void convertFromBedToPed(@Parameter(type = Type.STRING, direction = Direction.IN) String bedPrefix,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bedFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bimFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String famFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String pedFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String mapFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String logFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String chromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void convertFromPedToGen(@Parameter(type = Type.FILE, direction = Direction.IN) String pedFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String mapFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String genFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String sampleFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String logFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void createRsIdList(@Parameter(type = Type.FILE, direction = Direction.IN) String genOrBimFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String exclCgatFlag,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String pairsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String inputFormat,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${gtoolsMem}")
	void gtoolS(@Parameter(type = Type.FILE, direction = Direction.IN) String newGenFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String modSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String gtoolGenFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String gtoolSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String sampleExclFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String snpWtccFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String gtoolLogFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${qctoolSMem}")
	void qctoolS(@Parameter(type = Type.FILE, direction = Direction.IN) String imputeFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String inclusionRsIdFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafThresholdS,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredFileGz,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredLogFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${createListOfExcludedSnpsMem}")
	void createListOfExcludedSnps(@Parameter(type = Type.FILE, direction = Direction.IN) String phasingHapsFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String excludedSnpsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String exclCgatFlag,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "${phasingBedCU}", memorySize = "${phasingBedMem}")
	void phasingBed(@Parameter(type = Type.STRING, direction = Direction.IN) String chromo,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bedFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String bimFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String famFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String gmapFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phasingHapsFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phasingSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phasingLogFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String myPrefix,
			@Parameter(type = Type.STRING, direction = Direction.IN) String sex,
			@Parameter(type = Type.STRING, direction = Direction.IN) String phasingTool,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "${phasingCU}", memorySize = "${phasingMem}")
	void phasing(@Parameter(type = Type.STRING, direction = Direction.IN) String chromo,
			@Parameter(type = Type.FILE, direction = Direction.IN) String inputGenFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String inputSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String gmapFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phasingHapsFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phasingSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phasingLogFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String phasingTool,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void newSample(@Parameter(type = Type.FILE, direction = Direction.IN) String sampleFile,
			@Parameter(type = Type.FILE, direction = Direction.INOUT) String phasingSampleFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String responseVar,
			@Parameter(type = Type.STRING, direction = Direction.IN) String covariables,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${filterHaplotypesMem}")
	void filterHaplotypes(@Parameter(type = Type.FILE, direction = Direction.IN) String phasingHapsFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String phasingSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String excludedSnpsFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredLogFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredHapsVcfFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void samtoolsBgzip(@Parameter(type = Type.FILE, direction = Direction.INOUT) String input,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${imputeWithImputeMem}")
	void imputeWithImpute(@Parameter(type = Type.FILE, direction = Direction.IN) String gmapFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String knownHapFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String legendFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String shapeitHapsFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String shapeitSampleFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String lim1S,
			@Parameter(type = Type.STRING, direction = Direction.IN) String lim2S,
			@Parameter(type = Type.FILE, direction = Direction.IN) String pairsFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputeFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputeFileInfo,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputeFileSummary,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputeFileWarnings,
			@Parameter(type = Type.STRING, direction = Direction.IN) String theChromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${imputeWithMinimacMem}")
	void imputeWithMinimac(@Parameter(type = Type.FILE, direction = Direction.IN) String knownHapFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredHapsFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredListOfSnpsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String imputedMMFileName,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputedMMInfoFile,
			// @Parameter(type = Type.FILE, direction = Direction.OUT) String
			// imputedMMDraftFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputedMMErateFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputedMMRecFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputedMMDoseFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String imputedMMLogFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String theChromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String lim1S,
			@Parameter(type = Type.STRING, direction = Direction.IN) String lim2S,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${filterByAllMem}")
	void filterByAll(@Parameter(type = Type.FILE, direction = Direction.IN) String inputFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String outputFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String outputCondensedFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpanelName,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String infoThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweCohortThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweCasesThreshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafControlsThreshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${filterByInfoMem}")
	void filterByInfo(@Parameter(type = Type.FILE, direction = Direction.IN) String imputeFileInfo,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String threshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${jointFilteredByAllFilesMem}")
	void jointFilteredByAllFiles(@Parameter(type = Type.FILE, direction = Direction.IN) String filteredByAllA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredByAllB,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredByAllC,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpanelName,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpanelFlag,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${jointCondensedFilesMem}")
	void jointCondensedFiles(@Parameter(type = Type.FILE, direction = Direction.IN) String inputAFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String inputBFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String outputFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void collectSummary(@Parameter(type = Type.STRING, direction = Direction.IN) String chr,
			@Parameter(type = Type.FILE, direction = Direction.IN) String firstImputeFileInfo,
			@Parameter(type = Type.FILE, direction = Direction.IN) String snptestOutFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String reduceFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String infoThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweCohortThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweCasesThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweControlsThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${generateQQManhattanPlotsMem}")
	void generateQQManhattanPlots(@Parameter(type = Type.FILE, direction = Direction.IN) String lastCondensedFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String qqPlotFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String manhattanPlotFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String qqPlotTiffFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String manhattanPlotTiffFile,
			// @Parameter(type = Type.FILE, direction = Direction.OUT) String
			// correctedPvaluesFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${filterByMafMem}")
	void filterByMaf(@Parameter(type = Type.FILE, direction = Direction.IN) String summaryFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String filteredSummaryFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${snptestMem}")
	void snptest(@Parameter(type = Type.FILE, direction = Direction.IN) String mergedGenFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String mergedSampleFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String snptestOutFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String snptestLogFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String responseVar,
			@Parameter(type = Type.STRING, direction = Direction.IN) String covariables,
			@Parameter(type = Type.STRING, direction = Direction.IN) String theChromo,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${initPhenoMatrixMem}")
	void initPhenoMatrix(@Parameter(type = Type.FILE, direction = Direction.IN) String topHitsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String ttName,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpName,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phenomeFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${addPhenoMatrixMem}")
	void addToPhenoMatrix(@Parameter(type = Type.FILE, direction = Direction.IN) String phenomeAFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String topHitsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String ttName,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpName,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phenomeBFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${filloutPhenoMatrixMem}")
	void filloutPhenoMatrix(@Parameter(type = Type.FILE, direction = Direction.IN) String phenomeAFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredByAllFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredByAllXFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String thereisX,
			@Parameter(type = Type.STRING, direction = Direction.IN) String ttName,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpName,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phenomeBFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${finalizePhenoMatrixMem}")
	void finalizePhenoMatrix(@Parameter(type = Type.FILE, direction = Direction.IN) String phenomeAFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String phenomeBFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String ttName,
			@Parameter(type = Type.STRING, direction = Direction.IN) String rpName,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String phenomeCFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${mergeTwoChunksMem}")
	void mergeTwoChunks(@Parameter(type = Type.FILE, direction = Direction.IN) String reduceFileA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String reduceFileB,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String reduceFileC,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${mergeTwoChunksInTheFirstMem}")
	void mergeTwoChunksInTheFirst(@Parameter(type = Type.FILE, direction = Direction.INOUT) String reduceFileA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String reduceFileB,
			@Parameter(type = Type.STRING, direction = Direction.IN) String chrS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String type,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${combinePanelsMem}")
	void combinePanels(@Parameter(type = Type.FILE, direction = Direction.IN) String resultsFileA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String resultsFileB,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String resultsFileC,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void combinePanelsComplex(@Parameter(type = Type.FILE, direction = Direction.INOUT) String resultsFileA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String resultsFileB,
			@Parameter(type = Type.INT, direction = Direction.IN) int lim1,
			@Parameter(type = Type.INT, direction = Direction.IN) int lim2,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${combineCondensedFilesMem}")
	void combineCondensedFiles(@Parameter(type = Type.FILE, direction = Direction.IN) String filteredA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String filteredX,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String combinedCondensedFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String infoThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweCohortThresholdS,
			@Parameter(type = Type.STRING, direction = Direction.IN) String hweCasesThreshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String mafControlsThreshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${combinePanelsComplex1Mem}")
	void combinePanelsComplex1(@Parameter(type = Type.FILE, direction = Direction.IN) String resultsFileA,
			@Parameter(type = Type.FILE, direction = Direction.IN) String resultsFileB,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String resultsFileC,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${generateTopHitsMem}")
	void generateTopHits(@Parameter(type = Type.FILE, direction = Direction.IN) String resultsFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String outputTopHitsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String pvaThreshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "${generateTopHitsAllMem}")
	void generateTopHitsAll(@Parameter(type = Type.FILE, direction = Direction.IN) String resultsAFile,
			@Parameter(type = Type.FILE, direction = Direction.IN) String resultsBFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String outputTopHitsFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String pvaThreshold,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void gzipgunzip(@Parameter(type = Type.FILE, direction = Direction.IN) String inputFile,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String outFile,
			@Parameter(type = Type.STRING, direction = Direction.IN) String compressType,
			@Parameter(type = Type.STRING, direction = Direction.IN) String cmdToStore);

	@Method(declaringClass = "guidance.GuidanceImpl")
	@Constraints(computingUnits = "1", memorySize = "1.0f")
	void copyFile(@Parameter(type = Type.FILE, direction = Direction.IN) String fileA,
			@Parameter(type = Type.FILE, direction = Direction.OUT) String fileB);

}
