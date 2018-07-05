#!/bin/bash

### PHASE 1 ###

export phasingMem="50.0"
export phasingCU="48"

export phasingBedMem="50.0"
export phasingBedCU="48"


### PHASE 2 ###

export qctoolMem="16.0"

export qctoolSMem="16.0"

export gtoolsMem="6.0"

export samtoolsBgzipMem="6.0"

export imputeWithImputeLowMem="12.0"
export imputeWithImputeMediumMem="12.0"
export imputeWithImputeHighMem="12.0"

export imputeWithMinimacLowMem="8.0"
export imputeWithMinimacMediumMem="8.0"
export imputeWithMinimacHighMem="8.0"

export filterByInfoImputeMem="12.0"

export filterByInfoMinimacMem="24.0"

### PHASE 3 ###

export createListOfExcludedSnpsMem="6.0"

export filterHaplotypesMem="6.0"

export filterByAllMem="16.0"

export jointFilteredByAllFilesMem="12.0"

export jointCondensedFilesMem="12.0"

export generateQQManhattanPlotsMem="90.0"

export generateTopHitsAllMem="90.0"

export generateTopHitsMem="90.0"

export filterByMafMem="12.0"

export snptestMem="12.0"

export initPhenoMatrixMem="12.0"

export addPhenoMatrixMem="12.0"

export filloutPhenoMatrixMem="12.0"

export finalizePhenoMatrixMem="12.0"

export mergeTwoChunksMem="12.0"

export mergeTwoChunksInTheFirstMem="12.0"

export combinePanelsMem="24.0"

export combineCondensedFilesMem="24.0"

export combinePanelsComplex1Mem="24.0"
