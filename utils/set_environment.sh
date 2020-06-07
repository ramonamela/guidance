#!/bin/bash

### PHASE 1 ###

export phasingMem="4.0"
export phasingCU="2"

export phasingBedMem="4.0"
export phasingBedCU="2"


### PHASE 2 ###

export qctoolMem="1.0"

export qctoolSMem="1.0"

export gtoolsMem="1.0"

export samtoolsBgzipMem="1.0"

export imputeWithImputeLowMem="1.0"
export imputeWithImputeMediumMem="2.0"
export imputeWithImputeHighMem="4.0"

export imputeWithMinimacLowMem="1.0"
export imputeWithMinimacMediumMem="2.0"
export imputeWithMinimacHighMem="4.0"

export filterByInfoImputeMem="1.0"

export filterByInfoMinimacMem="1.0"

### PHASE 3 ###

export createListOfExcludedSnpsMem="1.0"

export filterHaplotypesMem="1.0"

export filterByAllMem="1.0"

export jointFilteredByAllFilesMem="1.0"

export jointCondensedFilesMem="1.0"

export generateTopHitsAllMem="1.0"

export generateTopHitsMem="1.0"

export filterByMafMem="1.0"

export snptestMem="1.0"

export mergeTwoChunksMem="1.0"

export mergeTwoChunksInTheFirstMem="1.0"

export combinePanelsMem="1.0"

export combineCondensedFilesMem="1.0"

export combinePanelsComplex1Mem="1.0"

### PHASE 4 ###

export generateCondensedTopHitsCU="2"

export generateCondensedTopHitsMem="4.0"

export generateQQManhattanPlotsCU="2"

export generateQQManhattanPlotsMem="8.0"

export phenoMergeMem="6.0"
