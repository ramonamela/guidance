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

package guidance.utils;

public class ChromoInfo {

    private static final int MAX_NUMBER_OF_CHROMOSOMES = 23;

    private static final int[] MIN_SIZE;
    private static final int[] MAX_SIZE;

    static {
        // The minimum size is always 1
        MIN_SIZE = new int[MAX_NUMBER_OF_CHROMOSOMES];
        for (int kk = 0; kk < MAX_NUMBER_OF_CHROMOSOMES; kk++) {
            MIN_SIZE[kk] = 1;
        }

        // The maximum size (maxSize) of each chromosome is well known, but it does not have
        // a rule to automatically determine it, therefore, we put the values individually.
        MAX_SIZE = new int[] { 252_000_000, // 0
                246_000_000, // 1
                201_000_000, // 2
                192_000_000, // 3
                183_000_000, // 4
                174_000_000, // 5
                162_000_000, // 6
                147_000_000, // 7
                144_000_000, // 8
                138_000_000, // 9
                138_000_000, // 10
                135_000_000, // 11
                117_000_000, // 12
                108_000_000, // 13
                105_000_000, // 14
                93_000_000, // 15
                84_000_000, // 16
                81_000_000, // 17
                60_000_000, // 18
                66_000_000, // 19
                51_000_000, // 20
                54_000_000, // 21
                156_000_000 // 22
        };
    }


    /**
     * Class constructor
     * 
     */
    public ChromoInfo() {
        // Nothing to do since all attributes are static
    }

    /**
     * Method to access the minSize of a chromosome
     * 
     * @param chromoNumber
     * @return
     */
    public int getMinSize(int chromoNumber) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromoNumber);

        return MIN_SIZE[chromoNumber - 1];
    }

    /**
     * Method to access the maxSize of a chromosome
     * 
     * @param chromoNumber
     * @return
     */
    public int getMaxSize(int chromoNumber) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromoNumber);

        return MAX_SIZE[chromoNumber - 1];
    }

    /**
     * Method to access biomart information
     * 
     * @param chromoNumber
     */
    public void printChromoInfo(int chromoNumber) {
        // Check that chromo index is within the bounds
        checkChromoIndex(chromoNumber);

        System.out.println("Gen file information for the chromosome " + chromoNumber);
        System.out.println("Max size     : " + MAX_SIZE[chromoNumber - 1]);
    }

    private void checkChromoIndex(int chromo) {
        if ((chromo < 1) || (chromo > MAX_NUMBER_OF_CHROMOSOMES)) {
            System.err.println("[chromoInfo] Errro, chromo " + chromo + "does not exist");
            System.exit(1);
        }
    }

}
