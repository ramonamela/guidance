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

package guidance.files;

/**
 * @brief GenericClass. This clase includes 5 attributes
 * @author Friman Sanchez
 * @date 2014-10-15
 */
public class GenericFile {

    private final String dir;
    private final String name;
    private final String fullName;
    private final String generatedBy;

    private String finalStatus;


    /**
     * A first constructor for GenericFile class
     * 
     * @param myDir
     * @param myName
     * @param myFinalStatus
     */
    public GenericFile(String myDir, String myName, String myFinalStatus) {
        this.dir = myDir;
        this.name = myName;
        this.fullName = myDir + "/" + myName;

        this.generatedBy = null;

        this.finalStatus = myFinalStatus;
    }

    /**
     * A second constructor for GenericFile class
     * 
     * @param myDir
     * @param myName
     * @param myFinalStatus
     * @param myGeneratedBy
     */
    public GenericFile(String myDir, String myName, String myFinalStatus, String myGeneratedBy) {
        this.dir = myDir;
        this.name = myName;
        this.fullName = myDir + "/" + myName;

        this.generatedBy = myGeneratedBy;

        this.finalStatus = myFinalStatus;
    }

    /**
     * A method to set finalStatus field
     * 
     * @param myFinalStatus
     */
    public void setFinalStatus(String myFinalStatus) {
        this.finalStatus = myFinalStatus;
    }

    /**
     * A method to get dir field
     * 
     * @return
     */
    public String getDir() {
        return this.dir;
    }

    /**
     * A method to get Name field
     * 
     * @return
     */
    public String getName() {
        return this.name;
    }

    /**
     * A method to get finalStatus field
     * 
     * @return
     */
    public String getFinalStatus() {
        return this.finalStatus;
    }

    /**
     * A method to get fullName field
     * 
     * @return
     */
    public String getFullName() {
        return this.fullName;
    }

    /**
     * A method to get generatedBy field
     * 
     * @return
     */
    public String getGeneratedBy() {
        return this.generatedBy;
    }

}
