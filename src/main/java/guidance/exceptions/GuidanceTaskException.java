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

package guidance.exceptions;

/**
 * Class to represent exceptions while checking environment variables
 *
 */
public class GuidanceTaskException extends Exception {

    /**
     * Exception serial version ID
     */
    private static final long serialVersionUID = 1L;


    /**
     * Default Environment exception
     * 
     */
    public GuidanceTaskException() {
        super();
    }

    /**
     * Environment exception with error message
     * 
     * @param msg
     */
    public GuidanceTaskException(String msg) {
        super(msg);
    }

    /**
     * Environment exception with nested exception
     * 
     * @param e
     */
    public GuidanceTaskException(Exception e) {
        super(e);
    }

    /**
     * Environment exception with error message and nested exception
     * 
     * @param msg
     * @param e
     */
    public GuidanceTaskException(String msg, Exception e) {
        super(msg, e);
    }

}
