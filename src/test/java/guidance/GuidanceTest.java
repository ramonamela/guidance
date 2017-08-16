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

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class GuidanceTest {

    @Test
    public void nameFromFile() {
        String filePath = "/a/full/path/";
        String fileName = "fileName";
        String fileExtension = ".gz";

        String fullFile = filePath + fileName + fileExtension;

        String genFileName = fullFile.substring(0, fullFile.length() - 3);
        String genFileExt = fullFile.substring(fullFile.length() - 3);

        assertEquals(filePath + fileName, genFileName);
        assertEquals(fileExtension, genFileExt);
    }

    @Test
    public void testTabFormat() {
        StringBuilder sb = new StringBuilder();
        sb.append("field1");
        sb.append("\t");
        sb.append("field1");
        for (String field : sb.toString().split("\t")) {
            assertEquals("field1", field);
        }

        String header = "field1\tfield1";
        for (String field : header.split("\t")) {
            assertEquals("field1", field);
        }
    }

}
