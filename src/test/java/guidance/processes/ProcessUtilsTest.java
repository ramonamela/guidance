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

package guidance.processes;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import guidance.processes.ProcessUtils;


public class ProcessUtilsTest {

    private static final String TEST_OUTPUT_FILE_1 = "test1.out";
    private static final String TEST_ERROR_FILE_1 = "test1.err";
    private static final String[] TEST_FILES = new String[] { TEST_OUTPUT_FILE_1, TEST_ERROR_FILE_1 };


    @BeforeClass
    public static void runBeforeClass() {
        deleteTestFiles();
    }

    @AfterClass
    public static void runAfterClass() {
        deleteTestFiles();
    }

    @Before
    public void setUp() {
        deleteTestFiles();
    }

    @After
    public void tearDown() {
        deleteTestFiles();
    }

    private static void deleteTestFiles() {
        for (String fileName : TEST_FILES) {
            File file = new File(fileName);
            if (file.exists()) {
                file.delete();
            }
        }
    }

    @Test
    public void binaryExitValue() throws IOException {
        String cmd = "ls";
        int exitValue = ProcessUtils.execute(cmd, TEST_OUTPUT_FILE_1, TEST_ERROR_FILE_1);
        assertEquals(0, exitValue);
    }

    @Test
    public void binaryOutput() throws IOException {
        String cmd = "echo hello";
        int exitValue = ProcessUtils.execute(cmd, TEST_OUTPUT_FILE_1, TEST_ERROR_FILE_1);

        assertEquals(0, exitValue);
        try (BufferedReader br = new BufferedReader(new FileReader(TEST_OUTPUT_FILE_1))) {
            String line = br.readLine();
            assertEquals("hello", line);
        }
    }

}
