# This file contains the functional tests for Protopart.

python protopart.py -i ../Input\ Files/four_cube.gmy -o functionaltest_data/output.txt --all 4 > functionaltest_data/tmp_out 
diff functionaltest_data/sample_out.txt functionaltest_data/tmp_out
