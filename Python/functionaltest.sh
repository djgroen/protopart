python protopart.py -i ../Input\ Files/four_cube.gmy -o output.txt --all 4 > tmp_out 
diff sample_out.txt tmp_out
rm tmp_out
