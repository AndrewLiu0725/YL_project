## Tutorial
1. All the functions are packed into Doublet_Functions.py. This python file uses C to speedup, so compile the cforDoublet_Functions.c to cforDoublet_Functions.so first by the command line "gcc -fPIC -shared -o cforDoublet_Functions.so cforDoublet_Functions.c".
2. Use the Preprocessing.py to preprocess the data to improve the time complexity to analyze data.