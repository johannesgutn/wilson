This is the code used to calculate the Wilson line correlators introduced in https://arxiv.org/abs/2107.02542.
For more details check out the paper.

Running run.py saves four csv files. One is for the time, while the three others are the 
Wilson line correlators as a function of time, with various degrees of approximation.

configs.py contains the parameters used, and can be changed to whatever is suitable.
From here one can also edit the file name.
The number and direction of the Wilson lines can be changed by editing r,R,n1 and n2.


It needs the modules numpy and scipy to run
