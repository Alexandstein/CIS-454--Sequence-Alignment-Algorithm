THOUGHTS:
	For this project, I went with an Object-oriented approach using the Python standard
library as my framework. I went with a design pattern for the class resembling a single
machine that you feed an input and press a "button" to generate an output. (class Aligner)
As per the suggestion in the homework instructions, I went with a specific class for each 
of he dynamic programming cells (class Cell) which contain the alignment score and 
information on where the parent cell is located to that traceback is possible.
	Unfortunately it doesn't seem like my project is fully functional. I tested the code
on the genetic sequences in the class notes and got good results with the scoring matrix,
but the traceback system seems to be non-functional as it doesn't seem to be able to trace
back correctly to the root cell at (0,0) despite producing the correct scoring.
	Aspects I found interesting on the project are how tight and succint dynamic
programming code is able to be, and that the challenge of the project was deciding exactly
how to write those few lines of code that the rest of the code centers on.

USAGE:
	The program is executed by invoking the interpreter on the file. `python align.py`
You will be prompted for two sequences to be compared, and be asked whether or not you 
would like to change the gap and mismatch penalties.
	The output is an n*m matrix of scores, with the dimensions being the lengths of the 
sequences + 1. Another output is a string of m's and g's, with m's meaning match/mismatch
and the g's meaning gap. As stated above, traceback is broken and hasn't been able to be 
fixed so only m's show up.

