#from Bio.Seq import Seq
#from Bio import SeqIO
#from Bio.Alphabet.IUPAC import unambiguous_dna
#from Bio.SeqUtils import GC
import sys, getopt

class Cell:	
	'''
	Class to hold Dynamic Programming table data for Aligner
	Members:
	    (int 2-tuple) parent:
	        Points to parent cell's coordinates relative to current cell.
	    (int) value:
	        Alignment score of the cell.
	'''
	
	def __repr__(self):
		output = "Cell: %d| (%d,%d)" % (self.value, self.parent[0],self.parent[1])
		return output
	
	def __int__(self):
		return self.value
	
	def __init__(self, value=-1, parent=(0,0)):
		'''
		Constructor
		Args:
		    (int) value:
		        Score of cell
		    (int 2-tuple) parent:
		        Parent cell coords relative to itself. (0,0) means that the Cell is a
		        root and has no parents.
		'''
		self.value  = value
		self.parent = parent
	
	__str__ = __repr__

class Aligner:
	'''
	Class used to align sequences and get their relative similarities
	Members:
		(string 2-tuple) sequences:
			Genetic sequences to align.
		(n*n list of Cells) cells
			Cells to hold dynamic programming data structure. Dimensions are 
			length of the first sequence by length of the second sequence.
		(dict: string->int) penalties:
			Penalty values. Defaults are ['gap'] = 2 and ['mismatch'] = 1.
		(int 2-tuple) dims:
			Dimensions of the alignment array.
	'''
	def __str__(self):
		output = '%4s' % ''
		for i in range(len(self.cells)):
			output += '%4s' % self.sequences[0][i]
		output += '\n'
		for j in range(len(self.cells[0])):
			output += '%4s' % self.sequences[1][j]
			for i in range(len(self.cells)):	
				output += '%4d' % int(self.cells[i][j])
			output += '\n'
		return output
	
	def __init__(self, seq1, seq2):
		'''
		Constructor
		Args:
			(string) seq1, seq2:
				The two sequences that will be compared.
		'''
		self.sequences  = (' ' + seq1, ' ' + seq2)
		self.dims = (len(self.sequences[0]), len(self.sequences[1]))
		
		#Create an empty 2D array
		self.cells = [[]] * self.dims[0]
		for i in range(self.dims[0]):
			self.cells[i] = [0] * self.dims[1]
		
		#Set default penalties.
		self.penalties = {'gap':1, 'mismatch':1}
	
	def traceBack(self):
		'''
		Starts at the last cell and traces back to the root cell, adding up all of the
		values.
		Args:
			(void)
		Return:
			A string containing the alignment string.
		'''
		output = ''
		cursor = (self.dims[0] - 1, self.dims[1] - 1)
		parent = self.cells[cursor[0]][cursor[1]].parent
		while not (cursor[0] < 0 or cursor[1] < 0):
			print parent
			if(parent == (-1,0) or parent == (0,-1)):
				output = output + 'g'
			else:
				output = output + 'm'
			cursor = (cursor[0]+parent[0], cursor[1]+parent[1])
			
		return output
	
	def getValue(self, x, y):
		'''
		Gets a value from the data array.
		Args:
			(int) x:
				The index of the first sequence.
			(int) y:
				The index of the second sequence.
		Return:
			Int containing the alignment score of that cell. -1 is returned if that score
			has not been filled out yet.
		'''
		return self.cells[x][y].value
	
	def evaluateCell(self, x, y):
		'''
		Assigns a score to a cell by evaluating the sequences.
		Args:
			(int) x:
				The index of the first sequence.
			(int) y:
				The index of the second sequence.
		Return:
			Int containing the alignment score of that cell.
		'''
		#Get cost of c
		if (self.sequences[0][x] == self.sequences[1][y]):
			c = 0
		else:
			c = self.penalties['mismatch']
		
		#Check if value is filled in.
		try:
			return self.getValue(x,y)
		#Exception cause due to missing cell. Fill in cell.
		except:
			#Base cases
			if x == y == -1:
				return Cell(0)
			elif x == -1 or y == -1:
				return Cell(sys.maxint)
			#Recursive cases
			outputs = (
			           Cell(int(self.evaluateCell(x-1,y-1)) + c,(-1,-1)),
					   Cell(int(self.evaluateCell(x-1,y)) + self.penalties['gap'],(0,-1)),
					   Cell(int(self.evaluateCell(x,y-1)) + self.penalties['gap'],(-1,0))
					   )
			self.cells[x][y] = min(outputs, key=lambda cell: cell.value)
			return self.cells[x][y]
		return 0
		
	def align(self):
		'''
		Starts off the alignment process.
		Args:
			(void)
		Return:
			Int value for alignment score.
		'''
		
		self.evaluateCell(self.dims[0] - 1,self.dims[1] - 1)

if __name__ == '__main__':
	s1 = raw_input("Enter your first sequence: ")
	s2 = raw_input("Enter another sequence: ")
	
	alignment = Aligner(s1, s2)
	
	if raw_input("Change penalties? y? ") == 'y':
		goodInput = False
		while True:
			try:
				alignment.penalties['mismatch'] = int(raw_input("Mismatch penalty = "))
				break
			except:
				print "Bad input, please try again. Enter an integer please"
				continue
		while True:
			try:
				alignment.penalties['gap'] = int(raw_input("Gap penalty = "))
				goodInput = True
				break
			except:
				print "Bad input, please try again. Enter an integer please"
				continue
	
	alignment.align()
	print alignment.traceBack()
	print alignment
	