import haiku

if __name__ == '__main__':
	generator = haiku.HaikuGenerator('dictionary.dict')
	outFile = open('out.txt','w+')
	outFile.write(haiku.generate())
	outFile.close()