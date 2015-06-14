# tokenizer.py
# 

import io
import string

class TokenTypes(object):
	UNDEFINED = -2
	EOF = -1
	EOL = 0
	START_C_COMMENT = 1
	END_C_COMMENT = 2
	START_LINE_COMMENT = 3
	START_DICT = 4
	END_DICT = 5
	START_LIST = 6
	END_LIST = 7
	

class Token(object):
	def __init__(self, string_):
		self._str = string_
		if string_ == -1:
			self._type = TokenTypes.EOF
		elif string_ == '\n':
			self._type = TokenTypes.EOL
		elif string_[:1] == r"/*":
			self._type = TokenTypes.START_C_COMMENT
		elif string_[-1:] == r"*/":
			self._type = TokenTypes.END_C_COMMENT
		elif string_[:1] == r"//":
			self._type = TokenTypes.START_LINE_COMMENT
		elif string_ == "{":
			self._type = TokenTypes.START_DICT
		elif string_ == "}":
			self._type = TokenTypes.END_DICT
		elif string_ == "(":
			self._type = TokenTypes.START_LIST
		elif string_ == ")":
			self._type = TokenTypes.END_LIST
		else:
			print("UNDEFINED: " + repr(string_))
			self._type = TokenTypes.UNDEFINED
		
		
	def __str__(self):
		return self._str
	
	def token_id(self):
		pass

def tokenize(ibuffer):
	for line in ibuffer:
		for t in line.split():
			yield Token(t)
		yield Token("\n")
	yield Token(-1)

# Tokens
if __name__ == "__main__":
	f = open("../../temp/meridionalPumpProfile/constant/polyMesh/blockMeshDict")
	for t in tokenize(f):
		print(t)
	f.close()
