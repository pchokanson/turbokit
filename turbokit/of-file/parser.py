

import io
import string
import tokenizer
from tokenizer import *

def parse_file(f, type_="dictionary"):
	gen_fn = tokenizer.tokenize(f)
	if type_ == "dictionary":
		return parse_dict(gen_fn)

def next_datum(t_iter):
	token = next(t_iter)
	if token._type == TokenTypes.START_LINE_COMMENT:
		parse_line_comment(t_iter)
	elif token._type == TokenTypes.START_C_COMMENT:
		parse_block_comment(t_iter)
	return token

def parse_line_comment(t_iter):
	for token in t_iter:
		if token._type == TokenTypes.EOF:
			break

def parse_block_comment(t_iter):
	for token in t_iter:
		if token._type == TokenTypes.END_C_COMMENT:
			break

def parse_dict(t_iter):
	key = next_datum(t_iter)
	print("KEY: " + repr(key))

def parse_list(t_iter):
	pass

if __name__ == "__main__":
	f = open("../../temp/meridionalPumpProfile/system/controlDict")
	parse_file(f)
	f.close()
