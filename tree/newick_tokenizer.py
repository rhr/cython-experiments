import string, sys, re, shlex, types, itertools
from cStringIO import StringIO
import pyparsing
pyparsing.ParserElement.enablePackrat()
from pyparsing import Word, Literal, QuotedString, CaselessKeyword, \
     OneOrMore, Group, Optional, Suppress, Regex, Dict

LABELCHARS = '-.|/?#&'
META = re.compile(r'([^,=\s]+)\s*=\s*(\{[^=}]*\}|"[^"]*"|[^,]+)?')

def add_label_chars(chars):
    global LABELCHARS
    LABELCHARS += chars

class Tokenizer(shlex.shlex):
    """Provides tokens for parsing newick strings."""
    def __init__(self, infile):
        global LABELCHARS
        shlex.shlex.__init__(self, infile, posix=False)
        self.commenters = ''
        self.wordchars = self.wordchars+LABELCHARS 
        self.quotes = "'"

    def parse_embedded_comment(self):
        ws = self.whitespace
        self.whitespace = ""
        v = []
        while 1:
            token = self.get_token()
            if token == '':
                sys.stdout.write('EOF encountered mid-comment!\n')
                break
            elif token == ']':
                break
            elif token == '[':
                self.parse_embedded_comment()
            else:
                v.append(token)
        self.whitespace = ws
        return "".join(v)
        ## print "comment:", v

    def parse_ampersand_comment(s):
        word = Word(string.letters+string.digits+"%_")
        key = word.setResultsName("key") + Suppress("=")
        single_value = (Word(string.letters+string.digits+"-.") |
                        QuotedString("'") |
                        QuotedString('"'))
        range_value = Group(Suppress("{") +
                            single_value.setResultsName("min") +
                            Suppress(",") +
                            single_value.setResultsName("max") +
                            Suppress("}"))
        pair = (key + (single_value | range_value).setResultsName("value"))
        g = OneOrMore(pair)
        d = []
        for x in g.searchString(s):
            v = x.value
            if type(v) == str:
                try: v = float(v)
                except ValueError: pass
            else:
                try: v = map(float, v.asList())
                except ValueError: pass
            d.append((x.key, v))
        return d
        
