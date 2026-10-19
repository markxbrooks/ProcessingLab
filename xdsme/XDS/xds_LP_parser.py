#!/usr/bin/env python

def xdsLpParser(lptxt):
    from string import strip, split
    #eq2tup = lambda s: tuple(split(s,'='))
    splitnames = lambda s: split(s,'*\n')
    parts = lambda s, delim: map(strip, split(s, delim))
    useful = lambda ss: filter(lambda s: s and s[0]!=';', ss)
    #config = map(lambda _:(_[0], map(eq2tup, useful(parts(_[1],'\n')))),
    #               map(splitnames, useful(parts(txt,'['))) )
    config = map(splitnames, lptxt)
    pprint.pprint(config)

#xdsLpParser(txt)
import re
def re_show(pat, s):
    print re.compile(pat, re.M).sub("{\g<0>}", s.rstrip()),'\n'

txt = open("IDXREF.LP").read()
re_show(r'^ \*{5}.*\*{5}$',txt)
sections = re.compile(r'^ \*{5,} (.*) \*{5,}$',re.M)
sections.sub("{\g<0>}", txt.rstrip())
sections_match = sections.search(txt)
sections_match
