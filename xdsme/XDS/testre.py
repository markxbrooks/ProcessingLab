import re

fl = open("INTEGRATE.LP").readlines()

_rfl, _rI = r"[\ ]+([0-9\.]+)", r"[\ ]+([\d]+) "
scalegr =   r" "+_rI+r"  (\d)"+_rfl+r"  ....... "+_rI
scalegr += r" .....    .... "+_rI+2*(_rfl)
sc = re.compile(scalegr)



for l in fl:
     hit = sc.search(l)
     if hit:
         print l[:-1], hit.groups()#, map(float, hit.groups())


