from pychembldb import *
from math import log10

#Inhibition of recombinant Syk
#Bioorg. Med. Chem. Lett. (2009) 19:1944-1949
assay = chembldb.query(Assay).filter_by(chembl_id="CHEMBL1022010").one()

for act in assay.activities:
    if act.standard_relation == "=":
        print act.compound.molecule.structure.molfile
        print "\n>  <pIC50>\n{}\n".format(9 - log10(act.standard_value))
        print "$$$$"
