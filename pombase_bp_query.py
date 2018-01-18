import argparse
import logging

from ontobio.assoc_factory import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory

# logging.basicConfig(level="DEBUG")

POMBASE = "NCBITaxon:4896"

onto = OntologyFactory().create("go")
afactory = AssociationSetFactory()
aset = afactory.create(onto, "gene", "function", taxon=POMBASE)
bps = {}    ### Dictionary of PomBase gene ID's grouped by BP GO term

parser = argparse.ArgumentParser()
parser.add_argument("--outfile")
args = parser.parse_args()

def is_biological_process(go_term):
    ancestors = onto.ancestors(go_term)
    if "GO:0008150" in ancestors:
        return True
    else:
        return False

def is_molecular_function(go_term):
    ancestors = onto.ancestors(go_term)
    if "GO:0003674" in ancestors:
        return True
    else:
        return False

def get_ancestor_bps(mf_go_term):
    bp_ancestors = []
    for ancestor in onto.ancestors(object_id):
        if is_biological_process(ancestor):
            bp_ancestors.append(ancestor)
    return bp_ancestors

def print_results(filepath=None):
    if filepath is None:
        for key in sorted(bps, key=lambda key: len(bps[key]), reverse=False):
            print(str(len(bps[key])) + " - " + key + " - " + onto.label(key))
        print("Total BPs: " + str(len(bps)))
    else:
        with open(filepath, 'w') as f:
            f.write("Total BPs: " + str(len(bps)) + "\n")
            for key in sorted(bps, key=lambda key: len(bps[key]), reverse=True):
                f.write(str(len(bps[key])) + " - " + key + " - " + onto.label(key) + "\n")

#print(onto.ancestors("GO:0031543"))

### subject_id = PomBase:SP######.## identifier - e.g. "PomBase:SPBP19A11.06"
for subject_id in aset.association_map:
    objects_for_subject = aset.objects_for_subject(subject_id)
    for object_id in objects_for_subject:
        # Get ancestor BP terms of found MF term 
        if is_molecular_function(object_id):
            for bp in get_ancestor_bps(object_id):
                if is_biological_process(bp):
                    if bp in bps:
                        bps[bp].append(subject_id)
                    else:
                        bps[bp] = [subject_id]
            # print(str(len(get_ancestor_bps(object_id))) + " - " + object_id + " - " + onto.label(object_id))

print_results(args.outfile)
