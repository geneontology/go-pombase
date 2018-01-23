import argparse
import datetime
import logging
import math
import sys

from ontobio.assoc_factory import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory

# logging.basicConfig(level="DEBUG")

POMBASE = "NCBITaxon:4896"

start = datetime.datetime.now()
onto = OntologyFactory().create("go")
afactory = AssociationSetFactory()
aset = afactory.create(onto, "gene", "function", taxon=POMBASE)

progress_counter = 0
current_progress = None

parser = argparse.ArgumentParser()
parser.add_argument("--outfile")
parser.add_argument("--c_out")
parser.add_argument("--n")
parser.add_argument("--m")

args = parser.parse_args()



def get_ancestors(go_term):
    ### BFO:0000050 = part of
    ### I assume "subClassOf" = is a?
    all_ancestors = onto.ancestors(go_term)
    all_ancestors.append(go_term)
    subont = onto.subontology(all_ancestors)
    return subont.ancestors(go_term, relations=["subClassOf","BFO:0000050"])

def get_parents(go_term):
    all_parents = onto.parents(go_term)
    all_parents.append(go_term)
    subont = onto.subontology(all_parents)
    return subont.parents(go_term, relations=["subClassOf","BFO:0000050"])

def is_biological_process(go_term):
    ancestors = onto.ancestors(go_term)
    if "GO:0008150" in ancestors:
        return True
    else:
        return False

def get_ancestor_bps(mf_go_term):
    bp_ancestors = []
    # for ancestor in onto.ancestors(mf_go_term):
    for ancestor in get_ancestors(mf_go_term):
        if is_biological_process(ancestor):
            bp_ancestors.append(ancestor)
    return bp_ancestors

class TermAnnotationDictionary():

    def __init__(self):
        self.bps = {}    ### Dictionary of PomBase gene ID's grouped by BP GO term

    def print_results(self, filepath=None):
        if filepath is None:
            for key in sorted(self.bps, key=lambda key: len(self.bps[key]), reverse=False):
                print(str(len(self.bps[key])) + " - " + key + " - " + onto.label(key))
            print("Total BPs: " + str(len(self.bps)))
        else:
            with open(filepath, 'w') as f:
                f.write("Total BPs: " + str(len(self.bps)) + "\n")
                for key in sorted(self.bps, key=lambda key: len(self.bps[key]), reverse=True):
                    f.write(str(len(self.bps[key])) + " - " + key + " - " + onto.label(key) + "\n")

    def has_parent_with_direct_annotations_greater_than(self, go_term, n=30):
        for p in get_parents(go_term):
            if p in self.bps and len(self.bps[p]) > n:
                return True
        return False

    def select_candidate_bps(self, n=30, m=5):
        if n is None:
            n = 30
        if m is None:
            m = 5
        candidate_bps = {}
        for bp in self.bps:
            if len(self.bps[bp]) >= m and len(self.bps[bp]) <= n and self.has_parent_with_direct_annotations_greater_than(bp, n):
                candidate_bps[bp] = self.bps[bp]
        return candidate_bps

    def cluster_bps_with_similar_genes(self, m=5):
        set_list = []
        bp_dict = self.bps
        x_counter = 0
        x_progress = None
        bp_count = len(bp_dict)
        bpx_term, bpx_genes = bp_dict.popitem()
        while bpx_term is not None:
            for bpy in bp_dict:
                if bpx_term != bpy and not already_clustered(bpx_term, bpy, set_list):
                    if len(set(bpx_genes) & set(bp_dict[bpy])) >= min(m, len(bpx_genes), len(bp_dict[bpy])):
                        add_to_clusters(bpx_term, bpy, set_list)
            x_counter, x_progress = print_progress(x_counter, x_progress, bp_count)
            if len(bp_dict) == 0:
                break
            else:
                bpx_term, bpx_genes = bp_dict.popitem()
        return set_list

def already_clustered(term1, term2, bp_clusters):
    for cluster in bp_clusters:
        if term1 in cluster or term2 in cluster:
            return True
    return False

def add_to_clusters(term1, term2, bp_clusters):
    c_index = 0
    for cluster in bp_clusters:
        if term1 in cluster:
            bp_clusters[c_index].append(term2)
        elif term2 in cluster:
            bp_clusters[c_index].append(term1)
        else:
            bp_clusters.append([term1, term2])
        c_index += 1
    if len(bp_clusters) == 0:
        bp_clusters.append([term1, term2])

def pair_bp_sets_with_similar_genes(bp_dict, m=5):
    set_list = []
    x_counter = 0
    x_progress = None
    for bpx in bp_dict:
        for bpy in bp_dict:
            if bpx != bpy and (bpx, bpy) not in set_list:
                if len(set(bp_dict[bpx]) & set(bp_dict[bpy])) >= min(m, len(bp_dict[bpx]), len(bp_dict[bpy])):
                    set_list.append((bpx, bpy))
        x_counter, x_progress = print_progress(x_counter, x_progress, len(bp_dict))
    return set_list

def cluster_pair_list(result_sets):
    clusters = []
    for r_set in result_sets:
        c_index = 0
        cluster_found = False
        term1 = r_set[0]
        term2 = r_set[1]
        for cluster in clusters:
            if term1 in cluster and term2 in cluster:
                cluster_found = True
                break
            elif term1 in cluster:
                if term2 not in cluster:
                    clusters[c_index].append(term2)
                    cluster_found = True
                    break
            elif term2 in cluster:
                if term1 not in cluster:
                    clusters[c_index].append(term1)
                    cluster_found = True
                    break
            c_index += 1
        if not cluster_found:
            clusters.append([term1, term2])
    return clusters

def print_clusters(clusters, c_out=None):
    cluster_counter = 1
    if c_out is None:
        for c in clusters:
            print("---- Cluster " + str(cluster_counter) + ", Total: " + str(len(c)) + " ----")
            for t in c:
                print(t + " - " + onto.label(t))
            cluster_counter += 1
    else:
        with open(c_out, 'w') as f:
            for c in clusters:
                f.write("---- Cluster " + str(cluster_counter) + ", Total: " + str(len(c)) + " ----\n")
                for t in c:
                    f.write(t + " - " + onto.label(t) + "\n")
                cluster_counter += 1

def print_progress(progress_counter, current_progress, total):
    progress_counter += 1
    updated_progress = int(math.floor((progress_counter / total) * 100))
    if current_progress is None or updated_progress > current_progress:
        current_progress = updated_progress
        sys.stdout.flush()
        print(str(current_progress) + "% complete")
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K")
    return progress_counter, current_progress

def execution_time():
    return datetime.datetime.now() - start

tad = TermAnnotationDictionary()

### subject_id = PomBase:SP######.## identifier - e.g. "PomBase:SPBP19A11.06"
for subject_id in aset.association_map:
    objects_for_subject = aset.objects_for_subject(subject_id)
    for object_id in objects_for_subject:
        if is_biological_process(object_id):
            ancestor_bps = get_ancestor_bps(object_id)
            ancestor_bps.append(object_id)
            for bp in ancestor_bps:
                if bp in tad.bps:
                    if subject_id not in tad.bps[bp]:
                        tad.bps[bp].append(subject_id)
                else:
                    tad.bps[bp] = [subject_id]


# print_results(bps, args.outfile)
if args.n is not None:
    args.n = int(args.n)
if args.m is not None:
    args.m = int(args.m)

print("Initial BP count: " + str(len(tad.bps)))
new_bps = tad.select_candidate_bps(args.n, args.m)
print("Candidate BP count: " + str(len(new_bps)))
p_list = pair_bp_sets_with_similar_genes(new_bps, args.m)
print("BP pair count: " + str(len(p_list)))
c_list = cluster_pair_list(p_list)
print_clusters(c_list, args.c_out)