import datetime
from ontobio.golr.golr_query import GolrAssociationQuery, GolrFields
from ontobio.assoc_factory import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory
# from pombase_direct_bp_annots_query import TermAnnotationDictionary, is_molecular_function, is_cellular_component, setup_pombase, pair_bp_sets_with_similar_genes, uniqueify
from pombase_direct_bp_annots_query import TermAnnotationDictionary, GOTermAnalyzer

M = GolrFields
POMBASE = "NCBITaxon:4896"

my_select_fields = [
                M.ID,
                M.IS_DEFINED_BY,
                M.SOURCE,
                M.SUBJECT,
                M.SUBJECT_LABEL,
                M.SUBJECT_TAXON,
                M.SUBJECT_TAXON_LABEL,
                M.RELATION,
                M.RELATION_LABEL,
                M.OBJECT,
                M.OBJECT_LABEL,
                M.OBJECT_TAXON,
                M.OBJECT_TAXON_LABEL,
                M.EVIDENCE_GRAPH,

                M.EVIDENCE_WITH,
                M.EVIDENCE_TYPE,
                M.REFERENCE,
                M.ANNOTATION_EXTENSION_JSON # special!
            ]
mf_part_of_relations = ['part_of']
mf_gene_relations =['has input','has_direct_input','has_regulation_target','regulates activity of']
cc_relations = ['occurs in','occurs at']

class GeneConnectionSet():
    def __init__(self):
        self.gene_connections = []

    def append(self, gene_connection):
        self.gene_connections.append(gene_connection)

    def contains(self, gene_connection):
        for gc in self.gene_connections:
            if gc.equals(gene_connection):
                return True
        return False

    def merge(self, other_set):
        for connection in other_set.gene_connections:
            if not self.contains(connection):
                self.append(connection)

class GeneConnection():
    def __init__(self, gene_a, gene_b, direction, annotation=None):
        self.gp_a = gene_a
        self.gp_b = gene_b
        self.direction = direction
        self.annotation = annotation

    def equals(self, gene_connection):
        if self.gp_a == gene_connection.gp_a and self.gp_b == gene_connection.gp_b and self.direction == gene_connection.direction:
            return True
        else:
            return False

def get_specific_annots(annots, subject_id, object_id):
    found_annots = []
    for annot in annots:
        if annot["subject"]["id"] == subject_id and annot["object"]["id"] == object_id:
            found_annots.append(annot)
    return found_annots

def query_for_annots(subject_id=None, object_id=None):
    q = GolrAssociationQuery('gene', 'function', subject_taxon=POMBASE, select_fields=my_select_fields,
                            subject=subject_id, 
                            object=object_id,
                            rows=100000) # Only 39435 came back when I set this higher
    a = q.exec()["associations"]
    return a

class AnnotationDataExtracter():
    def __init__(self, term_analyzer):
        self.analyzer = term_analyzer # GOTermAnalyzer()

    def get_annots_by_relationship(self, annotations, relations, ids=None):
        annots = []
        for annot in annotations:
            relevant_ext = None
            if 'annotation_extensions' in annot:
                relevant_ext = self.relevant_extension(annot["annotation_extensions"], relations, ids)
            if relevant_ext is not None:
                annots.append(annot)
        return annots

    def relevant_extension(self, exts, relations, ids=None):
        for ext in exts:
            for rel in ext["relationship"]["relation"]:
                if ids is None:
                    if rel["label"] in relations:
                        return ext
                elif rel["label"] in relations and ext["relationship"]["id"] in ids:
                    return ext

    def relevant_mf_extension(self, annotation, bp, tad):
        if 'annotation_extensions' in annotation:
            # 1a - extension having "part_of(bp)"
            ext = self.relevant_extension(annotation["annotation_extensions"], mf_part_of_relations, [bp])
            if ext is not None:
                return ext

            # 1b - extension having "has input(other gene also under same bp)" or "has_regulation_target(same bp)"
            ext = self.relevant_extension(annotation["annotation_extensions"], mf_gene_relations, tad.bps[bp])
            if ext is not None:
                return ext

    def relevant_mf_annotation(self, annotations, bp, tad):
        # 1a - extension having "part_of(bp)"
        mf_annots = self.get_annots_by_relationship(annotations, mf_part_of_relations, [bp])
        if len(mf_annots) > 0:
            return mf_annots[0]

        # 1b - extension having "has input(other gene also under same bp)" or "has_regulation_target(same bp)"
        mf_annots = self.get_annots_by_relationship(annotations, mf_gene_relations, tad.bps[bp])
        if len(mf_annots) > 0:
            return mf_annots[0]
        
        # 1c - is there another MF annotation
        if len(annotations) > 0:
            return annotations[0]

    def relevant_cc_annotation(self, annotations):
        occurs_annots = self.get_annots_by_relationship(annotations, cc_relations)
        # Go back in, find 'occurs in' relationship, and check if term is CC
        for annot in occurs_annots:
            if self.access_extensions(annot, self.is_extension_object_a_cc):
                return annot

    def cc_ext_description(self, ext):
        descriptions = []
        for rel in ext["relationship"]["relation"]:
            if rel["label"] in cc_relations:
                relation = rel["label"]
                object_id = ext["relationship"]["id"]
                descriptions.append(relation + "(" + object_id + ") - " + self.analyzer.label(object_id))
        if len(descriptions) > 0:
            return descriptions

    def is_extension_object_a_cc(self, ext):
        if self.analyzer.is_cellular_component(ext["relationship"]["id"]):
            return True
        else:
            return False

    def access_extensions(self, annotation, check_func):
        return_values = []
        if "annotation_extensions" in annotation:
            for ext in annotation["annotation_extensions"]:
                result = check_func(ext)
                if type(result) is bool:
                    return result
                elif result is not None:
                    return_values.append(result)
        return return_values

    def direct_cc_annotations(self, annotations):
        cc_annots = []
        for annot in annotations:
            if self.analyzer.is_cellular_component(annot["object"]["id"]):
                cc_annots.append(annot)
        return cc_annots

    def get_with_gene_connections(self, annotations, bp, tad):
        connections = GeneConnectionSet()
        for annot in annotations:
            # annotation of the gene product to the GO term “protein binding” 
            if "evidence_with" in annot and annot["object"]["id"] == "GO:0005515":
                with_genes = annot["evidence_with"]
                for wg in with_genes:
                    connection = GeneConnection(annot["subject"]["id"], wg, "<->", annot)
                    if wg in tad.bps[bp] and not connections.contains(connection):
                        connections.append(connection)
        return connections

    # We really need a method that returns the annotation along with only the relevant ext/relation
    def get_gene_connections(self, annotations, bp, tad):
        connections = GeneConnectionSet()
        # 1b - extension having "has input(other gene also under same bp)" or "has_regulation_target(same bp)"
        annots = self.get_annots_by_relationship(annotations, mf_gene_relations, tad.bps[bp])
        for annot in annots:
            connected_genes = self.access_extensions(annot, self.get_gene_connection)
            for gene_list in connected_genes:
                for gene in gene_list:
                    connection = None
                    if gene[1] in ["has_regulation_target","regulates activity of"]:
                        connection = GeneConnection(annot["subject"]["id"], gene[0], "-->", annot)
                    elif gene[1] in ["has input","has_direct_input"]:
                        connection = GeneConnection(annot["subject"]["id"], gene[0], "<--", annot)
                    # connection = (annot["subject"]["id"],gene[0])
                    if connection and gene[0] in tad.bps[bp] and not connections.contains(connection):
                        connections.append(connection)
        connections.merge(self.get_with_gene_connections(annotations, bp, tad))
        return connections

    def get_gene_connection(self, ext):
        connections = []
        for rel in ext["relationship"]["relation"]:
            if rel["label"] in mf_gene_relations:
                connections.append([ext["relationship"]["id"], rel["label"]])
        if len(connections) > 0:
            return connections

    def print_results(self, gene_info, tad, outfile=None):
        is_first_in_list = True
        if outfile is None:
            for k in gene_info:
                if is_first_in_list:
                    print(gene_info[k]["bp"])
                    print("---------------------------")
                    is_first_in_list = False
                if "molecular_function" in gene_info[k]:
                    mf = gene_info[k]["molecular_function"]
                    print(k + " - " + mf["object"]["id"] + " - " + mf["object"]["label"])
                    mf_ext = self.relevant_mf_extension(mf, gene_info[k]["bp"], tad)
                    if mf_ext is not None:
                        relation = mf_ext["relationship"]["relation"][0]["label"]
                        object_id = mf_ext["relationship"]["id"]
                        print("  Ext: " + relation + "(" + object_id + ") - " + self.analyzer.label(object_id))
                if "cellular_component" in gene_info[k]:
                    cc = gene_info[k]["cellular_component"]
                    print(k + " - " + cc["object"]["id"] + " - " + cc["object"]["label"])
                    for d_list in self.access_extensions(cc, self.cc_ext_description):
                        print("  Ext: " + "; ".join(d_list))
        else:
            with open(outfile, 'a') as f:
                for k in gene_info:
                    if is_first_in_list:
                        f.write("BP: " + gene_info[k]["bp"] + " - " + self.analyzer.label(gene_info[k]["bp"]) + "\n")
                        f.write("---------------------------\n")
                        is_first_in_list = False
                    if "molecular_function" in gene_info[k]:
                        mf = gene_info[k]["molecular_function"]
                        f.write(k + " - " + mf["object"]["id"] + " - " + mf["object"]["label"] + "\n")
                        mf_ext = self.relevant_mf_extension(mf, gene_info[k]["bp"], tad)
                        if mf_ext is not None:
                            relation = mf_ext["relationship"]["relation"][0]["label"]
                            object_id = mf_ext["relationship"]["id"]
                            f.write("  Ext: " + relation + "(" + object_id + ") - " + self.analyzer.label(object_id) + "\n")
                    if "cellular_component" in gene_info[k]:
                        cc = gene_info[k]["cellular_component"]
                        f.write(k + " - " + cc["object"]["id"] + " - " + cc["object"]["label"] + "\n")
                        for d_list in self.access_extensions(cc, self.cc_ext_description):
                            f.write("  Ext: " + "; ".join(d_list) + "\n")
                f.write("\n")

def do_stuff():
    start = datetime.datetime.now()

    pombase_annots = query_for_annots()

    ontology = OntologyFactory().create("go")
    afactory = AssociationSetFactory()
    aset = afactory.create(ontology, "gene", "function", taxon=POMBASE)
    tad = TermAnnotationDictionary(ontology, aset)
    analyzer = GOTermAnalyzer(ontology)
    extracter = AnnotationDataExtracter(analyzer)
    # tad.get_our_nice_lists(n=30, m=10)
    # result_bp_list = uniqueify(tad.pair_list)

    # Currently used relations in extensions field:
    #   ['has_direct_input', 'part_of', 'exists_during', 'occurs in', 'has_regulation_target', 
    #    'coincident with', 'happens_during', 'activated_by', 'has input', 'directly negatively regulates', 
    #    'occurs at', 'inhibited_by', 'directly positively regulates', 'not_exists_during', 
    #    'regulates activity of', 'not_happens_during']

    bps = ["GO:0010971"]
    # bps = ["GO:0010971", "GO:1902751"]
    # bps = ["GO:0000002"] # mitochondrial genome maintenance
    # bps = ["GO:0006903"] # vesicle targeting
    # bps = result_bp_list
    for bp in bps:
        ok_to_print_results = False
        gene_info = {}
        for g in tad.bps[bp]:
            ### Find annots for g subject where go_term is MF - then check for extension part_of "GO:0010971"
            gene_annots = []
            for annot in pombase_annots:
                if annot["subject"]["id"] == g:
                    gene_annots.append(annot)
            mf_annots = []
            for annot in gene_annots:
                if analyzer.is_molecular_function(annot["object"]["id"]):
                    mf_annots.append(annot)
            gene_info[g] = {}
            gene_info[g]["bp"] = bp
            mf = extracter.relevant_mf_annotation(mf_annots, bp, tad) # 2a,b,c
            if mf is not None:
                gene_info[g]["molecular_function"] = mf
            cc = extracter.relevant_cc_annotation(mf_annots) # 2a
            if cc is not None:
                # ok_to_print_results = True
                gene_info[g]["cellular_component"] = cc
            cc_annots = []
            if "cellular_component" not in gene_info[g]:       
                cc_annots = extracter.direct_cc_annotations(gene_annots)
                if len(cc_annots) > 0:
                    gene_info[g]["cellular_component"] = cc_annots[0] # 2b
            gene_info[g]["connections"] = extracter.get_gene_connections(mf_annots, bp, tad)

        # print_results(gene_info, tad, "GO:0010971_model_attempt.txt")
        # conn_outfile = "GO:0010971_connections.txt"
        # if ok_to_print_results:
        #     print_results(gene_info, tad, "models_w_CC_extension_annotations.txt")
        # print_results(gene_info, tad)
        conn_outfile = None
        # with open(conn_outfile, 'w') as f:
        for g in gene_info:
            for connect in gene_info[g]["connections"].gene_connections:
                object_id = connect.annotation["object"]["id"]
                print(connect.gp_a + "(" + aset.label(connect.gp_a) + ") " + connect.direction + " " + connect.gp_b + "(" + aset.label(connect.gp_b) + ") - " + object_id + " (" + aset.label(object_id) + ")")
                # f.write(connect[0] + "(" + aset.label(connect[0]) + ") -- " + connect[1] + "(" + aset.label(connect[1]) + ")\n")
    print("Execution time: " + str(datetime.datetime.now() - start))

# do_stuff()