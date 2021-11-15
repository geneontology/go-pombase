import datetime
from ontobio.golr.golr_query import GolrAssociationQuery, GolrFields
from ontobio.assoc_factory import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory
from ontobio.rdfgen import relations
from ontobio.model.association import GoAssociation, Curie
from prefixcommons import curie_util
# from pombase_direct_bp_annots_query import TermAnnotationDictionary, is_molecular_function, is_cellular_component, setup_pombase, pair_bp_sets_with_similar_genes, uniqueify
from pombase_direct_bp_annots_query import TermAnnotationDictionary, GOTermAnalyzer, ProgressTracker
from gaf_annotation_set import GafAnnotationSet
import json
import typing


curie_util.default_curie_maps.append({"GOREL": "http://purl.obolibrary.org/obo/GOREL_"})
WITH_SUPPORT_FROM = Curie("RO", "0002614")  # placeholder; do not emit - actually 'is evidence with support from'
# WITH_SUPPORT_FROM = "RO:0002614"


def convert_relation_labels_to_curies(relation_labels: typing.List):
    relation_curies = []
    for lbl in relation_labels:
        relation_uri = relations.lookup_label(lbl)
        relation_curie = Curie.from_str(curie_util.contract_uri(relation_uri)[0])
        # if relation_curie:
        #     relation_curies.append(relation_curie)
        relation_curies.append(relation_curie)
    return relation_curies


mf_part_of_relations = ['part_of']
regulation_relations = ['has_regulation_target',
                    'regulates activity of', 'directly_positively_regulates',
                    'directly_negatively_regulates']
input_relations = ['has input','has_direct_input']
# mf_gene_relations = regulation_relations + input_relations
cc_relations = ['occurs in','occurs at']
# convert all these to CURIE cuz that's how ontobio talks now
mf_part_of_relations_curie = convert_relation_labels_to_curies(mf_part_of_relations)
regulation_relations_curie = convert_relation_labels_to_curies(regulation_relations)
input_relations_curie = convert_relation_labels_to_curies(input_relations)
cc_relations_curie = convert_relation_labels_to_curies(cc_relations)

# ontology = OntologyFactory().create("go")
afactory = AssociationSetFactory()

class ExtensionGolrFields(GolrFields):
    ANNOTATION_EXTENSION_JSON="annotation_extension_json"
    EVIDENCE_WITH="evidence_with"
    REFERENCE="reference"
    EVIDENCE_TYPE="evidence_type"

class ExtensionGolrAssociationQuery(GolrAssociationQuery):
    def __init__(self, subject_category=None, object_category=None, **kwargs):
        GolrAssociationQuery.__init__(self, subject_category, object_category, **kwargs)
        self.select_fields = self.get_select_fields()

    def translate_doc(self, d, field_mapping=None, map_identifiers=None, **kwargs):
        # assoc = super(GolrAssociationQuery, self).translate_doc(d, field_mapping, map_identifiers, kwargs)
        assoc = GolrAssociationQuery.translate_doc(self, d, field_mapping, map_identifiers, **kwargs)

        if M.ANNOTATION_EXTENSION_JSON in d:
            assoc['annotation_extensions'] = [json.loads(ext) for ext in d[M.ANNOTATION_EXTENSION_JSON]]
        if M.EVIDENCE_TYPE in d:
            assoc[M.EVIDENCE_TYPE] = d[M.EVIDENCE_TYPE]
        if M.EVIDENCE_WITH in d:
            assoc[M.EVIDENCE_WITH] = d[M.EVIDENCE_WITH]
        if M.REFERENCE in d:
            assoc[M.REFERENCE] = d[M.REFERENCE]

        return assoc

    def get_select_fields(self):
        M = ExtensionGolrFields()

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

        return my_select_fields

M = ExtensionGolrFields()
POMBASE = "NCBITaxon:4896"

class GeneConnectionSet():
    def __init__(self):
        self.gene_connections = []

    def append(self, gene_connection):
        self.gene_connections.append(gene_connection)

    def contains(self, gene_connection):
        for gc in self.gene_connections:
            if gc.equals(gene_connection):
                return True
            if gene_connection.relation == WITH_SUPPORT_FROM and (self.find(gene_connection.gp_a, gene_connection.gp_b) or self.find(gene_connection.gp_b, gene_connection.gp_a)):
                return True
        return False

    def merge(self, other_set):
        for connection in other_set.gene_connections:
            if connection.annotation.object.id == Curie("GO", "0005515"):  # Check if triple through extension already exists
                if self.find(connection.gp_a, connection.gp_b, Curie("RO", "0002400")):  # has_direct_input
                    continue
            if not self.contains(connection):
                self.append(connection)

    def find(self, gp_a, gp_b, relation=None):
        connections = []
        for connection in self.gene_connections:
            if connection.gp_a == gp_a and connection.gp_b == gp_b:
                if relation is None or connection.relation == relation:
                    connections.append(connection)
        return connections

class GeneConnection():
    def __init__(self, gene_a: Curie, gene_b: Curie, object_id: Curie, relation: Curie, annotation: GoAssociation=None):
        self.gp_a = gene_a
        self.gp_b = gene_b
        self.object_id = object_id
        self.relation = relation
        self.annotation = annotation

    def equals(self, gene_connection):
        if self.gp_a == gene_connection.gp_a and self.gp_b == gene_connection.gp_b and self.relation == gene_connection.relation:
            return True
        else:
            return False

    def print_connection(self, a_set):
        # Need a better way to pass in GP info to get labels
        # a_set = afactory.create(ontology, subject_category='gene', object_category='function')  # Labels for debugging
        return str(self.gp_a) + " (" + a_set.label(self.gp_a) + ") " + str(self.relation) + " " + str(self.gp_b) + " (" + a_set.label(str(self.gp_b)) + ") through " + str(self.object_id) + " (" + a_set.label(str(self.object_id)) + ")"

def get_specific_annots(annots, subject_id, object_id):
    found_annots = []
    for annot in annots:
        if annot["subject"]["id"] == subject_id and annot["object"]["id"] == object_id:
            found_annots.append(annot)
    return found_annots

def query_for_annots(subject_id=None, object_id=None):
    q = ExtensionGolrAssociationQuery('gene', 'function', subject_taxon=POMBASE,
                                subject=subject_id,
                                object=object_id,
                                rows=10)
    a = q.exec()["associations"]
    return a


class AnnotationDataExtracter:
    def __init__(self, term_annotation_dictionary: TermAnnotationDictionary, gaf_annotation_set: GafAnnotationSet, gene_relations=None):
        # self.analyzer = term_analyzer # GOTermAnalyzer()
        self.tad = term_annotation_dictionary
        self.gas = gaf_annotation_set
        if gene_relations is None:
            self.mf_gene_relations = regulation_relations_curie + input_relations_curie
        else:
            self.mf_gene_relations = gene_relations

    @classmethod
    def parse_gaf(AnnotationDataExtracter, gaf_filename, go_ontology, json_file=None):
        gas = GafAnnotationSet(gaf_filename, go_ontology, filter_evidence=True)
        gas.filter_evidence()

        tad = TermAnnotationDictionary(go_ontology, gas.association_set, json_file)
        return AnnotationDataExtracter(tad, gas)

    def genes_and_annots_for_bp(self, bp_term):
        return self.calculate_genes_and_annots_for_bp(bp_term)

    def calculate_genes_and_annots_for_bp(self, bp_term):
        if bp_term not in self.tad.bps:
            return None
        gene_info = {}
        progress = ProgressTracker(len(self.tad.bps[bp_term]), "get relevant annotations/connections for each gene")
        for g in self.tad.bps[bp_term]:
            ### Find annots for g subject where go_term is MF - then check for extension part_of "GO:0010971"
            gene_annots = self.gas.annotations_for_subject(g)
            mf_annots = []
            for annot in gene_annots:
                if self.tad.analyzer.is_molecular_function(str(annot.object.id)):
                    mf_annots.append(annot)
            gene_info[g] = {}
            gene_info[g]["bp"] = bp_term
            mf = self.relevant_mf_annotation(mf_annots, bp_term)  # 2a,b,c
            if mf is not None:
                gene_info[g]["molecular_function"] = mf
            cc = self.relevant_cc_annotation(mf_annots)  # 2a
            if cc is not None:
                # ok_to_print_results = True
                gene_info[g]["cellular_component"] = [cc]
            if "cellular_component" not in gene_info[g]:
                all_cc_annotations = self.cc_annotations(gene_annots)
                if all_cc_annotations and len(all_cc_annotations) == 1:
                    gene_info[g]["cellular_component"] = all_cc_annotations
            gene_info[g]["connections"] = self.get_gene_connections(mf_annots, bp_term)

            progress.print_progress()

        return gene_info

    def get_annots_by_relationship(self, annotations, relations, ids=None):
        annots = []
        for annot in annotations:
            relevant_ext = None
            if annot.object_extensions:
                relevant_ext = self.relevant_object_extension(annot.object_extensions, relations, ids)
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

    def relevant_object_extension(self, exts, relations, ids=None):
        # relations_curie = convert_relation_labels_to_curies(relations)
        for ext in exts:
            for rel in ext.elements:
                if ids is None:
                    if rel.relation in relations:
                        return ext
                elif rel.relation in relations and str(rel.term) in ids:
                    return ext

    def relevant_mf_extension(self, annotation, bp, tad):
        if annotation.object_extensions:
            extensions = annotation.object_extensions
            # 1a - extension having "part_of(bp)"
            ext = self.relevant_object_extension(extensions, mf_part_of_relations_curie, [bp])
            if ext is not None:
                return ext

            # 1b - extension having "has input(other gene also under same bp)" or "has_regulation_target(same bp)"
            ext = self.relevant_object_extension(extensions, self.mf_gene_relations, tad.bps[bp])
            if ext is not None:
                return ext

    def relevant_mf_annotation(self, annotations, bp):
        # 1a - extension having "part_of(bp)"
        mf_annots = self.get_annots_by_relationship(annotations, mf_part_of_relations_curie, [bp])
        if len(mf_annots) > 0:
            return mf_annots[0]

        # 1b - extension having "has input(other gene also under same bp)" or "has_regulation_target(same bp)"
        mf_annots = self.get_annots_by_relationship(annotations, self.mf_gene_relations, self.tad.bps[bp])
        if len(mf_annots) > 0:
            return mf_annots[0]
        
        # 1c - is there another MF annotation
        if len(annotations) > 0:
            return annotations[0]

    def relevant_cc_annotation(self, annotations):
        occurs_annots = self.get_annots_by_relationship(annotations, cc_relations_curie)
        # Go back in, find 'occurs in' relationship, and check if term is CC
        for annot in occurs_annots:
            if self.access_extensions(annot, self.is_extension_object_a_cc):
                return annot

    def cc_annotations(self, annotations):
        cc_annots = self.direct_cc_annotations(annotations)
        if len(cc_annots) > 0:
            return cc_annots # 2b

    def cc_ext_description(self, ext):
        descriptions = []
        if "relationship" in ext:
            for rel in ext["relationship"]["relation"]:
                if rel["label"] in cc_relations_curie:
                    relation = rel["label"]
                    object_id = ext["relationship"]["id"]
                    descriptions.append(relation + "(" + object_id + ") - " + self.tad.analyzer.label(object_id))
        elif "intersection_of" in ext:
            for rel in ext["intersection_of"]:
                if rel["property"] in cc_relations_curie:
                    relation = rel["property"]
                    object_id = ext["filler"]
                    descriptions.append(relation + "(" + object_id + ") - " + self.tad.analyzer.label(object_id))
        if len(descriptions) > 0:
            return descriptions

    def is_extension_object_a_cc(self, extensions):
        for ext in extensions.elements:
            if self.tad.analyzer.is_cellular_component(str(ext.term)):
                return True
        return False

    def access_extensions(self, annotation, check_func):
        return_values = []
        extensions = annotation.object_extensions
        for ext in extensions:
            result = check_func(ext)
            if type(result) is bool:
                return result
            elif result is not None:
                return_values.append(result)
        return return_values

    def direct_cc_annotations(self, annotations):
        cc_annots = []
        for annot in annotations:
            if self.tad.analyzer.is_cellular_component(str(annot.object.id)):
                cc_annots.append(annot)
        return cc_annots

    def get_with_gene_connections(self, annotations, bp):
        connections = GeneConnectionSet()
        for annot in annotations:
            # annotation of the gene product to the GO term “protein binding”
            if annot.evidence.with_support_from and str(annot.object.id) == "GO:0005515":
                with_genes = annot.evidence.with_support_from[0].elements  # TODO: Assuming there were no "|" separators in with/from
                for wg in with_genes:
                    connection = None
                    if str(annot.subject.id) != wg and str(wg) in self.tad.bps[bp]:
                        connection = GeneConnection(annot.subject.id,
                                                    wg,
                                                    annot.object.id,
                                                    WITH_SUPPORT_FROM,
                                                    annot)
                    if connection is not None and not connections.contains(connection):
                        connections.append(connection)
        return connections

    # We really need a method that returns the annotation along with only the relevant ext/relation
    def get_gene_connections(self, annotations, bp):
        connections = GeneConnectionSet()
        # 1b - extension having "has input(other gene also under same bp)" or "has_regulation_target(same bp)"
        annots = self.get_annots_by_relationship(annotations, self.mf_gene_relations, self.tad.bps[bp])
        for annot in annots:
            connected_genes = self.access_extensions(annot, self.get_gene_connection)
            for gene_list in connected_genes:
                for gene in gene_list:
                    connection = None
                    if gene[1] in regulation_relations_curie\
                            or (gene[1] in input_relations_curie and annot.subject.id != gene[0]):
                        connection = GeneConnection(annot.subject.id,
                                                    gene[0],
                                                    annot.object.id,
                                                    gene[1],
                                                    annot)
                    # connection = (annot["subject"]["id"],gene[0])
                    if connection and str(gene[0]) in self.tad.bps[bp] and not connections.contains(connection):
                        connections.append(connection)
        connections.merge(self.get_with_gene_connections(annotations, bp))
        return connections

    def get_gene_connection(self, extensions):
        connections = []
        for ext in extensions.elements:
            if ext.relation in self.mf_gene_relations:
                # return [gene_b, relation_curie]
                connections.append([ext.term, ext.relation])
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
                    # print(k + " - " + mf["object"]["id"] + " - " + mf["object"]["label"])
                    print(k + " - " + mf["object"]["id"] + " - " + self.object_label(mf["object"], tad.ontology))
                    mf_ext = self.relevant_mf_extension(mf, gene_info[k]["bp"], tad)
                    if mf_ext is not None:
                        # relation = mf_ext["relationship"]["relation"][0]["label"]
                        relation = self.relation_label(mf_ext)
                        # object_id = mf_ext["relationship"]["id"]
                        object_id = self.relation_object_id(mf_ext)
                        object_label = self.tad.analyzer.label(object_id)
                        object_label = object_label if object_label is not None else object_id
                        print("  Ext: " + relation + "(" + object_id + ") - " + object_label)
                if "cellular_component" in gene_info[k]:
                    cc = gene_info[k]["cellular_component"]
                    print(k + " - " + cc["object"]["id"] + " - " + self.object_label(cc["object"], tad.ontology))
                    for d_list in self.access_extensions(cc, self.cc_ext_description):
                        print("  Ext: " + "; ".join(d_list))
        else:
            with open(outfile, 'a') as f:
                for k in gene_info:
                    if is_first_in_list:
                        f.write("BP: " + gene_info[k]["bp"] + " - " + self.tad.analyzer.label(gene_info[k]["bp"]) + "\n")
                        f.write("---------------------------\n")
                        is_first_in_list = False
                    if "molecular_function" in gene_info[k]:
                        mf = gene_info[k]["molecular_function"]
                        f.write(k + " - " + mf["object"]["id"] + " - " + mf["object"]["label"] + "\n")
                        mf_ext = self.relevant_mf_extension(mf, gene_info[k]["bp"], tad)
                        if mf_ext is not None:
                            relation = mf_ext["relationship"]["relation"][0]["label"]
                            object_id = mf_ext["relationship"]["id"]
                            f.write("  Ext: " + relation + "(" + object_id + ") - " + self.tad.analyzer.label(object_id) + "\n")
                    if "cellular_component" in gene_info[k]:
                        cc = gene_info[k]["cellular_component"]
                        f.write(k + " - " + cc["object"]["id"] + " - " + cc["object"]["label"] + "\n")
                        for d_list in self.access_extensions(cc, self.cc_ext_description):
                            f.write("  Ext: " + "; ".join(d_list) + "\n")
                f.write("\n")

    def object_label(self, obj_dict, onto):
        if "label" in obj_dict:
            return obj_dict["label"]
        else:
            return onto.label(obj_dict["id"])

    ### Needed because golr extension json and gaf parsed structure differ
    def relation_label(self, obj_ext):
        if "relationship" in obj_ext:
            return obj_ext["relationship"]["relation"][0]["label"]
        else:
            return obj_ext["intersection_of"][0]["property"]

    def relation_object_id(self, obj_ext):
        if "relationship" in obj_ext:
            return obj_ext["relationship"]["id"]
        else:
            return obj_ext["intersection_of"][0]["filler"]

def genes_and_annots_for_bp(bp_term):
    # pombase_annots = query_for_annots()
    pombase_annots = []

    ontology = OntologyFactory().create("go")
    afactory = AssociationSetFactory()
    aset = afactory.create(ontology, "gene", "function", taxon=POMBASE)
    tad = TermAnnotationDictionary(ontology, aset)
    analyzer = GOTermAnalyzer(ontology)
    extracter = AnnotationDataExtracter(tad, aset)

    ok_to_print_results = False
    gene_info = {}
    progress = ProgressTracker(len(tad.bps[bp_term]), "get relevant annotations/connections for each gene")
    for g in tad.bps[bp_term]:
        ### Find annots for g subject where go_term is MF - then check for extension part_of "GO:0010971"
        gene_annots = query_for_annots(g)
        # for annot in pombase_annots:
        #     if annot["subject"]["id"] == g:
        #         gene_annots.append(annot)
        mf_annots = []
        for annot in gene_annots:
            if analyzer.is_molecular_function(annot["object"]["id"]):
                mf_annots.append(annot)
        gene_info[g] = {}
        gene_info[g]["bp"] = bp_term
        mf = extracter.relevant_mf_annotation(mf_annots, bp_term) # 2a,b,c
        if mf is not None:
            gene_info[g]["molecular_function"] = mf
        cc = extracter.relevant_cc_annotation(mf_annots) # 2a
        if cc is not None:
            gene_info[g]["cellular_component"] = [cc]
        if "cellular_component" not in gene_info[g]:
        #     cc_annots = extracter.direct_cc_annotations(gene_annots)
        #     if len(cc_annots) > 0:
        #         gene_info[g]["cellular_component"] = cc_annots[0] # 2b
            gene_info[g]["cellular_component"] = extracter.cc_annotations(gene_annots)
        gene_info[g]["connections"] = extracter.get_gene_connections(mf_annots, bp_term)

        progress.print_progress()

    return gene_info

def do_stuff():
    start = datetime.datetime.now()

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
        print(genes_and_annots_for_bp(bp))
        
    print("Execution time: " + str(datetime.datetime.now() - start))

# do_stuff()