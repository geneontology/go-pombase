from ontobio.rdfgen.assoc_rdfgen import CamRdfTransform, TurtleRdfWriter, RdfTransform, genid
from ontobio.vocabulary.relations import OboRO, Evidence
from ontobio.vocabulary.upper import UpperLevel
from ontobio.ontol_factory import OntologyFactory
from ontobio.assoc_factory import AssociationSetFactory
from prefixcommons.curie_util import expand_uri
from rdflib.namespace import OWL, RDF
from rdflib import Literal
from rdflib.term import URIRef
from rdflib.namespace import Namespace
import rdflib
import logging
import argparse
import datetime
import os.path as path
from pombase_direct_bp_annots_query import setup_pombase, GOTermAnalyzer
from pombase_golr_query import query_for_annots, GeneConnectionSet, GeneConnection
# from pombase_golr_query import genes_and_annots_for_bp
from gaf_query import genes_and_annots_for_bp, NoCacheEagerRemoteSparqlOntology

# logging.basicConfig(level=logging.INFO)

ro = OboRO()
evt = Evidence()
upt = UpperLevel()
LEGO = Namespace("http://geneontology.org/lego/")
LAYOUT = Namespace("http://geneontology.org/lego/hint/layout/")
PAV = Namespace('http://purl.org/pav/')
DC = Namespace("http://purl.org/dc/elements/1.1/")

# Stealing a lot of code for this from ontobio.rdfgen:
# https://github.com/biolink/ontobio

HAS_SUPPORTING_REFERENCE = URIRef(expand_uri(evt.has_supporting_reference, cmaps=[evt._prefixmap]))
ENABLED_BY = URIRef(expand_uri(ro.enabled_by))
ENABLES = URIRef(expand_uri(ro.enables))
INVOLVED_IN = URIRef(expand_uri(ro.involved_in))
PART_OF = URIRef(expand_uri(ro.part_of))
OCCURS_IN = URIRef(expand_uri(ro.occurs_in))
COLOCALIZES_WITH = URIRef(expand_uri(ro.colocalizes_with))
MOLECULAR_FUNCTION = URIRef(expand_uri(upt.molecular_function))
REGULATES = URIRef("http://purl.obolibrary.org/obo/RO_0002211")

now = datetime.datetime.now()
onto = NoCacheEagerRemoteSparqlOntology("go")
a_set = AssociationSetFactory().create(onto, file="gene_association.pombase", fmt="gaf")

non_pmid_list = []

def print_triple(s, p, o):
    s_label = a_set.label(s)
    o_label = a_set.label(o)
    if o_label is None:
        o_label = o
    print(s + " (" + s_label + ") " + p + " " + o + " (" + o_label + ")")

class Annoton():
    def __init__(self, gene_info, subject_id):
        self.enabled_by = subject_id
        self.molecular_function = self.get_aspect_object(gene_info, "molecular_function")
        self.cellular_component = self.get_aspect_object(gene_info, "cellular_component")
        self.biological_process = self.get_aspect_object(gene_info, "bp")
        self.connections = gene_info["connections"]
        self.individuals = {}

    def get_aspect_object(self, gene_info, aspect):
        if aspect in gene_info:
            return gene_info[aspect]

class GoCamModel():
    relations_dict = {
        "has_direct_input" : "RO:0002400",
        "has input" : "RO:0002233",
        "has_regulation_target" : "RO:0002211", # regulates
        "regulates_activity_of" : "RO:0002578", # directly regulates
        "with_support_from" : "RO:0002233" # has input
    }

    def __init__(self, filepath, connection_relations=None):
        self.filepath = filepath
        self.modeltitle = path.basename(self.filepath)
        # if self.modeltitle.endswith(".ttl"):
        if path.splitext(self.filepath)[1] != ".ttl":
            self.filepath += ".ttl"
        else:
            self.modeltitle = self.modeltitle[:-4]
        cam_writer = CamTurtleRdfWriter(self.modeltitle)
        self.writer = AnnotonCamRdfTransform(cam_writer)
        self.classes = []
        self.individuals = {}   # Maintain entity-to-IRI dictionary. Prevents dup individuals but we may want dups?
        if connection_relations is None:
            self.connection_relations = GoCamModel.relations_dict
        else:
            self.connection_relations = connection_relations
        self.declare_properties()

    def declare_properties(self):
        # AnnotionProperty
        self.writer.emit_type(URIRef("http://geneontology.org/lego/evidence"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/x"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/y"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/pav/providedBy"), OWL.AnnotationProperty)

    def declare_class(self, class_id):
        self.writer.emit_type(URIRef("http://identifiers.org/" + class_id), OWL.Class)
        self.classes.append(class_id)

    def declare_individual(self, entity_id):
        entity = genid(base=self.writer.writer.base + '/')
        self.writer.emit_type(entity, self.writer.uri(entity_id))
        self.writer.emit_type(entity, OWL.NamedIndividual)
        self.individuals[entity_id] = entity
        return entity

    def add_axiom(self, statement):
        (source_id, property_id, target_id) = statement
        stmt_id = self.find_bnode(statement)
        if stmt_id is None:
            stmt_id = self.writer.blanknode()
            self.writer.emit_type(stmt_id, OWL.Axiom)
        self.writer.emit(stmt_id, OWL.annotatedSource, source_id)
        self.writer.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.writer.emit(stmt_id, OWL.annotatedTarget, target_id)
        return stmt_id

    def add_evidence(self, axiom, evidence_code, references):
        ev = {'type' : evidence_code,
            'has_supporting_reference' : references}
        # Try finding existing evidence object containing same type and references
        ev_id = self.writer.find_or_create_evidence_id(ev)
        self.writer.emit(axiom, URIRef("http://geneontology.org/lego/evidence"), ev_id)

    def add_connection(self, gene_connection, source_annoton):
        # Switching from reusing existing activity node from annoton to creating new one for each connection - Maybe SPARQL first to check if annoton activity already used for connection?
        # Check annoton for existing activity.
        # if gene_connection.object_id in source_annoton.individuals:
        #     # If exists and activity has connection relation,
        #     # Look for two triples: (gene_connection.object_id, ENABLED_BY, source_annoton.enabled_by) and (gene_connection.object_id, connection_relations, anything)
        source_id = None
        uri_list = self.uri_list_for_individual(gene_connection.object_id)
        o_count = 0
        for u in uri_list:
            if gene_connection.relation in self.connection_relations:
                rel = URIRef(expand_uri(self.connection_relations[gene_connection.relation]))
                # Annot MF should be declared by now - don't declare object_id if object_id == annot MF?
                try:
                    annot_mf = source_annoton.molecular_function["object"]["id"]
                except:
                    annot_mf = ""
                if self.writer.writer.graph.__contains__((u,rel,None)) and gene_connection.object_id != annot_mf:
                    source_id = self.declare_individual(gene_connection.object_id)
                    source_annoton.individuals[gene_connection.object_id] = source_id
                    break
                # for t in self.writer.writer.graph.triples((u,ENABLED_BY,None)):

                # if t[1] == OCCURS_IN:
                #     print(gene_connection.object_id + " OCCURS_IN " + str(t[2]))
        if source_id is None:
            try:
                source_id = source_annoton.individuals[gene_connection.object_id]
            except KeyError:
                source_id = self.declare_individual(gene_connection.object_id)
                source_annoton.individuals[gene_connection.object_id] = source_id
        # Add enabled by stmt for object_id - this is essentially adding another annoton connecting gene-to-extension/with-MF to the model
        self.writer.emit(source_id, ENABLED_BY, source_annoton.individuals[source_annoton.enabled_by])
        self.writer.emit_axiom(source_id, ENABLED_BY, source_annoton.individuals[source_annoton.enabled_by])
        property_id = URIRef(expand_uri(self.connection_relations[gene_connection.relation]))
        target_id = self.individuals[gene_connection.gp_b]
        # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
        self.writer.emit(source_id, property_id, target_id)
        # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
        self.writer.emit_axiom(source_id, property_id, target_id)

    def uri_list_for_individual(self, individual):
        uri_list = []
        graph = self.writer.writer.graph
        for t in graph.triples((None,None,self.writer.uri(individual))):
            uri_list.append(t[0])
        return uri_list

    def find_bnode(self, triple):
        (subject,predicate,object_id) = triple
        s_triples = self.writer.writer.graph.triples((None, OWL.annotatedSource, subject))
        s_bnodes = [s for s,p,o in s_triples]
        p_triples = self.writer.writer.graph.triples((None, OWL.annotatedProperty, predicate))
        p_bnodes = [s for s,p,o in p_triples]
        o_triples = self.writer.writer.graph.triples((None, OWL.annotatedTarget, object_id))
        o_bnodes = [s for s,p,o in o_triples]
        bnodes = set(s_bnodes) & set(p_bnodes) & set(o_bnodes)
        if len(bnodes) > 0:
            return list(bnodes)[0]

class CamTurtleRdfWriter(TurtleRdfWriter):
    def __init__(self, modeltitle):
        self.base = genid(base="http://model.geneontology.org")
        self.graph = rdflib.Graph(identifier=self.base)
        self.graph.bind("owl", OWL)
        self.graph.bind("obo", "http://purl.obolibrary.org/obo/")
        self.graph.bind("dc", DC)

        self.graph.add((self.base, RDF.type, OWL.Ontology))

        # Model attributes TODO: Should move outside init
        self.graph.add((self.base, URIRef("http://purl.org/pav/providedBy"), Literal("http://geneontology.org")))        
        self.graph.add((self.base, DC.date, Literal(str(now.year) + "-" + str(now.month) + "-" + str(now.day))))
        self.graph.add((self.base, DC.title, Literal(modeltitle)))
        self.graph.add((self.base, DC.contributor, Literal("http://orcid.org/0000-0002-6659-0416"))) #TODO
        self.graph.add((self.base, URIRef("http://geneontology.org/lego/modelstate"), Literal("development")))
        self.graph.add((self.base, OWL.versionIRI, self.base))
        self.graph.add((self.base, OWL.imports, URIRef("http://purl.obolibrary.org/obo/go/extensions/go-lego.owl")))

class AnnotonCamRdfTransform(CamRdfTransform):
    def __init__(self, writer=None):
        CamRdfTransform.__init__(self, writer)
        self.annotons = []
        self.classes = []
        self.evidences = []
        self.ev_ids = []
        self.bp_id = None

    def translate_annoton(self, annoton, and_xps=None):
        # See https://github.com/biolink/ontobio/pull/136
        # if the association has an annotation extension, and this
        # is a union, then we treat each element in the union
        # as a distinct assertion/annotation, where each assertion
        # has its own conjunction of relational expressions
        
        # if and_xps is None and 'object_extensions' in association:
        #     x = association['object_extensions']
        #     for ix in x['union_of']:
        #         and_xps = ix['intersection_of']
        #         self.translate(association, and_xps)
            
        # sub = association['subject']
        # obj = association['object']
        # rel = association['relation']
        sub = annoton.enabled_by
        obj = annoton.molecular_function["object"]
        sub_uri = self.uri(sub)
        obj_uri = self.uri(obj)

        # E.g. instance of gene product class
        if sub not in annoton.individuals:
            enabler_id = genid(base=self.writer.base + '/')
            annoton.individuals[sub] = enabler_id
            self.emit_type(enabler_id, sub_uri)
            self.emit_type(enabler_id, OWL.NamedIndividual)
        else:
            enabler_id = annoton.individuals[sub]

        # E.g. instance of GO class
        if obj["id"] not in annoton.individuals:
            tgt_id = genid(base=self.writer.base + '/')
            annoton.individuals[obj["id"]] = tgt_id
            self.emit_type(tgt_id, obj_uri)
            self.emit_type(tgt_id, OWL.NamedIndividual)
        else:
            tgt_id = annoton.individuals[obj["id"]]

        stmt = None

        self.emit_type(tgt_id, obj_uri)
        enabled_by_stmt = self.emit(tgt_id, ENABLED_BY, enabler_id)
        part_of_stmt = self.emit(tgt_id, PART_OF, self.bp_id)

        self.translate_evidence(annoton.molecular_function, enabled_by_stmt)
        if annoton.cellular_component is not None:
            # cc_object_id = annoton.cellular_component["object"]["id"]
            # cc_uri = self.uri(cc_object_id)
            # if cc_object_id not in annoton.individuals:
            #     cc_id = genid(base=self.writer.base + '/')
            #     annoton.individuals[cc_object_id] = cc_id
            #     self.emit_type(cc_id, cc_uri)
            #     self.emit_type(cc_id, OWL.NamedIndividual)
            # occurs_in_stmt = self.emit(tgt_id, OCCURS_IN, cc_id)
            # self.translate_evidence(annoton.cellular_component, occurs_in_stmt)
            for cellular_component in annoton.cellular_component:
                cc_object_id = cellular_component["object"]["id"]
                cc_uri = self.uri(cc_object_id)
                if cc_object_id not in annoton.individuals:
                    cc_id = genid(base=self.writer.base + '/')
                    annoton.individuals[cc_object_id] = cc_id
                    self.emit_type(cc_id, cc_uri)
                    self.emit_type(cc_id, OWL.NamedIndividual)
                occurs_in_stmt = self.emit(tgt_id, OCCURS_IN, cc_id)
                self.translate_evidence(cellular_component, occurs_in_stmt)
        
        if and_xps is not None:
            for ext in and_xps:
                filler_inst = genid(base=self.writer.base)
                self.emit_type(filler_inst, self.uri(ext['filler']))
                p = self.lookup_relation(ext['property'])
                if p is None:
                    logging.warning("No such property {}".format(ext))
                else:
                    self.emit(tgt_id, p, filler_inst)

    def translate_evidence(self, association, stmt):
        """

        ``
        _:1 a Axiom
            owl:annotatedSource s
            owl:annotatedProperty p
            owl:annotatedTarget o
            evidence [ a ECO ; ...]
        ``

        """
        #TODO make smarter Annotation class
        if "evidence_type" in association:
            ev = {'type' : association["evidence_type"],
                'has_supporting_reference' : association["reference"]}
        else:
            ev = association["evidence"]
        # Try finding existing evidence object containing same type and references
        ev_id = self.find_or_create_evidence_id(ev)
        # ev_id = None
        # if 'id' in ev:
        #     ev_id = self.uri(ev['id'])
        # else:
        #     ev_id = genid(base=self.writer.base + '/')

        # stmt_id = self.blanknode() ## OWL reification: must be blank
        (s,p,o) = stmt
        stmt_id = self.find_bnode(stmt)
        if stmt_id is None:
            stmt_id = self.blanknode() ## OWL reification: must be blank
            self.emit_type(stmt_id, OWL.Axiom)

        self.emit(stmt_id, OWL.annotatedSource, s)
        self.emit(stmt_id, OWL.annotatedProperty, p)
        self.emit(stmt_id, OWL.annotatedTarget, o)

        ev_triple = (stmt_id, URIRef("http://geneontology.org/lego/evidence"), ev_id)
        self.emit(stmt_id, URIRef("http://geneontology.org/lego/evidence"), ev_id)
    
    def find_or_create_evidence_id(self, evidence):
        for existing_evidence in self.evidences:
            if evidence["type"] == existing_evidence["type"] and set(evidence["has_supporting_reference"]) == set(existing_evidence["has_supporting_reference"]):
                if "id" not in existing_evidence:
                    existing_evidence["id"] = genid(base=self.writer.base + '/')
                    self.ev_ids.append(existing_evidence["id"])
                return existing_evidence["id"]
        ev_id = genid(base=self.writer.base + '/')
        evidence["id"] = ev_id
        ev_cls = self.eco_class(self.uri(evidence['type']))
        self.emit_type(ev_id, OWL.NamedIndividual)
        self.emit_type(ev_id, ev_cls)
        if 'with_support_from' in evidence:
            for w in evidence['with_support_from']:
                self.emit(ev_id, self.uri(evt.evidence_with_support_from), self.uri(w))
        for ref in evidence['has_supporting_reference']:
            # o = self.uri(ref)
            o = Literal(ref) # Needs to go into Noctua like 'PMID:####' rather than full URL
            # if ref == expand_uri(ref):
            #     o = Literal(ref)
            self.emit(ev_id, HAS_SUPPORTING_REFERENCE, o)
        if 'with_support_from' in evidence:
            for ref in evidence['with_support_from']:
                self.emit(ev_id, self.uri(evt.evidence_with_support_from), self.uri(ref))
        self.evidences.append(evidence)
        return evidence["id"]

    def find_bnode(self, triple):
        (subject,predicate,object_id) = triple
        s_triples = self.writer.graph.triples((None, OWL.annotatedSource, subject))
        s_bnodes = [s for s,p,o in s_triples]
        p_triples = self.writer.graph.triples((None, OWL.annotatedProperty, predicate))
        p_bnodes = [s for s,p,o in p_triples]
        o_triples = self.writer.graph.triples((None, OWL.annotatedTarget, object_id))
        o_bnodes = [s for s,p,o in o_triples]
        bnodes = set(s_bnodes) & set(p_bnodes) & set(o_bnodes)
        if len(bnodes) > 0:
            return list(bnodes)[0]

    def emit_axiom(self, source_id, property_id, target_id):
        stmt_id = self.blanknode()
        self.emit_type(stmt_id, OWL.Axiom)
        self.emit(stmt_id, OWL.annotatedSource, source_id)
        self.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.emit(stmt_id, OWL.annotatedTarget, target_id)
        return stmt_id

    def find_annotons(self, enabled_by, annotons_list=None):
        found_annotons = []
        if annotons_list is not None:
            annotons = annotons_list
        else:
            annotons = self.annotons
        for annoton in annotons:
            if annoton.enabled_by == enabled_by:
                found_annotons.append(annoton)
        return found_annotons

    def add_individual(self, individual_id, annoton):
        obj_uri = self.uri(individual_id)
        if individual_id not in annoton.individuals:
            tgt_id = genid(base=self.writer.base + '/')
            annoton.individuals[individual_id] = tgt_id
            self.emit_type(tgt_id, obj_uri)
            self.emit_type(tgt_id, OWL.NamedIndividual)
        else:
            tgt_id = annoton.individuals[individual_id]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--bp_term", type=str, required=True,
                        help="Biological process GO term that GOCAM should model")
    parser.add_argument('-f', "--filename", type=str, required=False,
                        help="Destination filename - will end in '.ttl'")
    parser.add_argument('-g', "--gaf_source", type=str, required=True,
                        help="filename of GAF file to use as annotation source")
    parser.add_argument('-j', "--tad_json", type=str, required=False,
                        help="Existing json data file to load into term-to-gene dictionary. Speeds up performance.")

    args = parser.parse_args()

    bp = args.bp_term
    gene_info = genes_and_annots_for_bp(bp, args.gaf_source, json_file=args.tad_json)

    modeltitle = args.filename
    if modeltitle.endswith(".ttl"):
        modeltitle = modeltitle[:-4]
    cam_writer = CamTurtleRdfWriter(modeltitle)
    writer = AnnotonCamRdfTransform(cam_writer)
    writer.bp_id = genid(base=writer.writer.base + '/')
    writer.emit_type(writer.bp_id, writer.uri(bp))
    writer.emit_type(writer.bp_id, OWL.NamedIndividual)

    # AnnotionProperty
    writer.emit_type(URIRef("http://geneontology.org/lego/evidence"), OWL.AnnotationProperty)
    writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/x"), OWL.AnnotationProperty)
    writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/y"), OWL.AnnotationProperty)
    writer.emit_type(URIRef("http://purl.org/pav/providedBy"), OWL.AnnotationProperty)

    for gene in gene_info:
        annoton = Annoton(gene_info[gene], gene)
        if annoton.molecular_function is not None:
            writer.annotons.append(annoton)

    global_individuals_list = {}

    # translate lists of annotations
    for annoton in writer.annotons:
        # Class
        if annoton.enabled_by not in writer.classes:
            writer.emit_type(URIRef("http://identifiers.org/" + annoton.enabled_by), OWL.Class)
            writer.classes.append(annoton.enabled_by)

        writer.translate_annoton(annoton)

        global_individuals_list = {**global_individuals_list, **annoton.individuals}

    ### Connections - Now that all individuals should have been created
    connection_relations = {"has_direct_input" : "RO:0002400",
                            "has input" : "RO:0002233",
                            "has_regulation_target" : "RO:0002211", # regulates
                            "regulates_activity_of" : "RO:0002578", # directly regulates
                            "with_support_from" : "RO:0002233" # has input
                            }
    for annoton in writer.annotons:
        for connection in annoton.connections.gene_connections:
            if connection.relation in ["has_direct_input", "has input", "with_support_from"]:
                try:
                    source_id = annoton.individuals[connection.object_id]
                except KeyError:
                    writer.add_individual(connection.object_id, annoton)
                    source_id = annoton.individuals[connection.object_id]
                object_source_id = annoton.individuals[annoton.molecular_function["object"]["id"]]
                # Add enabled by stmt for object_id
                e_stmt = writer.emit(object_source_id, ENABLED_BY, annoton.individuals[annoton.enabled_by])
                writer.emit_axiom(object_source_id, ENABLED_BY, annoton.individuals[annoton.enabled_by])
                # writer.translate_evidence(connection.annotation, e_stmt)
                property_id = URIRef(expand_uri(connection_relations[connection.relation]))
                target_id = global_individuals_list[connection.gp_b]
                # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                writer.emit(source_id, property_id, target_id)
                # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                writer.emit_axiom(source_id, property_id, target_id)
            elif connection.relation in ["has_regulation_target", "regulates_activity_of"]:
                source_id = annoton.individuals[connection.gp_a]
                property_id = URIRef(expand_uri(connection_relations[connection.relation]))
                # find annoton(s) of regulation target gene product
                target_annotons = writer.find_annotons(connection.gp_b)
                for t_annoton in target_annotons:
                    mf_annotation = t_annoton.get_aspect_object(gene_info[connection.gp_b], "molecular_function")
                    if mf_annotation is not None:
                        target_id = global_individuals_list[mf_annotation["object"]["id"]]
                        # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                        writer.emit(source_id, property_id, target_id)
                        # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                        writer.emit_axiom(source_id, property_id, target_id)

    with open(modeltitle + ".ttl", 'wb') as f:
        writer.writer.serialize(destination=f)

if __name__ == "__main__":
    main()