from ontobio.rdfgen.assoc_rdfgen import CamRdfTransform, TurtleRdfWriter, RdfTransform, genid
from ontobio.vocabulary.relations import OboRO, Evidence
from ontobio.vocabulary.upper import UpperLevel
from prefixcommons.curie_util import expand_uri
from rdflib.namespace import OWL, RDF
from rdflib import Literal
from rdflib.term import URIRef
from rdflib.namespace import Namespace
import rdflib
import logging
from pombase_direct_bp_annots_query import setup_pombase, GOTermAnalyzer
from pombase_golr_query import query_for_annots, genes_and_annots_for_bp, GeneConnectionSet, GeneConnection

# logging.basicConfig(level=logging.INFO)

ro = OboRO()
evt = Evidence()
upt = UpperLevel()
LEGO = Namespace("http://geneontology.org/lego/")
LAYOUT = Namespace("http://geneontology.org/lego/hint/layout/")
PAV = Namespace('http://purl.org/pav/')
DC = Namespace("http://purl.org/dc/elements/1.1/")
onto, aset = setup_pombase()
analyzer = GOTermAnalyzer(onto)
# associations = query_for_annots("PomBase:SPBC11B10.09", "GO:0004693")

# print(gene_info["PomBase:SPAC12B10.10"])


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

def get_aspect(go_term):
    if analyzer.is_molecular_function(go_term):
        return 'F'
    elif analyzer.is_cellular_component(go_term):
        return 'C'
    elif analyzer.is_biological_process(go_term):
        return 'P'

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

class CamTurtleRdfWriter(TurtleRdfWriter):
    def __init__(self, modeltitle):
        self.base = genid(base="http://model.geneontology.org")
        self.graph = rdflib.Graph(identifier=self.base)
        self.graph.bind("owl", OWL)
        self.graph.bind("obo", "http://purl.obolibrary.org/obo/")
        # self.graph.bind("lego", URIRef(LEGO))
        # self.graph.bind("layout", LAYOUT)
        # self.graph.bind("pav", URIRef(PAV))
        self.graph.bind("dc", DC)

        self.graph.add((self.base, RDF.type, OWL.Ontology))

        # Model attributes TODO: Should move outside init
        self.graph.add((self.base, URIRef("http://purl.org/pav/providedBy"), Literal("http://geneontology.org")))        
        self.graph.add((self.base, DC.date, Literal("2018-02-06")))
        self.graph.add((self.base, DC.title, Literal(modeltitle)))
        self.graph.add((self.base, DC.contributor, Literal("http://orcid.org/0000-0002-6659-0416")))
        self.graph.add((self.base, URIRef("http://geneontology.org/lego/modelstate"), Literal("development")))
        self.graph.add((self.base, OWL.versionIRI, self.base))
        self.graph.add((self.base, OWL.imports, URIRef("http://purl.obolibrary.org/obo/go/extensions/go-lego.owl")))

# class AnnotonCamRdfTransform(RdfTransform):
class AnnotonCamRdfTransform(CamRdfTransform):
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
        # enabler_id = genid(base=self.writer.base)
        if sub not in annoton.individuals:
            enabler_id = genid(base=self.writer.base + '/')
            annoton.individuals[sub] = enabler_id
            self.emit_type(enabler_id, sub_uri)
            self.emit_type(enabler_id, OWL.NamedIndividual)
        else:
            enabler_id = annoton.individuals[sub]

        # E.g. instance of GO class
        # tgt_id = genid(base=self.writer.base)
        if obj["id"] not in annoton.individuals:
            tgt_id = genid(base=self.writer.base + '/')
            annoton.individuals[obj["id"]] = tgt_id
            self.emit_type(tgt_id, obj_uri)
            self.emit_type(tgt_id, OWL.NamedIndividual)
        else:
            tgt_id = annoton.individuals[obj["id"]]

        # aspect = association['aspect']
        stmt = None

        # todo: use relation
        # if aspect == 'F':
        #     stmt = self.emit(tgt_id, ENABLED_BY, enabler_id)
        # elif aspect == 'P':
        #     mf_id = genid(base=self.writer.base)
        #     self.emit_type(mf_id, MOLECULAR_FUNCTION)
        #     stmt = self.emit(mf_id, ENABLED_BY, enabler_id)
        #     stmt = self.emit(mf_id, PART_OF, tgt_id)
        # elif aspect == 'C':
        #     mf_id = genid(base=self.writer.base)
        #     self.emit_type(mf_id, MOLECULAR_FUNCTION)
        #     stmt = self.emit(mf_id, ENABLED_BY, enabler_id)
        #     stmt = self.emit(mf_id, OCCURS_IN, tgt_id)

        self.emit_type(tgt_id, obj_uri)
        enabled_by_stmt = self.emit(tgt_id, ENABLED_BY, enabler_id)
        part_of_stmt = self.emit(tgt_id, PART_OF, bp_id)

        self.translate_evidence(annoton.molecular_function, enabled_by_stmt)
        if annoton.cellular_component is not None:
            cc_object_id = annoton.cellular_component["object"]["id"]
            cc_uri = self.uri(cc_object_id)
            if cc_object_id not in annoton.individuals:
                cc_id = genid(base=self.writer.base + '/')
                annoton.individuals[cc_object_id] = cc_id
                self.emit_type(cc_id, cc_uri)
                self.emit_type(cc_id, OWL.NamedIndividual)
            occurs_in_stmt = self.emit(tgt_id, OCCURS_IN, cc_id)
            self.translate_evidence(annoton.cellular_component, occurs_in_stmt)
        
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
        ev = {'type' : association["evidence_type"],
              'has_supporting_reference' : association["reference"]}
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
        for existing_evidence in evidences:
            if evidence["type"] == existing_evidence["type"] and set(evidence["has_supporting_reference"]) == set(existing_evidence["has_supporting_reference"]):
                # print(existing_evidence["id"])
                if "id" not in existing_evidence:
                    existing_evidence["id"] = genid(base=self.writer.base + '/')
                    ev_ids.append(existing_evidence["id"])
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
        evidences.append(evidence)
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

modeltitle = "rdf_output_annotons_queries"
cam_writer = CamTurtleRdfWriter(modeltitle)
writer = AnnotonCamRdfTransform(cam_writer)
classes = []
# individuals = {}
evidences = []
ev_ids = []

# AnnotionProperty
writer.emit_type(URIRef("http://geneontology.org/lego/evidence"), OWL.AnnotationProperty)
writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/x"), OWL.AnnotationProperty)
writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/y"), OWL.AnnotationProperty)
writer.emit_type(URIRef("http://purl.org/pav/providedBy"), OWL.AnnotationProperty)

# associations = []
annotons = []
bp = "GO:0010971"
gene_info = genes_and_annots_for_bp(bp)
bp_id = genid(base=writer.writer.base + '/')
writer.emit_type(bp_id, writer.uri(bp))
writer.emit_type(bp_id, OWL.NamedIndividual)
for gene in gene_info:
    annoton = Annoton(gene_info[gene], gene)
    if annoton.molecular_function is not None:
        annotons.append(annoton)
    # if "molecular_function" in gene_info[gene]:
    #     associations.append(gene_info[gene]["molecular_function"])
    # if "cellular_component" in gene_info[gene]:
    #     associations.append(gene_info[gene]["cellular_component"])

# translate lists of annotations
for annoton in annotons:
    # Class
    if annoton.enabled_by not in classes:
        writer.emit_type(URIRef("http://identifiers.org/" + annoton.enabled_by), OWL.Class)
        classes.append(annoton.enabled_by)

    writer.translate_annoton(annoton)

with open(modeltitle + ".ttl", 'wb') as f:
    writer.writer.serialize(destination=f)