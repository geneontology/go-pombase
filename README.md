# go-pombase
Code for querying GO annotations for PomBase and creating GO-CAM models using ontobio API.

Working towards creating GO-CAMS by inputting a list of annotations (either GAF or directly ontobio association objects). Right now, this takes a GO biologcal process term as input, does some heuristic gene set calculation and generates a GO-CAM ttl for the BP term’s gene set. Separating this gene set logic from the annotation-to-GO-CAM logic is another goal.

## Running
```
pip install -r requirements.txt
```
As this is coded right now for a specific use case, this can be ran simply by inputting a GO BP term, source GAF filename, and a destination filename:
```
python3 generate_rdf.py -t "GO:0010971" -g "gene_association.pombase" -f "filename.ttl"
```
With the source GAF filename argument this now frees up the library to create GO-CAM models from any set of GAF's, not just ones pertaining to *S. pombe*. The example GAF can be downloaded from ftp://ftp.geneontology.org/pub/go/gene-associations/.

## Running for generating PomBase GO-CAM models
For my purpose right now I'm running generate_pombase_model.py specifying BP term (-t), output filename (-f), and GAF input file (-g):
```
python3 generate_pombase_model.py -t 'GO:0031929' -f 'TOR signaling.ttl' -g 'gene_association.pombase'
```
### Reusing computed gene-to-BP term dictionary data
You can also specify the data (-j) to use in the first step in order to speed up processing during repeated runs (~1.5 min -> 10 sec):
```
python3 generate_pombase_model.py -j 'tad_go_gafs.json' -t 'GO:0031929' -f 'TOR signaling.ttl' -g 'gene_association.pombase'
```
To dump out this data into a reusable JSON, you can run:
```
python3 pombase_direct_bp_annots_query.py -j 'json_outfile.json' -g 'gene_association.pombase'
```
With -j specifying the JSON output path.

## Generate BP cluster lists from 'candidate BP terms'
Extract lists of BP terms that share commonly annotated genes
```
python3 pombase_direct_bp_annots_query.py -t candidate_bp_terms.txt -c term_clusters.txt -s unclustered_terms.tsv -g pombase.gaf
```
### Criteria for selecting candidate BP terms
1. Term gene set size is less than or equal to `x`=60
2. One of these is true:
    1. Term gene set size between `m`=5 and `n`=30 AND the BP term has a parent with gene set size > `n`=30
    2. Term gene set size > `n`=30 AND term has no child terms
3. Term has no ancestor terms that also meet above criteria (selects for most generic term)

Note that GO graph traversal here only follows "is_a" paths. Also, the gene set for each term contains genes both directly and indirectly annotated to that term.

## Dependencies
Requires [ontobio](https://github.com/biolink/ontobio).
