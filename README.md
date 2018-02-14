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
With the source GAF filename argument this now frees up the library to create GO-CAM models from any set of GAF's, not just ones pertaining to *S. pombe*. The example GAF can be downloaded from [here](ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.pombase.gz).

## Dependencies
Requires [ontobio](https://github.com/biolink/ontobio).
