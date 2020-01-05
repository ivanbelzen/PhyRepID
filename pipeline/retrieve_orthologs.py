## IvB April 2018
## Query to get all orthologs between human and selected species
## Add all to dictionary 
## TODO: implement automatic retry if sparql query fails

from SPARQLWrapper import SPARQLWrapper, JSON
import json
import pipeline_methods as pre
	
orthologs_output_file = pre.orthologs_file

reference_species = 9606 #human

sparql=SPARQLWrapper('https://www.ebi.ac.uk/rdf/services/sparql')
sparql.setReturnFormat(JSON)

mapping = {} #for each gene, save species - ortholog combination, if multiple append to list.

if pre.file_notempty(orthologs_output_file):
	print("Orthologs already downloaded, delete orthologs.json to redo")
	quit()


species_list = pre.get_taxon_list('', [reference_species])	
species_str = 'taxon:'+' taxon:'.join([str(x) for x in species_list])

for sp in species_list:
	
	sparql.setQuery(
	"""
		PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
		PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
		PREFIX owl: <http://www.w3.org/2002/07/owl#>
		PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
		PREFIX foaf: <http://xmlns.com/foaf/0.1/>
		PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
		PREFIX obo: <http://purl.obolibrary.org/obo/>
		PREFIX sio: <http://semanticscience.org/resource/>
		PREFIX ensembl: <http://rdf.ebi.ac.uk/resource/ensembl/>
		PREFIX taxon: <http://identifiers.org/taxonomy/>

		SELECT DISTINCT ?gene ?ortholog ?orthologTaxon {
		
		?gene a obo:SO_0001217 .
		?gene sio:SIO_000558 ?ortholog .
		?gene obo:RO_0002162 ?taxon .
		?ortholog obo:RO_0002162 ?orthologTaxon .
		VALUES ?taxon {  taxon:"""+str(reference_species)+""" } .
		VALUES ?orthologTaxon { taxon:"""+str(sp)+""" }
		FILTER (?taxon != ?orthologTaxon)
		} ORDER BY ?gene ?ortholog ?orthologTaxon

		""")
		
		
	results = sparql.query().convert()
	for result in results['results']['bindings']:
		taxon = result["orthologTaxon"]["value"]
		gene_uri = result["gene"]["value"]
		ortholog_uri = result["ortholog"]["value"]
		if gene_uri not in mapping: mapping[gene_uri] = {}
		if taxon not in mapping[gene_uri]: mapping[gene_uri][taxon] = []
		mapping[gene_uri][taxon].append(ortholog_uri)
			
# for testing purposes; VALUES ?gene { ensembl:ENSG00000162692 ensembl:ENSG00000117399 } .

#print(json.dumps(mapping))
with open(orthologs_output_file, 'w') as output:
	output.write(json.dumps(mapping))

