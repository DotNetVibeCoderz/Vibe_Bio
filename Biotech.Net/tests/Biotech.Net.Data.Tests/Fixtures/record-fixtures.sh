#!/usr/bin/env bash
# Records the real service responses that Biotech.Net.Data.Tests replays offline.
#
# Run from this directory:  bash record-fixtures.sh
# Needs curl and network access. Set NCBI_EMAIL (and optionally NCBI_API_KEY) to identify yourself to NCBI;
# neither value is written into a fixture, because the e-utilities do not echo them back.
#
# Notes
# - Each response body is saved verbatim. Where a test needs response headers (UniProt's Link pagination
#   header), curl also writes "<name>.headers"; the test handler replays the Link / X-Total-Results lines.
# - Files are kept small (< 100 KB): a reviewed but lightly annotated UniProt entry (crambin, P01542) instead
#   of a heavily annotated one, field-restricted UniProt searches, and gzip PDB/mmCIF downloads.
# - The UniProt ID-mapping job id changes on every run; the tests match the job endpoints by path prefix.
# - Recorded responses reflect the databases on the recording date. Tests assert on stable facts
#   (accessions, lengths, sequences), not on counts that grow with each database release.
set -euo pipefail
cd "$(dirname "$0")"

UA="Biotech.Net-fixture-recorder/0.1 (+https://github.com/DotNetVibeCoderz/Vibe_Bio)"
NCBI="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EXTRA=""
if [ -n "${NCBI_EMAIL:-}" ]; then EXTRA="$EXTRA&email=${NCBI_EMAIL}"; fi
if [ -n "${NCBI_API_KEY:-}" ]; then EXTRA="$EXTRA&api_key=${NCBI_API_KEY}"; fi
EXTRA="$EXTRA&tool=Biotech.Net"

get() { # get <file> <url> [extra curl args]
  local out="$1" url="$2"; shift 2
  curl -sS --fail-with-body --retry 4 --retry-all-errors --retry-delay 2 -A "$UA" -o "$out" "$@" "$url"
  echo "$(wc -c < "$out") bytes  $out"
  sleep 0.4   # stay below NCBI's 3 requests/second without a key
}

mkdir -p ncbi uniprot pdb ensembl

# ---------------------------------------------------------------- NCBI E-utilities
get ncbi/einfo.json            "$NCBI/einfo.fcgi?retmode=json$EXTRA"
get ncbi/einfo-pubmed.json     "$NCBI/einfo.fcgi?db=pubmed&retmode=json$EXTRA"
get ncbi/esearch-nuccore.json  "$NCBI/esearch.fcgi?db=nuccore&term=HBB%5Bgene%5D+AND+Homo+sapiens%5Borgn%5D+AND+refseq%5Bfilter%5D&retmax=5&retmode=json$EXTRA"
get ncbi/esearch-pubmed.json   "$NCBI/esearch.fcgi?db=pubmed&term=crambin+structure&retmax=3&sort=relevance&retmode=json$EXTRA"
get ncbi/esummary-pubmed.json  "$NCBI/esummary.fcgi?db=pubmed&id=6895315,19562114&retmode=json$EXTRA"
get ncbi/efetch-protein.fasta  "$NCBI/efetch.fcgi?db=protein&id=NP_000509.1&rettype=fasta&retmode=text$EXTRA"
get ncbi/efetch-nuccore.gb     "$NCBI/efetch.fcgi?db=nuccore&id=NM_000518.5&rettype=gb&retmode=text$EXTRA"
get ncbi/elink-protein-gene.json "$NCBI/elink.fcgi?dbfrom=protein&db=gene&id=NP_000509.1&retmode=json$EXTRA"
# History server: one search kept on the server, then two pages of FASTA fetched with WebEnv + query_key.
get ncbi/esearch-history.json  "$NCBI/esearch.fcgi?db=protein&term=(HBA1%5Bgene%5D+OR+HBB%5Bgene%5D+OR+HBD%5Bgene%5D)+AND+Homo+sapiens%5Borgn%5D+AND+refseq%5Bfilter%5D&usehistory=y&retmax=0&retmode=json$EXTRA"
WEBENV=$(sed -n 's/.*"webenv":"\([^"]*\)".*/\1/p' ncbi/esearch-history.json)
QK=$(sed -n 's/.*"querykey":"\([^"]*\)".*/\1/p' ncbi/esearch-history.json)
get ncbi/efetch-history-0.fasta "$NCBI/efetch.fcgi?db=protein&query_key=$QK&WebEnv=$WEBENV&retstart=0&retmax=2&rettype=fasta&retmode=text$EXTRA"
get ncbi/efetch-history-2.fasta "$NCBI/efetch.fcgi?db=protein&query_key=$QK&WebEnv=$WEBENV&retstart=2&retmax=2&rettype=fasta&retmode=text$EXTRA"

# ---------------------------------------------------------------- UniProt
UP="https://rest.uniprot.org"
get uniprot/P01542.json  "$UP/uniprotkb/P01542?format=json"
get uniprot/P01542.fasta "$UP/uniprotkb/P01542?format=fasta"
FIELDS="accession,id,protein_name,gene_names,organism_name,length,sequence"
get uniprot/search-page1.json "$UP/uniprotkb/search?query=gene:HBB+AND+reviewed:true&fields=$FIELDS&format=json&size=2" -D uniprot/search-page1.headers
NEXT=$(sed -n 's/^[Ll]ink: <\([^>]*\)>; rel="next".*/\1/p' uniprot/search-page1.headers | tr -d '\r')
get uniprot/search-page2.json "$NEXT" -D uniprot/search-page2.headers
get uniprot/idmapping-run.json "$UP/idmapping/run" -X POST --data "from=UniProtKB_AC-ID&to=GeneID&ids=P69905,P68871,P01542"
JOB=$(sed -n 's/.*"jobId":"\([^"]*\)".*/\1/p' uniprot/idmapping-run.json)
curl -sS -A "$UA" -o uniprot/idmapping-status-first.json "$UP/idmapping/status/$JOB"; echo "first status: $(cat uniprot/idmapping-status-first.json | head -c 200)"
sleep 3
get uniprot/idmapping-results.json "$UP/idmapping/results/$JOB"

# ---------------------------------------------------------------- RCSB PDB + AlphaFold DB
get pdb/1CRN.pdb.gz       "https://files.rcsb.org/download/1CRN.pdb.gz"
get pdb/1CRN.cif.gz       "https://files.rcsb.org/download/1CRN.cif.gz"
get pdb/entry-1CRN.json   "https://data.rcsb.org/rest/v1/core/entry/1CRN"
get pdb/polymer-entity-1CRN-1.json "https://data.rcsb.org/rest/v1/core/polymer_entity/1CRN/1"
get pdb/graphql-1CRN.json "https://data.rcsb.org/graphql" -X POST -H "Content-Type: application/json" \
  --data '{"query":"{ entry(entry_id: \"1CRN\") { rcsb_id struct { title } exptl { method } } }"}'
get pdb/search-text.json  "https://search.rcsb.org/rcsbsearch/v2/query" -X POST -H "Content-Type: application/json" \
  --data '{"query":{"type":"terminal","service":"full_text","parameters":{"value":"crambin"}},"return_type":"entry","request_options":{"paginate":{"start":0,"rows":5},"results_content_type":["experimental"]}}'
get pdb/search-sequence.json "https://search.rcsb.org/rcsbsearch/v2/query" -X POST -H "Content-Type: application/json" \
  --data '{"query":{"type":"terminal","service":"sequence","parameters":{"value":"TTCCPSIVARSNFNVCRLPGTPEALCATYTGCIIIPGATCPGDYAN","sequence_type":"protein","evalue_cutoff":0.1,"identity_cutoff":0.9}},"return_type":"polymer_entity","request_options":{"paginate":{"start":0,"rows":5},"results_content_type":["experimental"]}}'
get pdb/alphafold-P01542.json "https://alphafold.ebi.ac.uk/api/prediction/P01542"
AFPDB=$(sed -n 's/.*"pdbUrl":"\([^"]*\)".*/\1/p' pdb/alphafold-P01542.json)
get pdb/AF-P01542-F1-model.pdb "$AFPDB"

# ---------------------------------------------------------------- Ensembl REST
EN="https://rest.ensembl.org"
J="content-type=application/json"
get ensembl/lookup-symbol-BRCA2.json "$EN/lookup/symbol/homo_sapiens/BRCA2?$J"
get ensembl/lookup-id-HBB-expand.json "$EN/lookup/id/ENSG00000244734?expand=1;$J"
get ensembl/sequence-protein.json "$EN/sequence/id/ENST00000335295?type=protein;$J"
get ensembl/sequence-cds.json "$EN/sequence/id/ENST00000335295?type=cds;$J"
get ensembl/sequence-region.json "$EN/sequence/region/human/11:5227002..5227071:-1?$J"
get ensembl/overlap-region.json "$EN/overlap/region/human/11:5225000-5230000?feature=gene;feature=transcript;$J"
get ensembl/vep-hgvs.json "$EN/vep/human/hgvs/ENST00000335295:c.20A%3ET?$J"
get ensembl/vep-region.json "$EN/vep/human/region/11:5227002-5227002:1/A?$J"
get ensembl/xrefs-id.json "$EN/xrefs/id/ENSG00000244734?external_db=HGNC;$J"
get ensembl/xrefs-symbol.json "$EN/xrefs/symbol/homo_sapiens/HBB?$J"
