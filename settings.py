configuration = {
    "email": "jon.demasi@colorado.edu",
    "direct_query": False,
}

# Directory to store generated files, point url /wc to this location on your webserver
WORDCLOUD_STORAGE = '/var/tmp/wordclouds'

# Database string to be used by the sqlalchemy ORM
DATABASE_STRING = "sqlite:///db.sqlite3"

# Words to ignore when generating word cloud
FILTER_LIST = [
    'polymorphism', 'polymorphisms', 'nucleotide', 'nucleotides', 'snp', 'snps', 'allele', 'alleles', 'gene',
    'genes', 'genotype', 'genotypes', 'genotyped',  'single', 'singles', 'genetic', 'genetics', 'study', 'studies',
    'variant', 'variants', 'analysis', 'analyses'
]
