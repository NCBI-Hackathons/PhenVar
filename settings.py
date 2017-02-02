configuration = {
    "email": "jon.demasi@colorado.edu",
    "direct_query": False,
}

# Directory to store generated files
WORDCLOUD_STORAGE = 'static/wc'

# Database string to be used by the sqlalchemy ORM
DATABASE_STRING = "sqlite:///db.sqlite3"

# Words to ignore when generating word cloud
FILTER_LIST = [
    'polymorphism', 'polymorphisms', 'nucleotide', 'nucleotides', 'snp', 'snps', 'allele', 'alleles', 'gene',
    'genes', 'genotype', 'genotypes', 'genotyped',  'single', 'singles', 'genetic', 'genetics', 'study', 'studies',
    'variant', 'variants', 'analysis', 'analyses''association', 'study', 'SNPs', 'gene', 'genes', 'polymorphisms',
    'risk', 'studies', 'variants', 'analysis', 'results', 'disease', 'associations', 'SNP', 'controls', 'patients',
    'data', 'susceptibility', 'cases', 'population', 'evidence', 'role', 'loci', 'individuals', 'variation', 'findings',
    'factors', 'allele', 'candidate', 'region', 'effects', 'effect', 'polymorphism', 'analyses', 'levels', 'cancer',
    'genotype', 'locus', 'expression', 'linkage', 'age', 'subjects', 'samples', 'CI', 'regions', 'populations', 'OR',
    'sample', 'cohort', 'development', 'factor', 'odds', 'number', 'GWAS', 'protein', 'chromosome', 'alleles', 'DNA',
    'replication', 'receptor', 'regression', 'significance', 'control', 'model', 'genotypes', 'haplotype', 'frequency',
    'variant', 'markers', 'confidence', 'phenotypes', 'Study', 'traits', 'diseases', 'differences', 'approach', 'years',
    'function', 'response', 'case-control', 'total', 'models', 'cell', 'mutations', 'families', 'type', 'addition',
    'women', 'blood', 'disorder', 'interval', 'interaction', 'family', 'ratio', 'disequilibrium', 'pathways',
    'interactions', 'trait', 'pathway', 'test', 'groups', 'phenotype', 'pathogenesis', 'rate', 'level', 'group',
    'haplotypes', 'meta-analysis', 'set', 'participants', 'frequencies', 'treatment', 'discovery', 'diabetes',
    'tests', 'variations', 'mutation', 'cohorts', 'nucleotide', 'method', 'mechanisms', 'disorders', 'cells',
    'activity', 'status', 'testing', 'metabolism', 'identification', '×', 'sequence', 'index', 'mass', 'ancestry',
    'information', 'system', 'basis', 'variability', 'regulation', 'influence', 'score', 'correction', 'body',
    'methods', 'syndrome', 'aim', 'genotyping', 'Genome-wide', 'power', 'use', 'presence', 'etiology', 'sex',
    'promoter', 'hypothesis', 'marker', 'Affymetrix', 'ratios', 'relationship', 'variance', 'growth', 'research', 'time'
]
