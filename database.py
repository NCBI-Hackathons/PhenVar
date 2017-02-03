from sqlalchemy import create_engine, Column, Integer, String, ForeignKey, Date
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
from settings import DATABASE_STRING, FILTER_LIST, MAPPING_FILE
from ncbiutils import get_pubmed_articles
import json
import nltk
import progressbar


engine = create_engine(
    DATABASE_STRING,
    connect_args={'check_same_thread': False},
)
Base = declarative_base()
Session = sessionmaker(bind=engine)
session = Session()


# Define tables
class Article(Base):
    __tablename__ = 'articles'
    pmid = Column(Integer, primary_key=True)
    abstract = Column(String)
    date_created = Column(Date)
    date_revised = Column(Date)

    # Processes abstract with nltk library and returns list of nouns
    def abstract_nouns(self):
        # Tokenize and tag abstracts
        tokens = nltk.word_tokenize(self.abstract.lower())
        tagged = nltk.pos_tag(tokens)
        #tokens = tokenizer.tokenize(self.abstract)
        #tagged = tagger.tag(tokens)
        nouns = [tagged_word[0] for tagged_word in tagged if tagged_word[1] in ("NN", "NNS", "NNP", "NNPS")
                     and tagged_word[0] not in FILTER_LIST
                     and len(tagged_word[0]) > 1]
        return nouns

    def unique_abstract_nouns(self):
        nouns = self.abstract_nouns()
        return list(set(nouns))

    # Returns dictionary of form {'noun': count, ...} for each noun in abstract
    def abstract_wordcount(self):
        wordcount_dict = {}
        nouns = self.abstract_nouns()
        unique_nouns = self.unique_abstract_nouns()
        for noun in unique_nouns:
            wordcount_dict[noun] = nouns.count(noun)
        return wordcount_dict


# Author can be either normal person, or collective (type = p or c respectively)
class Author(Base):
    __tablename__ = 'authors'
    id = Column(Integer, autoincrement=True, primary_key=True)
    last_name = Column(String, nullable=True)
    first_name = Column(String, nullable=True)
    initials = Column(String, nullable=True)
    collective_name = Column(String, nullable=True)


# Authors are linked to pubmed articles via author_citations table
class AuthorCitation(Base):
    __tablename__ = 'author_citations'
    author_id = Column(Integer, ForeignKey('authors.id'), primary_key=True)
    pmid = Column(Integer, ForeignKey('articles.pmid'), primary_key=True)


# Affiliation is used to store general affiliation by primary key to avoid storing lengthy strings excessively
class Affiliation(Base):
    __tablename__ = 'affiliations'
    id = Column(Integer, autoincrement=True, primary_key=True)
    info = Column(String)


# Authors are linked to affiliation via id
class AuthorAffiliation(Base):
    __tablename__ = 'author_affiliations'
    author_id = Column(Integer, ForeignKey('authors.id'), primary_key=True)
    affiliation_id = Column(Integer, ForeignKey('affiliations.id'), primary_key=True)


# RSID's are only stored if linked to pubmed articles
class RSIDCitation(Base):
    __tablename__ = 'rsid_citations'
    rsid = Column(Integer, primary_key=True)
    pmid = Column(Integer, ForeignKey('articles.pmid'), primary_key=True)

    # Generate list of tuples for eventual vcf output
    def vcf_data(self):
        article = session.query(Article).get(self.pmid)
        return [(self.rsid, self.pmid, noun) for noun in article.unique_abstract_nouns()]


def build_complete_noun_mapping():
    mapping = {}
    articles = session.query(Article).all()
    count = 0
    with progressbar.ProgressBar(max_value=len(articles)) as bar:
        for article in articles:
            count += 1
            mapping[str(article.pmid)] = sorted(article.abstract_nouns())
            bar.update(count)
    with open(MAPPING_FILE, 'w') as json_file:
        json.dump(mapping, json_file, sort_keys=True, indent=4)


def save_noun_mapping(mapping):
    with open(MAPPING_FILE, 'w') as json_file:
        json.dump(mapping, json_file, sort_keys=True, indent=4)


def load_noun_mapping():
    with open(MAPPING_FILE, 'r') as json_file:
        return json.load(json_file)


def add_articles_to_mapping(articles):
    mapping = load_noun_mapping()
    count = 0
    with progressbar.ProgressBar(max_value=len(articles)) as bar:
        for article in articles:
            count += 1
            mapping[str(article.pmid)] = sorted(article.abstract_nouns())
            bar.update(count)
    save_noun_mapping(mapping)


article_noun_mapping = load_noun_mapping()


# Create the database schema defined above
def initialize_database():
    Base.metadata.create_all(engine)


# From model and attributes, create model if it doesn't exist and return the model
def model_ensure_present(model, **kwargs):
    model_query = session.query(model).filter_by(**kwargs)
    if not model_query.count():
        session.add(
            model(**kwargs)
        )
        session.commit()
    return model_query.first()


# Add author models and affiliations from ncbiutils.Author object
def load_from_author(author):
    try:
        author_model = model_ensure_present(
            Author,
            last_name=author.last_name,
            first_name=author.first_name,
            initials=author.initials
        )
    except AttributeError:
        author_model = model_ensure_present(
            Author,
            collective_name=author.collective_name
        )
    for affiliation in author.affiliations:
        affiliation_model = model_ensure_present(
            Affiliation,
            info=affiliation
        )
        model_ensure_present(
            AuthorAffiliation,
            author_id=author_model.id,
            affiliation_id=affiliation_model.id
        )
    return author_model


# Load from pubmed article object
def load_from_pubmed_article(article):
    article_model = model_ensure_present(
        Article,
        pmid=article.pmid
    )
    article_model.abstract = article.abstract
    article_model.date_created = article.date_created
    article_model.date_revised = article.date_revised
    session.commit()
    # Cite rsids
    for rsid in article.rsids():
        model_ensure_present(
            RSIDCitation,
            rsid=rsid,
            pmid=article.pmid
        )
    # Create/cite authors
    for author in article.authors:
        author_model = load_from_author(author)
        model_ensure_present(
            AuthorCitation,
            pmid=article.pmid,
            author_id=author_model.id
        )
    return article_model


# Get all cited pubmed articles, load associated data
def load_all_data():
    pubmed_articles = get_pubmed_articles()
    count = 0
    print("Loading articles into database...")
    with progressbar.ProgressBar(max_value=len(pubmed_articles)) as bar:
        for article in pubmed_articles:
            count += 1
            load_from_pubmed_article(article)
            bar.update(count)


# If load_all_data breaks partway through, only loads data for articles that aren't in the database
def load_data_not_present():
    pubmed_articles = get_pubmed_articles()
    count=0
    with progressbar.ProgressBar(max_value=len(pubmed_articles)) as bar:
        for article in pubmed_articles:
            count += 1
            if not session.query(Article).filter_by(pmid=article.pmid).count():
                load_from_pubmed_article(article)
            bar.update(count)


# Get any articles revised since the latest date_revised in database and load them
def update_data():
    latest_change = session.query(func.max(Article.date_revised)).first()[0]
    updated_pubmed_articles = get_pubmed_articles(since_date=latest_change)
    article_models = []
    count = 0
    with progressbar.ProgressBar(max_value=len(updated_pubmed_articles)) as bar:
        for article in updated_pubmed_articles:
            count += 1
            article_model = load_from_pubmed_article(article)
            article_models.append(article_model)
            bar.update(count)
    print("Updating mapping file...")
    add_articles_to_mapping(article_models)

