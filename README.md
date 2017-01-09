# PhenVar
PhenVar is designed to take one or more rsids and generate a list of PubMed IDs to query and generate novel associations between publications. 
## Known Issues / Improvements
#### True issues and improvements to be moved to GH issue tracking in the future
Currently the articles are being downloaded one by one and the abstracts appended to a list (prior to any language processing).  They could be downloaded in bulk, but ultimately this will use more memory.  If the datasets got sizeable enough where this should be considered an option, it should be added as a configuration option. Or proper block sizing that could change to increase speed / efficiency vs. not.  
