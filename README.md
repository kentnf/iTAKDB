
Pipeline For update iTAK Databse
--------------------------------

### A. Parse protein and CDS sequence ID:
 
     $parse_protein_id.pl

		SeqID before parse: >g18373.t1|PACid:27562759
		SeqID after parse : >g18373.t1

     *Important*
     After that, please add Rice and Arabidopsis output to your folder for making clusters. 
     For some species, using PACid is better than using original transcript ID

### B. Create a list for plant peptide and CDS sequences. The format of the list is like below:

     Each line include five columns, they are:
		A: folder for species (rice)
		B: monocotyledon / dicotyledon / non-angiosperms
		C: Name of plant species

     The protein, CDS and other files are named according to A: folder for species
		A: CDS sequence				-- rice_cds
		B: protein				-- rice_pep
		C: representive protein			-- rice_rep_pep (optional)
		D: protein & transcirpt & gene file	-- rice_trans_gene
		E: gene position file			-- rice_gene_position (optional)
		F: chr Size file			-- rice_chrSize (optional)
		G: block file				-- rice_block (optional)
		H: Other block file			-- rice_other_block (optional)

     Example:
		A: publishedPlantGenome/rice  
		B: monocotyledon     
		C: Rice   

		Below info do not have order
		A: publishedPlantGenome/rice/rice_cds
		B: publishedPlantGenome/rice/rice_pep
		C: publishedPlantGenome/rice/rice_rep_pep
		D: publishedPlantGenome/rice/rice_trans_gene
		E: publishedPlantGenome/rice/rice_gene_position	
		F: publishedPlantGenome/rice/rice_chrSize
		G: publishedPlantGenome/rice/rice_block
		H: publishedPlantGenome/rice/rice_other_block
		* other in H should be name of other species

     Each line of protein & transcirpt & gene file include three columns:
		A: protein ID
		B: transcript ID
		C: gene ID

     Each line of gene position file:
		A: chrID
		B: start
		C: end
		E: strand
		F: geneID/mRNA

     Each line of chr Size file
		A: chrID
		B: chrName
		C: length
		D: color

### C. run iTAK.

     Command or make a bash for each specie:
     		$perl /path/iTAK.pl -i apple_pep

     Input:
     		The peptide sequence

     Output:
     		Each plant species will generates one output folder in publishedPlantGenome folder, and each folder will have 6 files:
		apple_pep_pkaln : alignment of protein kinase  
		apple_pep_pkcat : assign each gene to different PK families
		apple_pep_pkseq : protein kinase sequences
		apple_pep_tf_align : alignment of transcription factor
		apple_pep_tf_family : assign each gene to TF families
		apple_pep_tf_seq : TF sequences

### D. generate file for updating iTAK database

     $perl itak_DB_update.pl  list

     Input:
     The list created in step 1. 

     Output:
     The output includes files that categorized to five parts (folders) for iTAK databases.

     1. for_database -- include files that need to imported to mysql database
          
		category
		family_summary  
		protein_domain_table  
		protein_family_table  
		gene_table
		*protein_annotation (generated at step 5)

     2. for_download -- include files that could be download from the iTAK website

		All_PKs_list
		All_PKs_protein
		All_TFs_TRs_nucleotide
		All_PKs_nucleotide
		All_TFs_TRs_list
		All_TFs_TRs_protein

		protein_kinases [folder]
		transcription_factors_and_transcriptional_regulators [folder]

     3. for_blast -- include files that need to be blast online, and gene annotation
		
		All_PKs_nucleotide
		All_PKs_protein
		All_TFs_TRs_nucleotide
		All_TFs_TRs_protein

     4. for_cluster / for_tree -- include nwk files generated by clustalW2

     5. for_synteny -- include gene list for synteny analysis
     
### E. Blast protein sequence against Genebank for annotation

     1. prepare the protein sequence for blast
     		$perl iTAK_DB_prepare_blast.pl  pre_TFs_proteins  pre_PKs_proteins  new_TFs_proteins  new_PKs_proteins
     		* the output file is Prepared_Protein.fa

     2. upload protein sequence to codequest Server

     3. parse the blast result:
     		$perl  iTAK_DB_parse_blast.pl  blast_results  > update_protein_annotation_version
     		* the output file is protein_annotation_table file for importing to database

### F. Upload files to server

     1. for_database -- /var/www/cgi-bin/itak/mysql_database/versionN/
     2. for_download -- /var/ftp/pub/program/itak/database
     3. for_blast    -- /var/www/html/tool/itak/blast_database
     4. for_cluster  -- /var/www/html/tool/itak/phylogenetic_tree
        for_tree     -- /var/www/html/tool/itak/phylo_tree

### G. Load database to iTAK SQL database 

     command for import these files to database: 
     LOAD DATA LOCAL INFILE 'family_summary' INTO TABLE `family` FIELDS TERMINATED BY '\t' ESCAPED BY '\\' LINES TERMINATED BY '\n';

     Here is the command for importing each table:
		mysql> truncate table itak.category;
		mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/category" into table itak.categroy;

		mysql> truncate table itak.gene;
		mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/gene_table" into table itak.gene;

		mysql> truncate table itak.family_summary;
		mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/family_summary" into table itak.family_summary;

		mysql> truncate table itak.gene_domain;
		mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/protein_domain_table" into table itak.gene_domain;

		mysql> truncate table itak.gene_family;
		mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/protein_family_table" into table itak.gene_family;

		** mysql> truncate table itak.gene_annotation; (do not need truncate this table)
		mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/update_protein_annotation_version" into table itak.gene_annotation;

     Problem: some genes assigned to both PKs and TFs families: please check gene_family table

### H. Format the database

     Using formatdb command to format the blast files
     $formatdb -i All_PKs_nucleotide -p F
     $formatdb -i All_PKs_protein -p T
     $formatdb -i All_TFs_TRs_nucleotide -p F
     $formatdb -i All_TFs_TRs_protein -p T

### G. Update webpage

     1. db_home.cgi
		add plant info to data hash
		add plant name to tree base on http://genomevolution.org/wiki/index.php/Sequenced_plant_genomes

Ref: Generate and Draw phylogenetic tree automatically   
------------------------------------------------------
      A. Generate phylogenetic tree using MEGA

      		1. using MEGA-CC can generate phylogenetic tree automatically 
      		(http://update.megasoftware.net/download.php?email=yz357@cornell.edu&location=primary&version=win-5.1CC)
      		(http://www.megasoftware.net/download.php?email=yz357@cornell.edu&location=secondary&version=win-5.1CC)

      B. Draw phylogenetic tree
		1. convert newick format to phyloXML format using forester_1019.jar (http://code.google.com/p/forester/)
		2. display phyloXML using jsPhyloSVG (http://www.jsphylosvg.com/)

      *using perl to construct the pipeline
      *about phyloxml (http://www.phyloxml.org/), including many tools to display the trees.

Ref: Parameters for Blast
-------------------------
      For your iTAK protein blast, can you tell me what parameters you want to use, specifically 
      "Max Scores"
      "Max Alignments"
      "Significance Threshold"
      http://132.236.156.152/decypher/algo-tera-blast/tera-blastp_aa.shtml

Ref: How to backup/install iTAK database on a different computer
----------------------------------------------------------------
     Before backup/install, your destination computer already has apache and mysql installed.

     Step1. Copy /var/www/html/tools/itak/ and /var/www/cgi-bin/itak/ folders from server to your destination computer.

     Step2. Create itak database and import database file into it. 
	$ mysql -u root -p                                          # enter mysql
	mysql> create database itak;                          # create itak database
	$ mysql -uroot -p itak< itak_db_stru.mysql;      # import the structure of itak tables

	# load dataset to itak tables
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/category" into table itak.categroy;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/domain_table" into table itak.domain;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/family_table" into table itak.family;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/gene_table" into table itak.gene;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/family_summary" into table itak.family_summary;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/gene_domain_table" into table itak.gene_domain;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/gene_family_table" into table itak.gene_family;
	mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/versionN/gene_annotation" into table itak.gene_annotation;

# create user for itak database
	mysql> create user itak;
	mysql> grant select [insert, update, delete] on itak.* to itak@localhost identified by "itak"

     Step3. Change the limit of folders and files using chmod command.

     Step4. Edit the itak.pm base on the environment of destination computer.
	our $tmp = "/var/www/html/tmp";
	our $home = "/var/www/html/tool/itak";
	our $ftp_url = "ftp://bioinfo.bti.cornell.edu/pub/program/itak";
	our $http_url = "http://bioinfo.bti.cornell.edu/tool/itak";
	our $root_url = "http://bioinfo.bti.cornell.edu";

#### * if you have any problem on mysql db, please check the mysql db card in Evernote

	Ref: The structure of iTAK database
	mysql> use database;
	mysql> show tables;
	+-----------------+
	| Tables_in_itak  |
	+-----------------+
	| domain          |
	| family          |
	| family_summary  |
	| gene            |
	| gene_annotation |
	| gene_domain     |
	| gene_family     |
	+-----------------+

	use the mysql command to show structure of tables:

	mysql> desc category;
	+----------+-------------+------+-----+---------+-------+
	| Field    | Type        | Null | Key | Default | Extra |
	+----------+-------------+------+-----+---------+-------+
	| organism | varchar(50) | NO   | MUL | NULL    |       |
	| category | varchar(50) | NO   |     | NULL    |       |
	+----------+-------------+------+-----+---------+-------+

	mysql> desc domain;
	+-------------+-------------+------+-----+---------+-------+
	| Field       | Type        | Null | Key | Default | Extra |
	+-------------+-------------+------+-----+---------+-------+
	| domain_ID   | varchar(15) | NO   | PRI | NULL    |       |
	| domain_desc | text        | YES  |     | NULL    |       |
	| category    | varchar(15) | YES  |     | NULL    |       |
	+-------------+-------------+------+-----+---------+-------+
	This table do not need update!!!

	mysql> desc family;
	+-------------+--------------+------+-----+---------+-------+
	| Field       | Type         | Null | Key | Default | Extra |
	+-------------+--------------+------+-----+---------+-------+
	| family_acc  | varchar(25)  | NO   | PRI | NULL    |       |
	| description | varchar(255) | YES  |     | NULL    |       |
	| category    | varchar(50)  | YES  |     | NULL    |       |
	+-------------+--------------+------+-----+---------+-------+
	This table do not need update!!!

	mysql> desc family_summary;
	+-------------+-------------+------+-----+---------+-------+
	| Field       | Type        | Null | Key | Default | Extra |
	+-------------+-------------+------+-----+---------+-------+
	| family_acc  | varchar(25) | NO   | PRI | NULL    |       |
	| Number_gene | smallint(6) | YES  |     | NULL    |       |
	| organism    | varchar(50) | NO   | PRI | NULL    |       |
	+-------------+-------------+------+-----+---------+-------+

	mysql> desc gene;
	+-------------+-------------+------+-----+---------+-------+
	| Field       | Type        | Null | Key | Default | Extra |
	+-------------+-------------+------+-----+---------+-------+
	| gene_ID     | varchar(50) | NO   | PRI | NULL    |       |
	| protein_seq | text        | YES  |     | NULL    |       |
	| CDS_seq     | text        | YES  |     | NULL    |       |
	| organism    | varchar(50) | NO   | PRI |         |       |
	+-------------+-------------+------+-----+---------+-------+

	mysql> desc gene_annotation;
	+-------------+-------------+------+-----+---------+-------+
	| Field       | Type        | Null | Key | Default | Extra |
	+-------------+-------------+------+-----+---------+-------+
	| trans_ID    | varchar(50) | NO   | PRI | NULL    |       |
	| gene_ID     | varchar(50) | NO   |     | NULL    |       |
	| protein_seq | text        | YES  |     | NULL    |       |
	| CDS_seq     | text        | YES  |     | NULL    |       |
	| organism    | varchar(50) | NO   | PRI |         |       |
	+-------------+-------------+------+-----+---------+-------+

	+----------+-------------+------+-----+---------+-------+
	| Field    | Type        | Null | Key | Default | Extra |
	+----------+-------------+------+-----+---------+-------+
	| gene_ID  | varchar(50) | NO   | PRI | NULL    |       |
	| hit_acc  | varchar(50) | NO   | PRI | NULL    |       |
	| hit_desc | text        | NO   |     | NULL    |       |
	| score    | float       | NO   |     | NULL    |       |
	| evalue   | varchar(10) | NO   |     | NULL    |       |
	+----------+-------------+------+-----+---------+-------+

	mysql> desc gene_domain;
	+--------------+-------------+------+-----+---------+-------+
	| Field        | Type        | Null | Key | Default | Extra |
	+--------------+-------------+------+-----+---------+-------+
	| gene_ID      | varchar(50) | NO   | PRI |         |       |
	| domain_ID    | varchar(50) | NO   | PRI |         |       |
	| query_start  | smallint(6) | NO   | PRI | 0       |       |
	| query_end    | smallint(6) | NO   | PRI | 0       |       |
	| target_start | smallint(6) | YES  |     | NULL    |       |
	| target_end   | smallint(6) | YES  |     | NULL    |       |
	| align_query  | text        | YES  |     | NULL    |       |
	| align_string | text        | YES  |     | NULL    |       |
	| align_hit    | text        | YES  |     | NULL    |       |
	| score        | float       | YES  |     | NULL    |       |
	| evalue       | varchar(10) | YES  |     | NULL    |       |
	+--------------+-------------+------+-----+---------+-------+

	mysql> desc gene_family;
	+------------+-------------+------+-----+---------+-------+
	| Field      | Type        | Null | Key | Default | Extra |
	+------------+-------------+------+-----+---------+-------+
	| gene_ID    | varchar(50) | NO   | PRI | NULL    |       |
	| family_acc | varchar(25) | NO   | PRI | NULL    |       |
	+------------+-------------+------+-----+---------+-------+

