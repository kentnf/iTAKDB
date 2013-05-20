
Pipeline For update iTAK Databse
--------------------------------

1. Parse protein and CDS sequence ID:
 
    $parse_protein_id.pl

    SeqID before parse: >g18373.t1|PACid:27562759
    SeqID after parse : >g18373.t1

    *Important*
    After that, please add Rice and Arabidopsis output to your folder for making clusters. 

2. Create a list for Plant peptide and cds sequences. The format of the list is like below:

Each line include four columns, they are:
     A: Path of peptide sequence
     B: monocotyledon dicotyledon
     C: Name of plant species
     D: Path of CDS sequence

Example:     
publishedPlantGenome/apple_pep  dicotyledon     Apple   publishedPlantGenome/apple_cds

2. run iTAK.

Command or make a bash for each specie: 
perl itak_db_step1.pl  list

perl ../../iTAK-1.2.x64/iTAK.pl -i grape_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i maize_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i cucumber_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i date_palm_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i medicago_truncatula_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i papaya_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i pigeonpea_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i poplar_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i potato_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i rice_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i selmo1_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i sorghum_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i soybean_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i strawberry_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i TAIR10_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i thellugiella_parvula_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i tomato_pep -a 2
perl ../../iTAK-1.2.x64/iTAK.pl -i watermelon_pep -a 2

Input:
The list created in step 1. 

Output:
Each line (species) will generates one output folder in publishedPlantGenome folder. 
And each output will have 6 files:
apple_pep_pkaln : alignment of protein kinase  
apple_pep_pkcat : assign each gene to different PK families
apple_pep_pkseq : protein kinase sequences
apple_pep_tf_align : alignment of transcription factor
apple_pep_tf_family : assign each gene to TF families
apple_pep_tf_seq : TF sequences

3. prepare files for iTAK database

Command:
perl itak_db_step2.pl  list

******************************
The itak_db_step2.pl needs to be edit for some point. 
Point A: convert CDS ID to protein ID : For some genomes, the CDS id is little difference with Protein ID, so we must change the ID to get CDS sequences. The point A has two place for both PKs and TFs.



*******************************
Input:
The list created in step 1. 

Output:
The output includes files that categorized to five parts (folders) for iTAK databases.

     1. for_database  -- include files that need to imported to mysql database
          
          family_summary  
          gene_annotation_table  
          gene_domain_table  
          gene_family_table  
          gene_table

          Delete the old records using itak_db_delete.pl:
          A: prepare a file including the organism name that need to be delete.
          B: delete the records using itak_db_delete.pl
          perl  itak_db_delete.pl  list_of_organism 
          * this script will generate list_of_genes that need to be delete.

*Important*
Please delete Rice and Arabidopsis results from the Database First or together with update data. 

          command for import these files to database: 
     LOAD DATA LOCAL INFILE 'family_summary' INTO TABLE `family` FIELDS TERMINATED BY '\t' ESCAPED BY '\\' LINES TERMINATED BY '\n';

          Here is the command for importing each table: 
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_table" into table itak.gene;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/family_summary" into table itak.family_summary;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_domain_table" into table itak.gene_domain;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_family_table" into table itak.gene_family;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_annotation" into table itak.gene_annotation;


mysql> load data local infile "/root/to_db" into table itak.gene_annotation;
Query OK, 713599 rows affected, 19 warnings (23.08 sec)
Records: 713833  Deleted: 0  Skipped: 234  Warnings: 19

* 19..... score
* 234............transID, hitACC..........evalue..score
* ......
          
     2. for_blast  -- include files that need to be blast against with GenBank database
          
          All_PKs_nucleotide
          All_PKs_protein
          All_TFs_TRs_nucleotide
          All_TFs_TRs_protein

          blast protein sequence against Genebank using blastp, then parse the blast result:
          Perl  itak_db_parse_blast.pl  blast_results
          * the output file is gene_annotation file for importing to database 

          Here is the command for importing gene_annotation table: 
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_annotation" into table itak.gene_annotation;

          Problem: some genes assigned to both PKs and TFs families: 
          /home/kentnf/project/itak/itak_database/for_blast/combine_family.

          Combine files with old one.
          The new files just have new and updated sequences, it should be combined with old ones.
          Using command:  (how about re-generated)

          Using formatdb command to format the blast files
          formatdb -i All_PKs_nucleotide -p F
          formatdb -i All_PKs_protein -p T
          formatdb -i All_TFs_TRs_nucleotide -p F
          formatdb -i All_TFs_TRs_protein -p T

     3. for_download -- include files that could be download from the iTAK website

          All_PKs_list
          All_PKs_protein
          All_TFs_TRs_nucleotide
          All_PKs_nucleotide
          All_TFs_TRs_list
          All_TFs_TRs_protein

          protein_kinases [folder]
          transcription_factors_and_transcriptional_regulators [folder]

     4. for_cluster -- include files that prepared for clustering (Plants and Arabidopsis/Rice)

     5. for_tree -- include files that prepared for clustering  (Just Plants)

4. Clustering

     1. list all the peptide sequence files for each family.
     $ list   for_cluster/*/*/*.pep   >   list_cluster
     $ list   for_tree/*/*/*.pep   >  list_tree

     Format of list files (list_cluster & list_tree):
     for_cluster/dicotyledon/tf/Apple_ABI3VP1.pep      nocluster
     for_cluster/dicotyledon/tf/Apple_AP2-EREBP.pep  nocluster
     ......                                                                 ......

     2. clustering files in cluster and tree folders
     $ perl  itak_db_step3.pl  list_cluster
     $ perl  itak_db_step3.pl  list_tree

     Input
     list of clusters and trees for each peptide sequence file (step 4.2)

     Output
     *generate five type of files (.aln .dnd .nwk .pep .phb)
     *reproduce the list files including more info about processing of cluster

     Format of the new list: 
     for_tree/dicotyledon/tf/Apple_ABI3VP1.pep     80     clustered
     
5. Draw images for clustering at step 4

     Input
     list of clusters and trees for each peptide sequence file (output of clustering at step 4)

     Output
     for_mega_cluster : folder includes clustered files of list_cluster
     for_mega_tree     : folder includes clustered files of list_tree

& generate images using MEGA software, upload images, clustering Trees to iTAK web &



Plan 1 ..........................................one month.

..............................................................................
.....................ModelTest.............Model....ModelTest.2005....ModelTest3.7.........DOS.................................MrMTgui......ModelTest...........ModelTest.......PAUP..........Model........ModelTest...............Model......
1......ClustalX 1.83..........PAUP.....Nexus......PAUP.
2.........Modeltest.........paupblock............. 
modelblockPAUPb10....paup....file.open..............paup.ML
....model...........................
3..........paupblock............model.scores....
4...MrMTgui................MT path.PAUP path............
Modeltest.paup.............................
5..........select file...model.scores........
Modeltest...............
6..............hLRT..........AIC......
Likelihood settings
from best-fit model (TVM+G) selected by AIC in Modeltest 3.7 on Tue May 26 09:57:45 
2009]
BEGIN PAUP;
Lset Base=(0.4370 0.1035 0.1546) Nst=6 Rmat=(0.2914 3.9607 2.9280 0.8550 3.9607) 
Rates=gamma Shape=1.5528 Pinvar=0;
Begin paup.........paup...........................


.jModeltest .......

.............cvtree

Ref: Generate and Draw phylogenetic tree automatically   
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
For your iTAK protein blast, can you tell me what parameters you want to use, specifically "Max Scores", "Max Alignments" and "Significance Threshold"? http://132.236.156.152/decypher/algo-tera-blast/tera-blastp_aa.shtml

Ref: How to backup/install iTAK database on a different computer
Before backup/install, your destination computer already has apache and mysql installed.

Step1. Copy /var/www/html/tools/itak/ and /var/www/cgi-bin/itak/ folders from server to your destination computer.

Step2. Create itak database and import database file into it. 
$ mysql -u root -p                                          # enter mysql
mysql> create database itak;                          # create itak database
$ mysql -uroot -p itak< itak_db_stru.mysql;      # import the structure of itak tables

# load dataset to itak tables
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/domain_table" into table itak.domain;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/family_table" into table itak.family;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_table" into table itak.gene;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/family_summary" into table itak.family_summary;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_domain_table" into table itak.gene_domain;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_family_table" into table itak.gene_family;
mysql> load data local infile "/var/www/cgi-bin/itak/mysql_database/version2/gene_annotation" into table itak.gene_annotation;

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

* if you have any problem on mysql db, please check the mysql db card in Evernote

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

