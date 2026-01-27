-- Preprocess ZFIN ortholog files into a single combined TSV
-- This SQL script is executed by DuckDB via preprocess.py

create or replace table zfin_orthologs as
(
    select 'ZFIN:' || "ZFIN ID" as zfin_gene, 'FB:' || "Flybase ID" as ortholog_gene, evidence as evidence, string_agg("Pub ID",'|') as publications
    from read_csv('data/fly_orthos.txt', names=['ZFIN ID', 'ZFIN Symbol', 'ZFIN Name', 'Fly Symbol', 'Fly Name"', 'Flybase ID', 'Gene ID', 'Evidence', 'Pub ID'])
    group by all
    union
    select 'ZFIN:' || "ZFIN ID" as zfin_gene, "MGI ID" as ortholog_gene, evidence as evidence, string_agg("Pub ID",'|') as publications
    from read_csv('data/mouse_orthos.txt', names=['ZFIN ID', 'ZFIN Symbol', 'ZFIN Name', 'Mouse Symbol', 'Mouse Name"', 'MGI ID', 'Gene ID', 'Evidence', 'Pub ID'])
    group by all
    union
    select 'ZFIN:' || "ZFIN ID" as zfin_gene, 'HGNC:' || "HGNC ID" as ortholog_gene, evidence as evidence, string_agg("Pub ID",'|') as publications
    from read_csv('data/human_orthos.txt', names=['ZFIN ID', 'ZFIN Symbol', 'ZFIN Name', 'Human Symbol', 'Human Name"', 'OMIM ID', 'Gene ID', 'HGNC ID', 'Evidence', 'Pub ID'])
    where "HGNC ID" is not null
    group by all
);

copy (select * from zfin_orthologs) to 'data/zfin_orthologs.tsv' (HEADER, DELIMITER '\t');
