CREATE TABLE asarLite.taxLite (
  ID      serial,
  taxID   INTEGER,
  species varchar(2055),
  genus   varchar(2055),
  family  varchar(2055),
  class   varchar(2055),
  phylum  varchar(2055),
  dom     varchar(2055),
  name    varchar(2055) NOT NULL,
  usp     varchar(2055) NOT NULL,
  ordr    varchar(2055));
CREATE TABLE asarLite.md5 (
  md5sum char(32) NOT NULL,
  PRIMARY KEY (md5sum));
CREATE TABLE asarLite.KOrt (
  id         INTEGER NOT NULL,
  access     varchar(25) NOT NULL,
  funct varchar(4055) NOT NULL,
  description varchar(8255),
  PRIMARY KEY (id));
CREATE TABLE asarLite.keggpath (
  ID          INTEGER NOT NULL,
  name        varchar(4055) NOT NULL,
  description varchar(8255),
  PRIMARY KEY (ID));
CREATE TABLE asarLite.study (
  ID          serial,
  name        varchar(255) NOT NULL,
  mgrastid varchar(55),
  description varchar(8255),
  metaID      INTEGER);
CREATE TABLE asarLite.metaproject (
  ID          serial,
  name        varchar(255) NOT NULL,
  description varchar(8255));
CREATE TABLE asarLite.smpl (
  ID          serial,
  name        varchar(2055) NOT NULL,
  studyID     INTEGER NOT NULL,
  mgrastid varchar(255),
  description varchar(8255),
  eval        INTEGER DEFAULT 5,
  ident  INTEGER DEFAULT 60,
  length      INTEGER DEFAULT 15,
  orgSource   varchar(255),
  funSource   varchar(255));
CREATE TABLE asarLite.properties (
  id          SERIAL,
  studyID     INTEGER NOT NULL,
  name        varchar(255) NOT NULL,
  domain      varchar(4255),
  description varchar(8255));
CREATE TABLE asarLite.collection (
  ID          serial,
  name        varchar(255) NOT NULL UNIQUE,
  description varchar(8255));
CREATE TABLE asarLite.keggpath_KOrt (
  keggpathID INTEGER NOT NULL,
  KOrtid     INTEGER NOT NULL,
  PRIMARY KEY (keggpathID,
  KOrtid));
CREATE TABLE asarLite.propvalue (
  sampleID     INTEGER NOT NULL,
  propertiesid INTEGER NOT NULL,
  val          varchar(4255) NOT NULL,
  PRIMARY KEY (sampleID,
  propertiesid));
CREATE TABLE asarLite.collection_sample (
  collectionID INTEGER NOT NULL,
  sampleID     INTEGER NOT NULL,
  PRIMARY KEY (collectionID,
  sampleID));
CREATE TABLE asarLite.md5_taxLite (
  md5md5sum char(32) NOT NULL,
  taxLiteID INTEGER NOT NULL,
  PRIMARY KEY (md5md5sum,
  taxLiteID));
CREATE TABLE asarLite.md5_KOrt (
  md5md5sum char(32) NOT NULL,
  KOrtid    INTEGER NOT NULL,
  PRIMARY KEY (md5md5sum,
  KOrtid));
CREATE TABLE asarLite.abundance (
  md5md5sum char(32) NOT NULL,
  sampleID  INTEGER NOT NULL,
  abundance BIGINT NOT NULL,
  PRIMARY KEY (md5md5sum,
  sampleID));
CREATE TABLE asarLite.seedLite (
  ID   serial,
  accession varchar(32) NOT NULL,
  FUN1 varchar(4055) NOT NULL,
  FUN2 varchar(4055) NOT NULL,
  FUN3 varchar(4055) NOT NULL,
  FUN4 varchar(4055) NOT NULL,
  CONSTRAINT funcUnq
    UNIQUE (FUN1, FUN2, FUN3, FUN4));
CREATE TABLE asarLite.seedLite_md5 (
  seedLiteID INTEGER NOT NULL,
  md5md5sum  char(32) NOT NULL,
  PRIMARY KEY (seedLiteID,
  md5md5sum));
CREATE TABLE asarLite.tmp_md5 (
  md5md5sum  char(32) NOT NULL,
  cnt INTEGER NOT NULL);
CREATE TABLE asarLite.tmp_annot (
  mgrastid  varchar(255) NOT NULL,
  md5sum        char(32) NOT NULL,
  ab        BIGINT,
  accession     varchar(2055) NOT NULL);
CREATE TABLE asarLite.tmp_kort (
  ID INTEGER NOT NULL,
  access     varchar(25) NOT NULL,
  funct varchar(4055) NOT NULL,
  description varchar(8255),
  keggpathID INTEGER);
CREATE TABLE asarLite.tmp_seed (
  mgrastid  varchar(255),
  md5sum    char(32),
  ab        BIGINT,
  uFUN varchar(4055),
  usp     varchar(2055)
  );
CREATE TABLE asarLite.tmp_prepMD5 (
  md5    char(32),
  src varchar(255),
  fun     varchar(2055),
  accession     varchar(2055),
  ncbi_tax_id     varchar(2055),
  organism     varchar(2055),
  types varchar(255)
  );
CREATE TABLE asarLite.mv_smpl_stat (
  stname        varchar(255) NOT NULL,
  project varchar(55) NOT NULL,
  ID          INTEGER NOT NULL,
  name        varchar(2055) NOT NULL,
  smpl varchar(255) NOT NULL,
  cnt INTEGER NOT NULL,
  abundance BIGINT NOT NULL);
create view smpl_stat as
  select t.name as stname,
         t.mgrastid as project,
         s.id,
         s.name,
         s.mgrastid as smpl,
         count(distinct m.md5md5sum) as cnt,
         sum(abundance) as abundance
  from study t join
       smpl s on t.id=s.studyid left join
       abundance m on s.id=m.sampleid
  group by stname,project,s.id,s.name,smpl;
create view funtaxall as
select distinct a.md5md5sum as md5sum,
       s.mgrastid,
       a.abundance,
       usp,species,genus,family,ordr,class,phylum,dom,
       fun4 as ufun, fun3 as fun2, fun2 as fun3, fun4 as fun1
  from abundance a join smpl s on a.sampleid=s.id join
       seedlite_md5 sm on a.md5md5sum=sm.md5md5sum join
       md5_taxlite tm on a.md5md5sum=tm.md5md5sum  join
       seedlite sl on sl.id=sm.seedliteid join
       taxlite tl on tl.id=tm.taxliteid;
CREATE TABLE "asarlite"."mv_abund_stat" (
	"sampleid"  INTEGER       NOT NULL,
	"cnt"       INTEGER       NOT NULL,
	"abundance" BIGINT        NOT NULL
);
CREATE TABLE asarLite.mv_funtaxab (
  md5sum    char(32) NOT NULL,
  mgrastid  varchar(255) NOT NULL,
  abundance BIGINT NOT NULL,
  usp       varchar(2055) NOT NULL,
  ufun      varchar(32) NOT NULL);
create view funtaxab as
select distinct a.md5md5sum as md5sum,
 s.mgrastid,
 a.abundance,
 usp,
 sl.accession as ufun
 from abundance a join smpl s on a.sampleid=s.id join
 seedlite_md5 sm on a.md5md5sum=sm.md5md5sum join
 md5_taxlite tm on a.md5md5sum=tm.md5md5sum join
 seedlite sl on sl.id=sm.seedliteid join
 taxlite tl on tl.id=tm.taxliteid;
CREATE TABLE asarlite.mv_dkres (
       id         VARCHAR(255)  NOT NULL,
       md5        CHAR(32)      NOT NULL,
       annotation VARCHAR(8255) NOT NULL,
       ko         VARCHAR(25)   NOT NULL
);
create view dkres as
select distinct s.mgrastid as id,m.md5md5sum as md5,
       kort.description as annotation,
       kort.access as ko
   from smpl s join
        abundance a on s.id=a.sampleid join
        md5_kort m on m.md5md5sum=a.md5md5sum join
        kort on kort.id=m.kortid;
CREATE INDEX taxLite_species
  ON asarLite.taxLite (species);
CREATE INDEX taxLite_genus
  ON asarLite.taxLite (genus);
CREATE INDEX taxLite_family
  ON asarLite.taxLite (family);
CREATE INDEX taxLite_class
  ON asarLite.taxLite (class);
CREATE INDEX taxLite_phylum
  ON asarLite.taxLite (phylum);
CREATE INDEX taxLite_domain
  ON asarLite.taxLite (dom);
CREATE INDEX taxLite_name
  ON asarLite.taxLite (name);
CREATE INDEX taxLite_taxid
  ON asarLite.taxLite (taxid);
CREATE UNIQUE INDEX KOrt_access
  ON asarLite.KOrt (access);
CREATE INDEX KOrt_funct
  ON asarLite.KOrt (funct);
CREATE INDEX keggpath_name
  ON asarLite.keggpath (name);
CREATE UNIQUE INDEX study_name
  ON asarLite.study (name);
CREATE INDEX id
  ON asarLite.study (mgrastid);
CREATE UNIQUE INDEX metaproject_name
  ON asarLite.metaproject (name);
CREATE UNIQUE INDEX sample_name
  ON asarLite.smpl (name);
CREATE UNIQUE INDEX properties
  ON asarLite.properties (studyID, name);
CREATE INDEX seedLite_acc
  ON asarLite.seedLite (accession);
CREATE INDEX seedLite_FUN1
  ON asarLite.seedLite (FUN1);
CREATE INDEX seedLite_FUN2
  ON asarLite.seedLite (FUN2);
CREATE INDEX seedLite_FUN3
  ON asarLite.seedLite (FUN3);
CREATE INDEX seedLite_FUN4
  ON asarLite.seedLite (FUN4);
ALTER TABLE asarLite.study ADD CONSTRAINT FKstudy265557 FOREIGN KEY (metaID) REFERENCES asarLite.metaproject (ID);
ALTER TABLE asarLite.smpl ADD CONSTRAINT FKsample944577 FOREIGN KEY (studyID) REFERENCES asarLite.study (ID);
ALTER TABLE asarLite.properties ADD CONSTRAINT FKproperties566586 FOREIGN KEY (studyID) REFERENCES asarLite.study (ID);
ALTER TABLE asarLite.keggpath_KOrt ADD CONSTRAINT FKkeggpath_K508193 FOREIGN KEY (keggpathID) REFERENCES asarLite.keggpath (ID);
ALTER TABLE asarLite.keggpath_KOrt ADD CONSTRAINT FKkeggpath_K791205 FOREIGN KEY (KOrtid) REFERENCES asarLite.KOrt (id);
ALTER TABLE asarLite.propvalue ADD CONSTRAINT FKvalue866516 FOREIGN KEY (sampleID) REFERENCES asarLite.smpl (ID);
ALTER TABLE asarLite.propvalue ADD CONSTRAINT FKvalue844367 FOREIGN KEY (propertiesid) REFERENCES asarLite.properties (id);
ALTER TABLE asarLite.collection_sample ADD CONSTRAINT FKcollection570840 FOREIGN KEY (collectionID) REFERENCES asarLite.collection (ID);
ALTER TABLE asarLite.collection_sample ADD CONSTRAINT FKcollection272515 FOREIGN KEY (sampleID) REFERENCES asarLite.smpl (ID);
ALTER TABLE asarLite.md5_taxLite ADD CONSTRAINT FKmd5_taxLit27959 FOREIGN KEY (md5md5sum) REFERENCES asarLite.md5 (md5sum);
ALTER TABLE asarLite.md5_taxLite ADD CONSTRAINT FKmd5_taxLit822362 FOREIGN KEY (taxLiteID) REFERENCES asarLite.taxLite (ID);
ALTER TABLE asarLite.md5_KOrt ADD CONSTRAINT FKmd5_KOrt990010 FOREIGN KEY (md5md5sum) REFERENCES asarLite.md5 (md5sum);
ALTER TABLE asarLite.md5_KOrt ADD CONSTRAINT FKmd5_KOrt781708 FOREIGN KEY (KOrtid) REFERENCES asarLite.KOrt (id);
ALTER TABLE asarLite.abundance ADD CONSTRAINT FKabundance377377 FOREIGN KEY (md5md5sum) REFERENCES asarLite.md5 (md5sum);
ALTER TABLE asarLite.abundance ADD CONSTRAINT FKabundance495968 FOREIGN KEY (sampleID) REFERENCES asarLite.smpl (ID);
ALTER TABLE asarLite.seedLite_md5 ADD CONSTRAINT FKseedLite_m345617 FOREIGN KEY (seedLiteID) REFERENCES asarLite.seedLite (ID);
ALTER TABLE asarLite.seedLite_md5 ADD CONSTRAINT FKseedLite_m489908 FOREIGN KEY (md5md5sum) REFERENCES asarLite.md5 (md5sum);
