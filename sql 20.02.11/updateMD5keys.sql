alter table md5_taxlite drop CONSTRAINT "md5_taxlite_md5md5sum_taxliteid_pkey";
alter table md5_taxlite drop CONSTRAINT "fkmd5_taxlit27959";
alter table md5_kort drop CONSTRAINT "md5_kort_md5md5sum_kortid_pkey";
alter table md5_kort drop CONSTRAINT "fkmd5_kort990010";
alter table seedlite_md5 drop CONSTRAINT "seedlite_md5_seedliteid_md5md5sum_pkey";
alter table seedlite_md5 drop CONSTRAINT "fkseedlite_m489908";
alter table abundance drop CONSTRAINT "abundance_md5md5sum_sampleid_pkey";
alter table abundance drop CONSTRAINT "fkabundance377377";
alter table md5 drop CONSTRAINT "md5_md5sum_pkey";

drop view dkres;
drop view funtaxab;
drop view funtaxall;
drop view smpl_stat;

--alter table drop CONSTRAINT "";

alter table md5 add column id serial;
alter table abundance add column md5id INTEGER;
alter table md5_taxlite add column md5id INTEGER;
alter table md5_kort add column md5id INTEGER;
alter table seedlite_md5 add column md5id INTEGER;

update abundance set md5id=(select id from md5 where md5.md5sum=abundance.md5md5sum);
update md5_taxlite set md5id=(select id from md5 where md5.md5sum=md5_taxlite.md5md5sum);
update md5_kort set md5id=(select id from md5 where md5.md5sum=md5_kort.md5md5sum);
update seedlite_md5 set md5id=(select id from md5 where md5.md5sum=seedlite_md5.md5md5sum);

--alter table md5 ADD CONSTRAINT "md5_md5sum_pkey" PRIMARY KEY ("id");
ALTER TABLE "asarlite"."abundance" ADD CONSTRAINT "abundance_md5id_sampleid_pkey" PRIMARY KEY ("md5id", "sampleid");
ALTER TABLE "asarlite"."md5_taxlite" ADD CONSTRAINT "md5_taxlite_md5id_taxliteid_pkey" PRIMARY KEY ("md5id", "taxliteid");
ALTER TABLE "asarlite"."md5_kort" ADD CONSTRAINT "md5_kort_md5id_kortid_pkey" PRIMARY KEY ("md5id", "kortid");
ALTER TABLE "asarlite"."seedlite_md5" ADD CONSTRAINT "seedlite_md5_md5id_seedliteid_pkey" PRIMARY KEY ("md5id", "seedliteid");

ALTER TABLE asarLite.abundance ADD CONSTRAINT FKabundance377377 FOREIGN KEY (md5id) REFERENCES asarLite.md5 (id);
ALTER TABLE asarLite.md5_taxLite ADD CONSTRAINT FKmd5_taxLit27959 FOREIGN KEY (md5id) REFERENCES asarLite.md5 (id);
ALTER TABLE asarLite.md5_KOrt ADD CONSTRAINT FKmd5_KOrt990010 FOREIGN KEY (md5id) REFERENCES asarLite.md5 (id);
ALTER TABLE asarLite.seedLite_md5 ADD CONSTRAINT FKseedLite_m489908 FOREIGN KEY (md5id) REFERENCES asarLite.md5 (id);

alter table abundance drop COLUMN md5md5sum;
alter table md5_taxlite drop COLUMN md5md5sum;
alter table md5_kort drop COLUMN md5md5sum;
alter table seedlite_md5 drop COLUMN md5md5sum;


create view dkres as
select distinct s.mgrastid as id,d.md5sum as md5,
       kort.description as annotation,
       kort.access as ko
   from smpl s join
        abundance a on s.id=a.sampleid join
        md5 d on a.md5id=d.id join
        md5_kort m on m.md5id=a.md5id join
        kort on kort.id=m.kortid;

create view funtaxab as
select distinct d.md5sum as md5sum,
 s.mgrastid,
 a.abundance,
 usp,
 sl.accession as ufun
 from abundance a join smpl s on a.sampleid=s.id join
        md5 d on a.md5id=d.id join
        seedlite_md5 sm on a.md5id=sm.md5id join
        md5_taxlite tm on a.md5id=tm.md5id join
        seedlite sl on sl.id=sm.seedliteid join
        taxlite tl on tl.id=tm.taxliteid;

create view funtaxall as
select distinct d.md5sum as md5sum,
       s.mgrastid,
       a.abundance,
       usp,species,genus,family,ordr,class,phylum,dom,
       fun4 as ufun, fun3 as fun2, fun2 as fun3, fun4 as fun1
  from abundance a join smpl s on a.sampleid=s.id join
        md5 d on a.md5id=d.id join
        seedlite_md5 sm on a.md5id=sm.md5id join
        md5_taxlite tm on a.md5id=tm.md5id join
       seedlite sl on sl.id=sm.seedliteid join
       taxlite tl on tl.id=tm.taxliteid;
create view smpl_stat as
  select t.name as stname,
         t.mgrastid as project,
         s.id,
         s.name,
         s.mgrastid as smpl,
         count(distinct m.md5id) as cnt,
         sum(abundance) as abundance
  from study t join
       smpl s on t.id=s.studyid left join
       abundance m on s.id=m.sampleid
  group by stname,project,s.id,s.name,smpl;
