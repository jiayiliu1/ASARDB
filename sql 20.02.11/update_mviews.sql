 insert into asarLite.mv_abund_stat(sampleid,cnt,abundance)
 select sampleid,
 count(distinct m.md5id) as cnt,
 sum(abundance) as abundance
 from abundance m
 where sampleid in (
    select id as sampleid from smpl except
    select sampleid from mv_abund_stat)
 group by sampleid



insert into mv_smpl_stat (stname,project,id,name,smpl,cnt,abundance)
select distinct t.name as stname,
 t.mgrastid as project,
 s.id,
 s.name,
 s.mgrastid as smpl,
 m.cnt,
 m.abundance
from study t join
 smpl s on t.id=s.studyid join
 mv_abund_stat m on s.id=m.sampleid
 order by s.id;

insert into mv_dkres (id,md5,annotation,ko) 
select id,md5,annotation,ko 
from dkres 
where id in (
   select distinct s.mgrastid as id from smpl s 
   except 
   select distinct id from mv_dkres);

insert into mv_funtaxab (md5sum,mgrastid,abundance,usp,ufun) 
select md5sum,mgrastid,abundance,usp,ufun 
from funtaxab
where mgrastid in (
   select distinct s.mgrastid as mgrastid from smpl s 
   except 
   select distinct mgrastid from mv_funtaxab);