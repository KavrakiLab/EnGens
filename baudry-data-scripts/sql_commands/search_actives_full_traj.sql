use docking_results;

-- table with top10 adora2a_top10
select count(*) from adora2a_top10;
select count(*) from (select * from adora2a_top10 group by frame_num) as tmp;

-- view  with significant frames of top10 
create view top10_sig as
select top10.index, top10.frame_num, name, value from adora2a_top10 as top10
left join significant_frames_top10 as sig on sig.frame_num = top10.frame_num
WHERE sig.frame_num IS NOT NULL;

select * from top10_sig;

-- actives number in top10_sig
create table unique_actives_sig_top10
as
select * from 
top10_sig
left join ( select * from adora2a_actives group by mol_code) as adora2a_actives
 on 
-- join on molecule code (without "_2/_3/_4")
extract_molecule_withoutline(top10_sig.name)=adora2a_actives.mol_code
WHERE mol_code IS NOT NULL;

-- count of actives within significant frames
select count(*) from(
select count(*) from unique_actives_sig_top10 group by extract_molecule_withline(name)) as tmp;

-- count of actives without the significant frames

-- extract the actives

-- actives number in top10_sig
create table unique_actives_top10
as
select * from 
adora2a_top10 as top10
left join ( select * from adora2a_actives group by mol_code) as adora2a_actives
 on 
-- join on molecule code (without "_2/_3/_4")
extract_molecule_withoutline(top10.name)=adora2a_actives.mol_code
WHERE mol_code IS NOT NULL;

select * from (
select *, count(*) as active_cnt from unique_actives_top10 group by frame_num order by active_cnt desc
) as tmp
where tmp.frame_num = 226400;

-- conbine with magic

select * from unique_actives_top10 as top10
left join adora2a_magic_num as mag on top10.frame_num = mag.f_num
where mag.flag = 1;

