use docking_results;

select * from (select *, count(adora2a_top10.index) as cnt  from adora2a_top10 
group by adora2a_top10.frame_num) as tmp where not cnt=1174;

select count(*) from adora2a_top10_of_all; -- 4172062 entries are top 10% of all docking scores

-- unique actives in top 0.5
create view unique_actives_top05
as 
select * from 
(-- take top 0.5% and extract mol_codes
select *, 
extract_molecule_withline(adora2a_top10_of_all.name) as mol_full, 
extract_molecule_withoutline(adora2a_top10_of_all.name) as mol_code_part
from adora2a_top10_of_all 
 order by adora2a_top10_of_all.index limit 0, 208603
 ) -- 208606 is top 0.5% of all docks with mol codes 
as moltop05 
-- join on mol code (without _2 ) from all actives
inner join adora2a_actives on moltop05.mol_code_part=adora2a_actives.mol_code
-- group by mol code (with _2) to get unique actives
group by moltop05.mol_full
;
-- select from view --
select * from unique_actives_top05;
-- res: 654 with full, 373 excluding full

-- identify significant frames here 
-- get only significant actives
select * from (
-- count actives per frame
select *, count(*) as active_count from
-- get all actives from top 0.5%
(select * from unique_actives_top05) as actives_top05
group by actives_top05.frame_num
) as tmp where tmp.active_count >9;

-- try to figure out the MAGIC!
-- in top 10 for each complex
-- how many actives are there?

-- take only frames that have significant active count (more than 9)
select * from (
-- join with actives to extract only active counts
-- group by frame to count the molecules (full) within 
select *, count(moltop05.mol_full) as actives_cnt from 
-- select top 0.5 from top 10
-- extract the molecule (full and the code only)
(select *, 
extract_molecule_withline(adora2a_top10_of_all.name) as mol_full, 
extract_molecule_withoutline(adora2a_top10_of_all.name) as mol_code_ext
from adora2a_top10_of_all 
 order by adora2a_top10_of_all.index limit 0, 208603) -- 208606 is top 0.5% of all docks
as moltop05 
inner join adora2a_actives on moltop05.mol_code_ext=adora2a_actives.mol_code 
group by  moltop05.frame_num ) 
as tmp where tmp.actives_cnt > 9;

select *, extract_molecule_withline(tmp.name), extract_molecule_withoutline(tmp.name)
as mol, locate('_2',tmp.name) from (select * from adora2a_results limit 0, 31091) as tmp
where locate('_2', tmp.name)>0 
limit 0, 1000;

-- find significant frames?
-- they are significant if in the top 0.05% [...] we have
-- more that random actives 9...
 

select * from adora2a_actives;
