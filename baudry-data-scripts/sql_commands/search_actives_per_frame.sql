use docking_results;

show tables;

-- number of all molecules in 10% for each frame
select *, count(adora2a_top05.index) as cnt  from adora2a_top05
group by adora2a_top05.frame_num;

-- restrain top 10% for each frame to unique actives 
drop table if exists unique_actives_per_frame05;

create table unique_actives_per_frame02
as
select * from 
adora2a_top05
left join ( select * from adora2a_actives group by mol_code) as adora2a_actives
 on 
-- join on molecule code (without "_2/_3/_4")
extract_molecule_withoutline(adora2a_top05.name)=adora2a_actives.mol_code
WHERE mol_code IS NOT NULL;

-- unique actives per frame
select *, min(value) as minval, count(*) as active_cnt from unique_actives_per_frame1
group by unique_actives_per_frame1.frame_num  order by minval;

-- significant frames sign. threshold 107
select * from (
select *, min(value) as minval, count(*) as active_cnt from unique_actives_per_frame
group by unique_actives_per_frame.frame_num  order by minval ) as tmp
where  active_cnt >= 107 and frame_num <= 600000;

-- put the significant frames in a table
create table significant_frames_top10
as 
select tmp.frame_num from (
select *, min(value) as minval, count(*) as active_cnt from unique_actives_per_frame
group by unique_actives_per_frame.frame_num  order by minval ) as tmp
where  active_cnt >= 107 and frame_num <= 600000;

select * from significant_frames_top10 order by cast(frame_num as unsigned);

-- combine with magic
create view magic_not_sig 
as
select * from 
(select * from 
(
select *, min(value) as minval, count(*) as active_cnt from unique_actives_per_frame
group by unique_actives_per_frame.frame_num  order by minval ) as tmp
right join adora2a_magic_num as magic on magic.f_num = tmp.frame_num) as tmp2
where tmp2.flag=1 and active_cnt<107 and frame_num <= 600000;

select * from magic_not_sig;

select *, extract_molecule_withline(a.name) from adora2a_top10 as a limit 1, 1000
