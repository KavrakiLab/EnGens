use docking_results;

create table adora2a_top1 as
	SELECT
		A.index, A.frame_num, A.name, A.value
	FROM
			(
			Select ares.index, ares.frame_num, ares.name, ares.value
				,ROW_NUMBER() OVER(PARTITION BY ares.frame_num ORDER BY ares.index) AS RN
			from adora2a_results as ares
			) A
	WHERE A.RN < ceil(11743*0.01);
    
    create table adora2a_top5 as
	SELECT
		A.index, A.frame_num, A.name, A.value
	FROM
			(
			Select ares.index, ares.frame_num, ares.name, ares.value
				,ROW_NUMBER() OVER(PARTITION BY ares.frame_num ORDER BY ares.index) AS RN
			from adora2a_results as ares
			) A
	WHERE A.RN < ceil(11743*0.05);
    
        create table adora2a_top10 as
	SELECT
		A.index, A.frame_num, A.name, A.value
	FROM
			(
			Select ares.index, ares.frame_num, ares.name, ares.value
				,ROW_NUMBER() OVER(PARTITION BY ares.frame_num ORDER BY ares.index) AS RN
			from adora2a_results as ares
			) A
	WHERE A.RN < ceil(11743*0.1);
    
drop table adora2a_top10_of_all;
create table adora2a_top10_of_all as
	Select ares.index, ares.frame_num, ares.name, ares.value
	from adora2a_results as ares order by ares.index 
	limit 0, 4172062;
    
select * from adora2a_top10_of_all limit 0, 1000;
select count(adora2a_results.index) from adora2a_results;
    