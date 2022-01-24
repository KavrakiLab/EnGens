use docking_results;
select *, a.name, extract_conformation(a.name) as conf, 
extract_fname(a.name) as fname,
extract_molecule(a.name) as mol
from adora2a_results as a
limit 100, 120;

select * from adora2a_results limit 41720600, 41720617;

select * from adora2a_results limit 1, 100;
select * from adora2a_frames limit 1, 100;

SHOW KEYS FROM adora2a_results WHERE Key_name = 'PRIMARY';

ALTER TABLE adora2a_results
ADD CONSTRAINT K_Name UNIQUE KEY (name);

ALTER TABLE adora2a_results
ADD COLUMN frame_num VARCHAR(50);

ALTER TABLE adora2a_results
ADD COLUMN n INT;
ALTER TABLE adora2a_results MODIFY COLUMN n INT not null auto_increment primary key;

SET SQL_SAFE_UPDATES = 0;

update docking_results.`adora2a_results`
SET docking_results.`adora2a_results`.frame_num = 
(
select docking_results.`adora2a_res`.frame_num from docking_results.`adora2a_res`
where docking_results.`adora2a_res`.rank_dock=docking_results.`adora2a_results`.index
);

ALTER TABLE adora2a_frames MODIFY frame_num VARCHAR(50);

ALTER TABLE adora2a_results ADD CONSTRAINT fk_fname_id 
FOREIGN KEY (frame_num) REFERENCES adora2a_frames(`frame_num`);

CREATE INDEX idconf ON adora2a_frames 
((extract_conformation(name)));

drop view adora2a_top05;

create table adora2a_top05 as
	SELECT
		A.index, A.frame_num, A.name, A.value
	FROM
			(
			Select ares.index, ares.frame_num, ares.name, ares.value
				,ROW_NUMBER() OVER(PARTITION BY ares.frame_num ORDER BY ares.index) AS RN
			from adora2a_results as ares
			) A
	WHERE A.RN < ceil(11743*0.005);

select *, extract_molecule2(name), locate(name, 'chembl') from adora2a_top05;

SET @@AUTOCOMMIT=0;
LOCK TABLES adora2a_frames WRITE, 
adora2a_results READ;
LOCK TABLES a2 WRITE;

INSERT INTO adora2a_frames
            select extract_conformation(adora2a_results.name) as frame_num, 
			extract_fname(adora2a_results.name) as fname,
			count(*) as dock_cnt
			from adora2a_results
			group by extract_conformation(adora2a_results.name);
UNLOCK TABLES;

select * from adora2a_frames;
delete from adora2a_frames where frame_num <1000000;

SET GLOBAL log_bin_trust_function_creators = 1;
