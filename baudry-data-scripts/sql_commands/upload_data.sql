use docking_results;
drop table adora2a_res;

create table adora2a_res (rank_dock INT NOT NULL auto_increment, primary key(rank_dock),
                              frame_num INT NOT NULL);

load data local
infile 'C:\\ProgramData\\MySQL\\MySQL Server 8.0\\Uploads\\adora2a_sr_003.txt'
into table adora2a_res
lines starting by 'frame'  terminated by '_' (frame_num);

select * from adora2a_res;

ALTER TABLE adora2a_res ADD CONSTRAINT fk_dock_id 
FOREIGN KEY (rank_dock) REFERENCES adora2a_results(`index`);

load data local
infile 'C:\\ProgramData\\MySQL\\MySQL Server 8.0\\Uploads\\adora2a_sr_003.txt'
into table adora2a_results
fields terminated by '\t' lines terminated by '\n' ();

drop table adora2a_frames;
CREATE TABLE adora2a_frames (frame_num INT NOT NULL,
			PRIMARY KEY (frame_num));
load data local
infile 'C:\\ProgramData\\MySQL\\MySQL Server 8.0\\Uploads\\adora2a_sr_003.txt'
into table adora2a_frames
lines starting by 'frame'  terminated by '_' (frame_num);

select count(*) from adora2a_frames;           

drop table adora2a_actives;
create table adora2a_actives(id MEDIUMINT NOT NULL AUTO_INCREMENT,
			mol_code VARCHAR(50) NOT NULL,
			PRIMARY KEY (id));

load data local
infile 'C:\\ProgramData\\MySQL\\MySQL Server 8.0\\Uploads\\ADORA2A_actives_final.mol2'
into table adora2a_actives
lines starting by '@<TRIPOS>MOLECULE\n' terminated by '\n'(mol_code);

select * from adora2a_actives;

drop table adora2a_magic;
CREATE TABLE adora2a_magic (ind INT NOT NULL,
			PRIMARY KEY (ind), name VARCHAR(50), flag INT);
            
load data local
infile 'C:\\ProgramData\\MySQL\\MySQL Server 8.0\\Uploads\\adora2a_magic.csv'
into table adora2a_magic
fields terminated by ','
lines terminated by '\r\n'
IGNORE 1 LINES;

drop view adora2a_magic_num;
create view adora2a_magic_num
as
select ind, name as magname, flag as flag, extract_fnum_magic(name) as f_num from adora2a_magic;




