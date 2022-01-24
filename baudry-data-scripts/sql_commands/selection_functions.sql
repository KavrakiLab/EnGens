use docking_results;
DROP FUNCTION  IF EXISTS extract_conformation;
DROP FUNCTION  IF EXISTS extract_fname;
DROP FUNCTION  IF EXISTS extract_molecule;

DROP FUNCTION IF EXISTS extract_molecule_withoutline;
DROP FUNCTION IF EXISTS extract_fnum_magic;

DELIMITER $$
CREATE FUNCTION extract_fnum_magic(s CHAR(50))
RETURNS INT
BEGIN
    RETURN substr(s,  locate("_", s)+1, length(s)-locate("_", s));
END$$
DELIMITER ;


DELIMITER $$
CREATE FUNCTION extract_molecule_withoutline(s CHAR(50)) 
RETURNS VARCHAR(50)
BEGIN
	DECLARE index1 INT;
	DECLARE index2 INT;
    DECLARE str_tmp VARCHAR(50);
    DECLARE res VARCHAR(50);
    
    SET index1 = locate("_", s);
    SET index2 = locate(".", s);
    SET str_tmp = substr(s, index1+1, index2-index1-1);
    
    IF locate("_", str_tmp) > 0 THEN
		SET res = substr(str_tmp, 1, locate("_", str_tmp)-1);
	ELSE
		SET res = str_tmp;
    END IF;
    
    return res;
END$$


CREATE FUNCTION extract_conformation(s CHAR(50)) 
RETURNS integer
BEGIN
	DECLARE tmpstr VARCHAR(50);
	DECLARE res integer;
	DECLARE prefix VARCHAR(50) DEFAULT "at-frame";
	DECLARE suffix VARCHAR(50) DEFAULT ".pdbqt";
    
    IF instr(s, prefix) THEN
		SET tmpstr = substr(s,  
        length(prefix)+1, 
        locate("_", s)-length(prefix)-1);
	ELSEIF instr(s, "frame") THEN
		SET prefix = "frame";
		SET tmpstr = substr(s,  
        length(prefix)+1, 
        locate("_", s)-length(prefix)-1);
	END IF;
    SET res = cast(tmpstr as unsigned);
    RETURN res;
END$$

CREATE FUNCTION extract_fname(s CHAR(50)) 
RETURNS VARCHAR(50)
BEGIN
    RETURN substr(s,  1, locate("_", s)-1);
END$$

DELIMITER ;
