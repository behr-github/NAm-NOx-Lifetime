-- the only way to actually write a "script" in SQL is to make it a "procedure"
-- which may actually be stored and run on the server. We need a delimiter to 
-- indicate when the procedure is done that is distinct from the normal delimiter.
delimiter ;;

-- this avoids an error about trying to update a table without using a primary
-- key in the WHERE clause
SET SQL_SAFE_UPDATES = 0;;

-- Set the value after the "use" command to the database that the MOVES
-- output went into
use moves_output_core_counties_2005to2014;;
drop procedure if exists aggregate_emis;;
create procedure aggregate_emis()
begin
drop table if exists TotalEmissions;
create table TotalEmissions( 
	moves_id int,
    moves_run_file varchar(2048),
	emis_year int, 
    emis_month int,
    species_id int,
    species_name varchar(128),
    total_emis_kg double
);

insert into totalemissions(moves_id, emis_year, emis_month, species_id, total_emis_kg)
	SELECT MOVESRunID, yearID, monthID, pollutantID, SUM(emissionQuant) FROM movesoutput GROUP BY MOVESRunID, monthID, pollutantID;
    
-- Not the best idea, but saves having to read the horrifically slow decoded output table
update totalemissions set species_name = 'NOx' WHERE species_id = 3;
update totalemissions set species_name = 'NO' WHERE species_id = 32;
update totalemissions set species_name = 'NO2' WHERE species_id = 33;

SET @run = (SELECT MIN(MOVESRunID) FROM movesoutput);
SET @runmax = (SELECT MAX(MOVESRunID) FROM movesoutput);
WHILE @run <= @runmax DO

	SET @movesfile = (SELECT runSpecFileName FROM movesrun WHERE MOVESRunID = @run);
	UPDATE totalemissions SET moves_run_file = @movesfile WHERE moves_id = @run;
	SET @run = @run + 1;

END WHILE;


END;;

call aggregate_emis();
